-- TODO:
-- 29.04.25 - сделать разметку антител и обрезку по концу J фрагментов (они же концы вариабельных доменов)
-- What does it mean?
-- 06.05.25 - реализовать парсинг .fasta, шлифовка кода и тестирование работоспособности (аналитика)
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module DatabaseProcessor (processDatabase) where

import BioTypes (Atom (..), Distance, Residue (..), ResidueName, distance, parseAtom)
import Control.Applicative (liftA2)
import Control.Monad (join, zipWithM_)
import Data.Aeson (KeyValue ((.=)), ToJSON (toEncoding), defaultOptions, eitherDecode, encode, genericToEncoding, object)
import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as LC8
import Data.Char (ord)
import Data.Function (on)
import qualified Data.HashMap.Lazy as Map
import Data.List (find, groupBy, nub, (\\))
import Data.List.Split (splitOn)
import Data.Maybe (catMaybes, fromJust)
import qualified DatabaseDownloader as DD
import GHC.Generics
import System.Directory (listDirectory)
import System.FilePath ((</>))

processDatabase :: FilePath -> IO () -- Main function here
processDatabase dir = listDirectory dir >>= mapM_ (processPDBDir dir)

toShortName :: L.ByteString -> L.ByteString
toShortName = id -- TODO: Make this

parseResidueTable :: L.ByteString -> ResidueName -> [MissingAtom]
parseResidueTable _contents _name = [] -- TODO: Make it...

-- parseResidueTable = mapLookup . map parseResidue . groupBy (groupResidue "") . filter informational . LC8.lines
--  where
--   informational s = not (";" `LC8.isPrefixOf` s) && not (LC8.null s)
--   groupResidue indent a b = (LC8.pack (indent ++ "[") `LC8.isPrefixOf` a) && (LC8.pack (indent ++ " ") `LC8.isPrefixOf` b)
--   mapLookup = (.) fromJust . (flip Map.lookup . Map.fromList) -- TODO: understand this
--   parseResidue (n : d) = (toShortName $ LC8.words n !! 1, getParseAtoms d)
--   parseResidue _ = error "Invalid residue description"
--   parseAtoms (_ : a) = map ((!! 0) . (!! 1) . LC8.words) a -- TODO : !! 0 wont work for LC8
--   getParseAtoms = parseAtoms . tail . findAtoms . groupBy (groupResidue " ")
--
processPDBDir :: FilePath -> String -> IO ()
processPDBDir dir name = processPDB . parseResidueTable <$> L.readFile "aminoacids.rtp" <*> readFiles >>= writeFiles
 where
  filePaths = map ((dir </> name) </>) ["description.json", name ++ ".pdb", name ++ ".fasta", "numbering"]
  readFiles = mapM L.readFile filePaths
  writeFiles = zipWithM_ L.writeFile filePaths

-- TODO: Add numbering file here as well
processPDB :: (ResidueName -> [MissingAtom]) -> [L.ByteString] -> [L.ByteString] -- See filePaths above
processPDB residueTable [jsonContent, pdbContent, fastaContent, numberingContent] = case eitherDecode jsonContent of
  (Left err) -> [makeReason $ "Error parsing JSON file: " ++ err]
  (Right initialRecord) -> chooseBestReflection description newPDB fastaContent
   where
    description =
      Description
        { resolution = DD.resolution initialRecord
        , method = DD.method initialRecord
        , hChains = splitOn "," $ DD.hChain initialRecord
        , lChains = splitOn "," $ DD.lChain initialRecord
        , numbering = [] -- parseNumbering numberingContent TODO: Do
        , reflections = [] -- catMaybes $ zipWith analyzeReflection (hChains description) (lChains description) -- TODO: do
        }
    newPDB = cleanupPDB description pdbContent
    pdbAtoms = map (parseAtom . LC8.unpack) $ init $ LC8.lines newPDB -- init for the END
    analyzeReflection :: L.ByteString -> L.ByteString -> Maybe Protein
    analyzeReflection h l =
      ( \residues ->
          Protein
            { missingResidues = getMissingResidues residues foundResidues
            , missingAtoms = getMissingAtoms residueTable foundResidues pdbAtoms
            , atomClashes = getAtomClashes pdbAtoms
            }
      )
        <$> getChainFromFasta h fastaContent -- ???
     where
      foundResidues = pdbResidues $ filter (isChain $ head (LC8.unpack h)) pdbAtoms -- TODO: Pasted head here mindlessly...
      isChain chain atom = chain == residueChain (atomResidue atom)
processPDB _ _ = [makeReason "Missing files"]

makeReason :: String -> L.ByteString
makeReason bs = encode $ object ["reason" .= bs]

-- Makes a file set for the best reflection, or a json with a missing reason otherwise
chooseBestReflection :: PDBDescription -> L.ByteString -> L.ByteString -> [L.ByteString]
chooseBestReflection _ _ _ = [makeReason "TODO: Choose best reflection"]

data Protein = Protein
  { missingResidues :: [(Int, Int)]
  , missingAtoms :: [MissingAtom]
  , atomClashes :: [(Int, Int, Distance)]
  }
  deriving (Generic)
data PDBDescription = Description
  { resolution :: Double
  , method :: String
  , hChains :: [String]
  , lChains :: [String]
  , numbering :: [(Int, String)]
  , reflections :: [Protein]
  }
  deriving (Generic)

instance ToJSON MissingAtom where
  toEncoding = genericToEncoding defaultOptions

instance ToJSON Protein where
  toEncoding = genericToEncoding defaultOptions

-- How will this actually encode...
instance ToJSON PDBDescription where
  toEncoding = genericToEncoding defaultOptions

pdbResidues :: [Atom] -> [Residue]
pdbResidues = nub . map atomResidue -- n^2

-- TODO: In DatabaseDownloader fastas should have adequate names beforehand...
-- Assuming header looks like "> X ..."
getChainFromFasta :: L.ByteString -> L.ByteString -> Maybe [Residue]
getChainFromFasta chainName = fmap toResidues . find isRightChain . map flattenRecord . groupBy fastaRecord . filter (not . L.null) . LC8.lines
 where
  fastaRecord a b = L.head a == fromIntegral (ord '>') && L.head b /= fromIntegral (ord '>') -- I'm sorry
  flattenRecord (h : rs) = (h, L.concat rs) -- Shouldn't be necessary
  flattenRecord x = error $ "Invalid FASTA:\n" ++ unlines (map LC8.unpack x)
  isRightChain (h, _) = LC8.words h !! 1 == chainName
  toResidue index name = Residue{residueId = index, residueChain = head $ LC8.unpack chainName, residueName = name}
  toResidues = zipWith toResidue [0 ..] . map (: []) . LC8.unpack . snd -- TODO: This should be index 0? Or what? See docs

-- Get ranges of missing residues
getMissingResidues :: [Residue] -> [Residue] -> [(Int, Int)]
getMissingResidues = (groupConseq .) . ((\\) `on` map residueId)

-- Group consecutive elements
groupConseq :: [Int] -> [(Int, Int)]
groupConseq = foldr makeRange []
 where
  makeRange a ((l, h) : xs) | l - a == 1 = (a, h) : xs
  makeRange a xs = (a, a) : xs -- TODO: Check that in the morning...

-- TODO: Make a better name.
data MissingAtom = MissingAtom
  { missingAtomResidue :: String
  , missingAtomResidueId :: Int
  , missingAtomType :: Char
  }
  deriving (Generic, Eq)

getMissingAtoms :: (ResidueName -> [MissingAtom]) -> [Residue] -> [Atom] -> [MissingAtom]
getMissingAtoms getAtoms residues atoms = concatMap (getAtoms . residueName) residues \\ map toMissing atoms

toMissing :: Atom -> MissingAtom
toMissing a = MissingAtom (residueName $ atomResidue a) (residueId $ atomResidue a) (atomType a)

-- TODO: Add complications...
-- only count is actualy needed
getAtomClashes :: [Atom] -> [(Int, Int, Distance)]
getAtomClashes = filter tooClose . filter different . map idsDist . join (liftA2 (,))
 where
  tooClose (_, _, dist) = dist < 2.5
  idsDist (a, b) = (atomId a, atomId b, distance a b)
  different (a, b, _) = a /= b

-- TODO: keep non-CDRs from right chains only
-- use numbering description
cleanupPDB :: PDBDescription -> L.ByteString -> L.ByteString
cleanupPDB description = LC8.unlines . filter isAtomOrEnd . LC8.lines
 where
  isAtomOrEnd line = L.take 4 line == "ATOM" || L.take 3 line == "END"
