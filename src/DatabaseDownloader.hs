{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module DatabaseDownloader (downloadDatabase, downloadSummary, PDBRecord (..)) where

import BioTypes (Fasta, FastaRecord (fAliases, fSequence), ResidueSequence, parseFasta)
import Control.Applicative (liftA2, optional)
import Control.Arrow (first)
import Control.Lens ((<&>))
import Data.Aeson (FromJSON (parseJSON), ToJSON (toJSON), eitherDecode, genericParseJSON, genericToJSON, object, (.=))
import qualified Data.Aeson as A
import Data.Aeson.Encode.Pretty
import Data.Bifunctor (bimap)
import qualified Data.ByteString.Lazy as L
import Data.ByteString.Lazy.Char8 (unpack)
import Data.Char (isDigit)
import Data.Csv hiding (encode, (.=))
import Data.Data (Data)
import Data.Function (on)
import Data.List (find, groupBy, intercalate)
import Data.Maybe (fromJust, fromMaybe)
import qualified Data.Vector as V
import Debug.Trace (traceShowId)
import GHC.Generics (Generic)
import Network.HTTP.Conduit (simpleHttp)
import System.Directory (createDirectoryIfMissing)
import System.FilePath ((</>))

downloadDatabase :: FilePath -> (PDBRecord -> Bool) -> L.ByteString -> IO ()
downloadDatabase path predicate summaryFile = do
  records <- readTsvRecords summaryFile
  let filtered = filter predicate records
      ids = map pdb filtered
  putStrLn $ "Selected " ++ show (length ids) ++ " PDB entries based on predicate."
  putStrLn "Downloading PDB files..."
  downloadPDBs path ids
  putStrLn "Downloading FASTA files..."
  downloadFastas path ids
  putStrLn "Writing initial JSON descriptions..."
  writeJSONs path $ map mergeChains $ groupBy ((==) `on` pdb) filtered
  putStrLn "Downloading chain numberings..."
  downloadNumberings path ids
  putStrLn "Download complete."

summaryURL :: String
summaryURL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"

downloadSummary :: FilePath -> IO ()
downloadSummary file = simpleHttp summaryURL >>= L.writeFile file

data PDBRecord = Record
  { pdb :: String
  , resolution :: Double
  , method :: String
  , scfv :: Bool
  , heavySpecies :: String
  , hChain :: String
  , lChain :: String
  }
  deriving (Generic, Data, Show)

instance FromJSON PDBRecord -- TODO: WTF is that
instance ToJSON PDBRecord

instance FromNamedRecord PDBRecord where
  parseNamedRecord r =
    Record
      <$> r .: "pdb"
      <*> (optional (r .: "resolution") <&> fromMaybe 0.0)
      <*> r .: "method"
      <*> (r .: "scfv" <&> read)
      <*> r .: "heavy_species"
      <*> r .: "Hchain"
      <*> r .: "Lchain"

-- instance ToJSON PDBRecord where
--   toJSON record =
--     object
--       [ "pdb" .= pdb record
--       , "resolution" .= resolution record
--       , "method" .= method record
--       , "scfv" .= scfv record
--       , "heavy_species" .= heavySpecies record
--       , "Hchain" .= hChain record
--       , "Lchain" .= lChain record
--       ]

-- TODO: Download in bulk
downloadPDBs :: FilePath -> [String] -> IO ()
downloadPDBs path = mapM_ downloadWithDir
 where
  url id = "https://files.rcsb.org/download/" ++ id ++ ".pdb"
  filePath id = path </> id </> (id ++ ".pdb")
  downloadWithDir id = do
    createDirectoryIfMissing True (path </> id)
    simpleHttp (url id) >>= L.writeFile (filePath id)

downloadFastas :: FilePath -> [String] -> IO ()
downloadFastas path = mapM_ downloadWithDir
 where
  url id = "https://www.rcsb.org/fasta/entry/" ++ id ++ "/download"
  filePath id = path </> id </> (id ++ ".fasta")
  downloadWithDir id = do
    createDirectoryIfMissing True (path </> id)
    simpleHttp (url id) >>= L.writeFile (filePath id)

fetchNumbering :: ResidueSequence -> IO L.ByteString
fetchNumbering s = simpleHttp $ "http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq=" ++ s ++ "&scheme=-m" -- Martin

data NumberingJSON = NumberingJSON
  { heavy :: String -- TODO: make something more sensible than string
  , light :: String
  }
  deriving (Generic, Data, Show)

instance FromJSON NumberingJSON where
  parseJSON = genericParseJSON A.defaultOptions
instance ToJSON NumberingJSON where
  toJSON = genericToJSON A.defaultOptions

downloadNumberings :: FilePath -> [String] -> IO ()
downloadNumberings = mapM_ . downloadNumbering

downloadNumbering :: FilePath -> String -> IO ()
downloadNumbering path id = getFasta >>= makeNumbering >>= L.writeFile numberingFile
 where
  numberingFile = path </> id </> "numbering.txt"
  jsonFile = path </> id </> "description.json"
  getFasta = parseFasta <$> readFile (path </> id </> (id ++ ".fasta"))
  getChainNames file = case eitherDecode file :: Either String PDBRecord of
    Left err -> error $ "Error parsing JSON while downloading numbering for " ++ id ++ ": " ++ err
    Right record -> (head $ hChain record, head $ lChain record)
  makeNumbering :: Fasta -> IO L.ByteString
  makeNumbering fasta = L.readFile jsonFile >>= (uncurry (liftA2 (<>)) . bimap getNum getNum) . getChainNames
   where
    getNum c = fetchNumbering $ fSequence $ fromJust $ find (elem c . fAliases) fasta

writeJSONs :: FilePath -> [PDBRecord] -> IO ()
writeJSONs path = mapM_ $ \record -> L.writeFile (path </> pdb record </> "description.json") (encodePretty record)

readTsvRecords :: L.ByteString -> IO [PDBRecord]
readTsvRecords tsvFile = case decodeByNameWith options tsvFile of
  Left err -> error $ "Error parsing TSV summary file '" ++ err
  Right (_, records) -> return $ V.toList records
 where
  options = defaultDecodeOptions{decDelimiter = fromIntegral (fromEnum '\t')}

mergeChains :: [PDBRecord] -> PDBRecord
mergeChains [] = error "Cannot merge an empty list of records"
mergeChains (r : rs) =
  Record
    { pdb = pdb r
    , resolution = resolution r
    , method = method r
    , scfv = scfv r
    , heavySpecies = heavySpecies r
    , hChain = intercalate "," $ map hChain (r : rs)
    , lChain = intercalate "," $ map lChain (r : rs)
    }
