module BioTypes where

import Data.List (groupBy)
import Data.List.Split (splitOn)
import Data.Ord (comparing)
import Debug.Trace (traceId, traceShowId)
import Text.Regex.TDFA

type Position = (Double, Double, Double)
type Distance = Double
type ResidueName = String -- TODO: make enum

data Residue = Residue
  { residueName :: ResidueName
  , residueChain :: Char
  , residueId :: Int
  }
  deriving (Show, Eq)

data Atom = Atom
  { atomId :: Int
  , atomName :: String
  , atomPos :: Position
  , atomType :: Char -- TODO: Maybe use enum
  , atomRadius :: Double
  , atomResidue :: Residue
  }
  deriving (Show, Eq)

instance Ord Atom where
  compare = comparing atomId

type ResidueSequence = String
data FastaRecord = FastaRecord
  { fAliases :: [Char]
  , fSequence :: ResidueSequence
  }
  deriving (Show)
type Fasta = [FastaRecord]

parseFasta :: String -> Fasta
parseFasta = map parseRecord . groupBy sameRecord . lines
 where
  sameRecord a b = head a == '>' && head b /= '>'
  parseRecord (header : residueSequence) = FastaRecord{fAliases = aliases, fSequence = concat residueSequence}
   where
    aliases = map handleAlias . splitOn ", " $ extractAliases header
    extractAliases :: String -> String
    extractAliases input = case input =~ "(Chain|Chains) ([^|]+)" of
      (_ : _ : gem : _) : _ -> gem
      e -> error $ "Invalid FASTA header: " ++ header ++ "\nfound: " ++ show e
    handleAlias :: String -> Char
    handleAlias [a] = a
    handleAlias [_, '[', 'a', 'u', 't', 'h', ' ', b, ']'] = b
    handleAlias _ = error $ "Invalid FASTA header: " ++ header
  parseRecord _ = error "Invalid FASTA:"

radiusMap :: Char -> Double
radiusMap 'C' = 1.70
radiusMap 'H' = 1.10
radiusMap 'N' = 1.55
radiusMap 'O' = 1.52
radiusMap 'P' = 1.80
radiusMap 'S' = 1.80
radiusMap _ = 1.5

distance :: Atom -> Atom -> Distance
distance a1 a2 = sqrt $ (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
 where
  (x1, y1, z1) = atomPos a1
  (x2, y2, z2) = atomPos a2

parseAtom :: String -> Atom
parseAtom s =
  Atom
    { atomId = read (take 5 $ drop 6 s)
    , atomName = filter (/= ' ') (take 4 $ drop 12 s)
    , atomPos = (read (take 8 $ drop 30 s), read (take 8 $ drop 38 s), read (take 8 $ drop 46 s))
    , atomType = s !! 77
    , atomRadius = radiusMap (s !! 77)
    , atomResidue =
        Residue
          { residueName = filter (/= ' ') (take 3 $ drop 17 s)
          , residueChain = s !! 21
          , residueId = read (take 4 $ drop 22 s)
          }
    }
