{-# LANGUAGE DeriveDataTypeable #-}
{-# LANGUAGE DeriveGeneric #-}
{-# LANGUAGE OverloadedStrings #-}
{-# OPTIONS_GHC -Wno-name-shadowing #-}

module DatabaseDownloader (downloadDatabase, downloadSummary, PDBRecord (..)) where

import BioTypes (ResidueSequence)
import Control.Applicative (optional)
import Control.Lens ((<&>))
import Data.Aeson (FromJSON, ToJSON (toJSON), encode, object, (.=))
import Data.Aeson.Encode.Pretty
import qualified Data.ByteString.Lazy as L
import Data.Csv hiding (encode, (.=))
import Data.Data (Data)
import Data.Function (on)
import Data.List (groupBy, intercalate)
import Data.Maybe (fromMaybe)
import qualified Data.Vector as V
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

instance ToJSON PDBRecord where
  toJSON record =
    object
      [ "pdb" .= pdb record
      , "resolution" .= resolution record
      , "method" .= method record
      , "scfv" .= scfv record
      , "heavy_species" .= heavySpecies record
      , "Hchain" .= hChain record
      , "Lchain" .= lChain record
      ]

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

downloadNumbering :: ResidueSequence -> IO L.ByteString
downloadNumbering s = simpleHttp $ "http://www.bioinf.org.uk/abs/abnum/abnum.cgi?plain=1&aaseq=" ++ s ++ "&scheme=Martin"

-- TODO: Just what is numbering for a pdb??
downloadNumberings :: FilePath -> [String] -> IO ()
downloadNumberings path = mapM_ $ \id -> getFasta id >>= downloadNumbering >>= L.writeFile (filePath id)
 where
  filePath id = path </> id </> (id ++ "numbering")
  getFasta id = readFile $ path </> id </> (id ++ ".fasta") -- TODO: Get residueSequence here.

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
