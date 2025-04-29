{-# LANGUAGE LambdaCase #-}

module Main (main) where

import qualified Data.ByteString.Lazy as L
import qualified Data.ByteString.Lazy.Char8 as LC8
import DatabaseDownloader (PDBRecord (..), downloadDatabase, downloadSummary)
import DatabaseProcessor (processDatabase)
import System.Directory.Internal.Prelude (getArgs)

usage :: String
usage =
  unlines
    [ "Usage: myprogram <command> [options]"
    , ""
    , "Commands:"
    , "  database <subcommand>       Manage the database."
    , "    subcommands:"
    , "      download-summary        Download database summary(required)"
    , "      download <path> [n]     Download [first n] records from the database."
    , "      process  <path>         Process the database."
    , ""
    , "  aminoacids                  Manage amino acids (not implemented yet)."
    ]

summaryFile :: String
summaryFile = "sabdab-summary.tsv"

isGoodPDB :: PDBRecord -> Bool
isGoodPDB pdb_record =
  heavySpecies pdb_record == "homo sapiens"
    && resolution pdb_record <= 3.0
    && (method pdb_record `elem` ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"])
    && not (scfv pdb_record)

handleDatabase :: [String] -> IO ()
handleDatabase = \case
  ["download-summary"] -> downloadSummary summaryFile
  ["download", path] -> L.readFile summaryFile >>= downloadDatabase path isGoodPDB
  ["download", path, n] -> L.readFile summaryFile >>= downloadDatabase path isGoodPDB . LC8.unlines . take (read n) . LC8.lines
  ["process", path] -> processDatabase path
  _ -> putStrLn usage

main :: IO ()
main =
  getArgs >>= \case
    "database" : a -> handleDatabase a
    "aminoacids" : _ -> putStrLn "Not implemented yet"
    _ -> putStrLn usage
