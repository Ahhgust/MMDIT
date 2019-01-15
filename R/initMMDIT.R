#!/usr/bin/env Rscript
# Written by August Woerner
# This creates the (empty) database
# makeDB function
#
# And it loads it with the appropriate data

# database dependencies
library(RSQLite)
library(sqldf)
library(XLConnect)

# quality of life dependencies
suppressPackageStartupMessages( library(tidyverse) )

# for reading fasta files (reference genome)
library(seqinr)

# CONSTANTS
# data directory
dataDir <- 'data'

# final database filename
dbName <- 'mmdit.sqlite3'

# reference genome
refGenome <- 'rCRS.fasta'

precisionIDBed <- "PrecisionID_mtDNA_WG_targets.bed"

sqliteURL <- 'https://www.sqlite.org/src/tarball/sqlite.zip.gz?r=release'

editdistFile <- 'spellfix.c'


# helper function. 
# Normally we send a query and then clear it. this wraps it up into 1 statement
dbSendQueryAndClear <- function(db, query) {
  dbClearResult( dbSendQuery(db, query) )
}

# This function takes in a path (e.g., data/mmdit.sqlite3)
# and creates a database with the requisite tables in it
# if overwrite is FALSE, then it won't overwrite an existing database
# in that path.
# otherwise, it will clobber what was there (use overwrite=TRUE with caution!!)
makeDB <- function(path, dbFilename, overwrite=FALSE) {
  
  dbFile <- paste(path, dbFilename, sep="/")
  
  if (file.exists(dbFile)) { #DB already exists
    if (overwrite==FALSE) {
      stop(
          paste("DB file " , dbFile, " already exists. Let's not overwrite things, shall we?")
          )
      return(NULL)
    }   
  }

  db <- dbConnect( RSQLite::SQLite(), dbFile, loadable.extensions = TRUE)
  
  # create the database tables:
  
  # information on the amplicons
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS amps") 
  dbSendQueryAndClear(db,
              "CREATE TABLE amps 
              (ampid INTEGER PRIMARY KEY ASC, 
              kit TEXT NOT NULL, 
              start INTEGER NOT NULL, 
              stop INTEGER NOT NULL)" ) 

  
  # and the sample IDs (and their corresponding populations)
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS populations") 
  dbSendQueryAndClear(db,
              "CREATE TABLE populations 
              (sampleid TEXT PRIMARY KEY, 
              pop TEXT NOT NULL)" )
   
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS haplotypes")
  dbSendQueryAndClear(db,
              "CREATE TABLE haplotypes
              (hapid INTEGER PRIMARY KEY ASC, 
               ampid INTEGER NOT NULL,
               sampleid TEXT NOT NULL,
               haplotype TEXT NOT NULL, 
               haplotype_count INTEGER DEFAULT 1,
                 FOREIGN KEY (ampid) REFERENCES amps(ampid),
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid))" )

  #TODO: THIS IS WRONG!
  # information on the amplicon dropout
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS amp_dropout")
  dbSendQueryAndClear(db,
              "CREATE TABLE amp_dropout 
              (ampid INTEGER NOT NULL,
              dropout_prob REAL NOT NULL, 
              sampleid TEXT NOT NULL,
                 FOREIGN KEY (ampid) REFERENCES amps(ampid),
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                 PRIMARY KEY (ampid, sampleid)
              )" )
 
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS haplogroups")
  dbSendQueryAndClear(db,
              "CREATE TABLE haplogroups 
              (sampleid TEXT NOT NULL,
                 haplogroup TEXT NOT NULL,
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                 PRIMARY KEY (sampleid)
              )" )
  
  # the individuals that we have full mitochondrial genomes for
  # provenance must describe where the data came from
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS fullMtgenomes") 
  dbSendQueryAndClear(db, 
                      "CREATE TABLE fullMtgenomes
                      (sampleid TEXT NOT NULL,
                       provenance TEXT NOT NULL,
                      FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                      PRIMARY KEY (sampleid)
                      )" )
  
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS fullMtgenomesDiffseqs") 
  dbSendQueryAndClear(db, 
            "CREATE TABLE fullMtgenomesDiffseqs
              (sampleid TEXT NOT NULL,
                 diffPosition INTEGER NOT NULL,  
                 diffEvent TEXT NOT NULL,
                 diffBasecall TEXT NOT NULL,
                 rowid INTEGER NOT NULL,
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                 FOREIGN KEY (sampleid) REFERENCES fullMtgenomes(sampleid),
                 PRIMARY KEY (rowid)
              )" )
  
   

  return(db)
}

# call this first.
# this makes a table of just the reference sequence
# this is hard-coded to the rCRS sequence
# though this can be undone (doing so would require completely redoing the database)
loadReferenceGenome <- function(db, fastaFile) {
  # a trivial table of just the reference genome
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS mtgenome")
  dbSendQueryAndClear(db,
              "CREATE TABLE mtgenome
              (mitoid INTEGER PRIMARY KEY ASC 
                sequence TEXT,
              seqlen INTEGER NOT NULL)")
   
   fa <- seqinr::read.fasta(fastaFile, as.string=TRUE, seqtype="DNA", seqonly = TRUE)
   len <- nchar(fa)
   dbSendQueryAndClear(db,
                 sprintf(
                   "INSERT INTO mtgenome (sequence, seqlen) VALUES
                    (%d, \"%s\", %d)", 
                    1, fa[[1]], len[[1]])
               )
   return(db)
}

# this loads in amplicon data
# the data are assumed to be in bed format. Since this is all one genome
# I just use columns 2 and 3 from it to give me coordinates of the PCR products
# (sans primers)
# If the alignments were to a modified reference (duplicating bits to accommodate circular alignment)
# this routine converts those back to the original units.
# thus the last amplicon may "wrap around" (eg., 16541-80 is a valid locus)
loadAmpData <- function(db, ampFile, sep="\t", kitName, append=TRUE) {
   amps <- suppressMessages( readr::read_delim(ampFile, delim=sep) )
   
   genomeLength <- dbGetQuery(db, "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]
   if (max(amps[,3]) > genomeLength) { # this can happen. We assume it's wrap-around coordinates
      amps[,3] <- amps[,3] %% genomeLength
   }
   amps.df <- tibble(kit=kitName, start=pull(amps,2), stop=pull(amps,3)) # coerce into a  tibble
   dbWriteTable(db,
                "amps", amps.df, append=append, overwrite=FALSE)
   return(db)
}






# makes the empty database
db <- makeDB(dataDir, dbName, overwrite=TRUE )
# add the reference genome (can be only 1)
loadReferenceGenome(db, paste(dataDir, refGenome, sep="/") )
# add amplicons for the precision ID panel
loadAmpData(db, paste(dataDir, precisionIDBed, sep="/"), "\t", "PrecisionID" )


           

# close the connection...
dbDisconnect(db)


