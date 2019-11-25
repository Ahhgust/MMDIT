#!/usr/bin/env Rscript
# Written by August Woerner
# This creates the (empty) database
# makeDB function
#
# And it loads it with the appropriate data

# database dependencies

library(devtools)

# Note to self
# Install bit, bit64 and blob FIRST!
# the dependencies are not handled (correctly)...
#devtools::install_github("Ahhgust/RSQLite", ref="devel", upgrade=TRUE)
#devtools::install_github("Ahhgust/Haplotypical", ref="master", upgrade=TRUE)


library(Rcpp)
library(DBI)
library(RSQLite)
library(sqldf)
library(XLConnect)
library(Haplotypical)


# quality of life dependencies
suppressPackageStartupMessages( library(tibble) )
suppressPackageStartupMessages( library(readr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(magrittr) )



# for reading fasta files (reference genome)
library(seqinr)

# CONSTANTS
# data directory
dataDir <- 'mito_data'

# final database filename
dbName <- 'mmdit.sqlite3'

# reference genome
refGenome <- 'rCRS.fasta'

precisionIDBed <- "PrecisionID_mtDNA_WG_targets.bed"
default_kit <- "PrecisionID"

# helper function.
# Normally we send a query and then clear it. this wraps it up into 1 statement
dbSendQueryAndClear <- function(db, query) {
  DBI::dbClearResult( DBI::dbSendQuery(db, query) )
}

#' Creates (an empty) MMDIT database
#'
#' This *creates* an EMPTY MMDIT database
#' @param path; a directory (e.g., data)
#' @param a filename (e.g., mmdit.sqlite3)
#' @param overwrite (if there's already a file with that name, should I overwrite it?)
#' and creates a database with the requisite tables in it
#' if overwrite is FALSE, then it won't overwrite an existing database
#' in that path.
#' otherwise, it will clobber what was there (use overwrite=TRUE with caution!!)
#'
#' @return the handle to the database (throughout all other documentation, this is db)
#' @examples
#' # Make a new MMDIT database...
#'
#' # db <- makeDB(dataDir, dbName, overwrite=TRUE )
#' # add the reference genome (can be only 1)
#' # loadReferenceGenome(db, paste(dataDir, refGenome, sep="/") )
#'
makeDB <- function(path, dbFilename, overwrite=FALSE) {


  if (! dir.exists(path) ) {
    dir.create(path)
    if (! dir.exists(path)) {
      stop(
        paste("Failed to create directory: ", path , "\nThis is probably a file permission issue...")
      )
      return(NULL)
    }
  }

  dbFile <- paste(path, dbFilename, sep="/")

  if (file.exists(dbFile)) { #DB already exists
    if (overwrite==FALSE) {
      stop(
          paste("DB file " , dbFile, " already exists. Let's not overwrite things, shall we?")
          )
      return(NULL)
    }
  }

  db <- DBI::dbConnect( RSQLite::SQLite(), dbFile, loadable.extensions = TRUE)
  RSQLite::initExtension(db)


  # create the database tables:

  # information on the amplicons
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS amps")
  dbSendQueryAndClear(db,
              "CREATE TABLE amps
              (ampid INTEGER PRIMARY KEY ASC,
              kit TEXT NOT NULL,
              start INTEGER NOT NULL,
              stop INTEGER NOT NULL)" )

  # and the sample IDs (and their corresponding populations and information where they came from (provenance))
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS populations")
  dbSendQueryAndClear(db,
              "CREATE TABLE populations
              (sampleid TEXT PRIMARY KEY,
              whence TEXT,
              pop TEXT NOT NULL)" )


  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS full_mito_diffseqs")
  dbSendQueryAndClear(db,
            "CREATE TABLE full_mito_diffseqs
              (genomeid INTEGER NOT NULL,
              sampleid TEXT NOT NULL,
                 position INTEGER NOT NULL,
                 event TEXT NOT NULL,
                 basecall TEXT NOT NULL,
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                 PRIMARY KEY(genomeid)
              )" )


  # a trivial table of just the reference genome
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS mtgenome")
  dbSendQueryAndClear(db,
                      "CREATE TABLE mtgenome
              (mitoid INTEGER PRIMARY KEY ASC,
                sequence TEXT,
                seqlen INTEGER NOT NULL)")




  return(db)
}

#' Takes editdist3 distance into the unit-edit distance
#'
#' @param db the database (DBI object)
#'
#' the table that lets you specify the snp/indel costs in the editdist3 function
#' as per:
#' https://www.sqlite.org/spellfix1.html#the_editdist3_function
initEditDists <- function(db) {
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS editcost")
  dbSendQueryAndClear(db,
                    "CREATE TABLE editcost
                      (iLang INT,
                      cFrom TEXT,
                      cTo   TEXT,
                      iCost INT)")

  # makes this into a true unit editdist (all costs are 1)
  tib <- tibble::tribble(
    ~iLang, ~cFrom, ~cTo, ~iCost,
       0,     '',    '?',  1,
       0,     '?',    '',  1,
       0,     '?',  '?',   1)


  DBI::dbWriteTable(db,
    "editcost",
    tib, append=TRUE)
}
#' Use this to add the rCRS sequence to MMDIT
#' It adds the record to the table mtgenome
#' this table has *1* row in it. The reference genome.
#' If you call this function multiple times, the last time
#' will be the record that "sticks"
#'
#' Technically any genome would work here, BUT
#' the difference encodings MUST be with respect to this reference genome!!
#' i.e., in all likelihood this should *ONLY* be the rCRS sequence.
#'
#' @param db (database handle)
#' @param fastafile (must be both the name of the file and the path to open it. This is passed directly to seqinr::read.fasta
#'
#' @examples
#'
#' # loadReferenceGenome(db, paste(dataDir, refGenome, sep="/") )
#'
loadReferenceGenome <- function(db, fastaFile) {

   fa <- seqinr::read.fasta(fastaFile, as.string=TRUE, seqtype="DNA", seqonly = TRUE)
   len <- nchar(fa)
   dbSendQueryAndClear(db,
                 sprintf(
                   "INSERT INTO mtgenome (mitoid, sequence, seqlen) VALUES
                    (%d, \"%s\", %d)",
                    1, fa[[1]], len[[1]])
               )
   return(db)
}

#' Creates temp table of haplotypes
#'
#' this creates a temporary table called haplotypes
#' The haplotypes table converts the difference encodings into strings (after applying the mask)
#'
#' #' @param db (database handle)
makeHaplotypeTable <- function(db) {
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS haplotypes")

  dbSendQueryAndClear(db,
                      "CREATE TEMP TABLE haplotypes
              (hapid INTEGER NOT NULL,
               ampid INTEGER NOT NULL,
               sampleid TEXT NOT NULL,
               haplotype TEXT NOT NULL,
               haplotype_count INTEGER DEFAULT 1,
                 FOREIGN KEY (ampid) REFERENCES amps(ampid),
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid),
                 PRIMARY KEY(hapid)
               )" )
  return(db)
}

#' this loads in amplicon data (i.e., the coordinates of the amps)
#' into the MITO db.
#' If the alignments were to a modified reference (duplicating bits to accommodate circular alignment)
#' this routine converts those back to the original units.
#' thus the last amplicon may "wrap around" (eg., 16541-80 is a valid locus)
#'
#' @param db (the database connection)
#' @param ampfile (a BED file of amplicons; columns 2 and 3 are used only)
#' @param kitname (the NAME of the kit that this corresponds to (e.g., PrecisionID) )
#' @param sep (\"\\t\" for a tab-separated file. Anything else is a little wacky... no longer a bed file)
#' @param append (whether or not the amplicons are APPENDED to the list of current amplicons; set to FALSE to overwrite)
#'
#' @return the database handle (not necessary to keep this)
#'
#' @examples
#' # One way to add another kit:
#'
#' # loadAmpData(db, "data/Mybedfile.bed", kit="SomeOtherKit", append=TRUE)
#'
#'
loadAmpData <- function(db, ampFile,  kitName, sep="\t", append=TRUE) {
   amps <- suppressMessages( readr::read_delim(ampFile, delim=sep) )

   overwrite<-FALSE
   if (append == FALSE) {
     overwrite<-TRUE
   }

   if (! overwrite) {
     redundantCount <-
       DBI::dbGetQuery(db, sprintf(
       "SELECT count(kit) FROM amps WHERE kit == '%s'", kitName))

    if (redundantCount[[1]] > 0) {
        stop(
          sprintf("Kit %s is already in the database...\nYou can't add it twice!", kitName)
        )
      return(db)
    }
   }
   genomeLength <- DBI::dbGetQuery(db, "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]
   if (max(amps[,3]) > genomeLength) { # this can happen. We assume it's wrap-around coordinates
      amps[,3] <- amps[,3] %% genomeLength
   }
   amps.df <- tibble::tibble(kit=kitName, start=dplyr::pull(amps,2), stop=dplyr::pull(amps,3)) # coerce into a  tibble
   DBI::dbWriteTable(db,
                "amps", amps.df, append=append, overwrite=overwrite)
   return(db)
}


#' Amplicon sequences of reference
#'
#' Generates the amplicon sequences
#' based on the reference sequence.
#' For use with MMDIT only!
#'
#' @param db the database connection
#' @param kit the kit name (must be a valid KIT in the amps table)
#' @return returns data frame with 4 columns: the amp ID, the start coordinate (0-based), the stop coordinate (1-based), and the Seq(uence) of the rcrs for the amp
#'
#' @example
#' # amps <- getRefAmpSeqs(db, kit="PrecisionID")
#'
getRefAmpSeqs <- function(db, kit=default_kit) {
  loc <- DBI::dbGetQuery(db,
                    sprintf(
                    "SELECT ampid, start, stop
                    FROM
                     amps
                    WHERE
                    kit= '%s'", kit) )

  if (nrow(loc) == 0) {
     stop(
       paste(
       "No data for kit ", kit , " was found...")
     )
  }
  # get the mitochondrial genome (concatenated with itself)
  # the self-concatentation allows for the easy lookup of amplicons that span the
  # first and last bases of the linearized mitochondrial genome.
  mtgenomeSeq <- stringr::str_dup(
    DBI::dbGetQuery(db,
                            "SELECT sequence FROM mtgenome LIMIT 1")[[1]],
    2)

  mtgenomeLen <- DBI::dbGetQuery(db,
                            "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]

  # recall substr using 1-based indexing
  # and the coordinates are bed format (0-based start, 1-based stop)
  loc$Seq <-
      stringr::str_sub(mtgenomeSeq,
                       loc$start + 1,
                       ifelse( loc$start < loc$stop,
                         loc$stop,
                         loc$stop + mtgenomeLen))

  loc
}

#' Trivial testing to see if the editdist3 function works...
#' Takes no arguments; returns a data frame (if editdist3 works!)...
#' If this *doesn't* work you need to make sure you're
#' using the development version of rsqlite as downloaded from my github at:
#'
#' devtools::install_github("Ahhgust/RSQLite", ref="devel", upgrade=TRUE)
#'
testEditdist <- function() {
  db <- RSQLite::datasetsDb()
  RSQLite::initExtension(db)
  # 100 == insertion, 150 is used for substitution
  DBI::dbGetQuery(db,"SELECT *, editdist3(row_names, 'Merc 230') AS editdist
             FROM mtcars WHERE editdist3(row_names, 'Merc 230') < 251") -> tib
  DBI::dbDisconnect(db)
  return(tib)
}

getSeqdiffCodes <- function() {
  c(X=1, I=2,D=3)
}

initDB <- function() {
# makes the empty database
  db <- makeDB(dataDir, dbName, overwrite=TRUE )
# add the reference genome (can be only 1)
  loadReferenceGenome(db, paste(dataDir, refGenome, sep="/") )
# add amplicons for the precision ID panel
  loadAmpData(db, paste(dataDir, precisionIDBed, sep="/"), kitName=default_kit, sep="\t" )

# make the edit distance a true editdist
  initEditDists(db)

  refAmps <- getRefAmpSeqs(db)
  return(db)
}
# close the connection...
#dbDisconnect(db)


