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
library(stringr)
library(Haplotypical)
library(stringdist)

# for reading fasta files (reference genome)
library(seqinr)

# quality of life dependencies
suppressPackageStartupMessages( library(tibble) )
suppressPackageStartupMessages( library(readr) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(magrittr) )





# CONSTANTS
# data directory
dataDir <- 'mito_data'

# final database filename
dbName <- 'mmdit.sqlite3'

# reference genome
refGenome <- 'rCRS.fasta'

precisionIDBed <- "PrecisionID_mtDNA_WG_targets.bed"
default_kit <- "PrecisionID"

hmtFilename <- "Empop.tsv"

# helper function.
# Normally we send a query and then clear it. this wraps it up into 1 statement
dbSendQueryAndClear <- function(db, query) {
  DBI::dbClearResult( DBI::dbSendQuery(db, query) )
}

#' Creates (an empty) MMDIT database
#'
#' This *creates* an EMPTY MMDIT database
#' @param path a directory (e.g., data)
#' @param dbFilename (e.g., mmdit.sqlite3)
#' @param overwrite (if there's already a file with that name, should I overwrite it?)
#' and creates a database with the requisite tables in it
#' if overwrite is FALSE, then it won't overwrite an existing database
#' in that path.
#' otherwise, it will clobber what was there (use overwrite=TRUE with caution!!)
#'
#' @return the handle to the database (throughout all other documentation, this is db)
#' @export
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

# primary key is the rowid (default in SQlite)
  dbSendQueryAndClear(db, "DROP TABLE IF EXISTS full_mito_diffseqs")
  dbSendQueryAndClear(db,
            "CREATE TABLE full_mito_diffseqs
              (
              sampleid TEXT NOT NULL,
                 position INTEGER NOT NULL,
                 event TEXT NOT NULL,
                 basecall TEXT NOT NULL,
                 FOREIGN KEY (sampleid) REFERENCES populations(sampleid)
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


#' Loads the database...
#'
#' This *loads* the MMDIT database
#' The haplotypes table will be empty (but fillable)
#' @param path a directory (e.g., data)
#' @param dbFilename (e.g., mmdit.sqlite3)
#'
#' @export
loadMMDIT <- function(path=dataDir, dbFilename=dbName) {
  dbFile <- paste(path, dbFilename, sep="/")

  if (file.exists(dbFile)) {
    db <- DBI::dbConnect( RSQLite::SQLite(), dbFile, loadable.extensions = TRUE)
    RSQLite::initExtension(db)
    return(db)
  }
  return(NULL)
}

#' Takes editdist3 distance into the unit-edit distance
#'
#' @param db the database (DBI object)
#' @param stops integer vector of iLangs (stop coords)
#' @param tbName name of cost table (used will spellfix1)
#' the table that lets you specify the snp/indel costs in the editdist3 function
#' as per:
#' https://www.sqlite.org/spellfix1.html#the_editdist3_function
initEditDists <- function(db, stops=c(1), tbName="EDITCOST") {
  dbSendQueryAndClear(db,
                      paste0("DROP TABLE IF EXISTS ",
                             tbName) )
  dbSendQueryAndClear(db,
                        paste0(
                        "CREATE TABLE " ,
                          tbName,
                        " (iLang INT,
                          cFrom TEXT,
                          cTo   TEXT,
                          iCost INT)"
                        )
                      )

  # makes this into a true unit editdist (all costs are 1)
  # the same edit dists are made for each iLang
    tib <- data.frame(
      iLang=rep(unique(stops), each=3),
      cFrom=c('', '?','?'),
      cTo=c('?','', '?'),
      iCost=c(1))

  DBI::dbWriteTable(db,
    tbName,
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
#' @param fastaFile (must be both the name of the file and the path to open it. This is passed directly to seqinr::read.fasta
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
#' @importFrom magrittr %>%
#' @param db (database handle)
#' @param sampleid (individual sample IDs from MMDIT)
#' @param stopCoord (stop coordinate of amplicon)
#' @param seq (the DNA sequence of this individual for this amplicon)
#' @param tableName (the name of the table to be made; default="haplotypes")
#' @param vtableName (the name of the virtual spellfix table to be made; default="sequences")
#'
makeHaplotypeTable <- function(db, sampleid, stopCoord, seq, tableName="haplotypes", vtableName="sequences") {

# initialize table of edit distances
  # changes the default to the unit edit distance
  # and sets up iLangs
  # which allow multiple amplicons to be put into the same table
  # the default edit distance table name is EDITCOST
  # initEditDists(db, stopCoord)

  dbSendQueryAndClear(db,
                      paste0("DROP TABLE IF EXISTS " , tableName))

  dbSendQueryAndClear(db,
                      paste0("DROP TABLE IF EXISTS " , vtableName))


#  dbSendQueryAndClear(db,
 #       paste0("CREATE VIRTUAL TABLE ", vtableName ,
  #            " USING spellfix1(edit_cost_table=EDITCOST)" )
   #            )

  tib <- tibble::tibble(sampleid=sampleid,
                        stop=as.integer(stopCoord),
                        seq=seq) # seqid TBD

  # make a tibble of each unique haplotype for each unique amplicon
  sngl <- dplyr::group_by(tib, stop, seq) %>%
            dplyr::filter(dplyr::row_number()==1) %>%
              dplyr::ungroup()

  # and make an index (an int) that identifies each hap
  sngl$seqid <- 1:nrow(sngl)

  # now tib has a seq id
  dplyr::left_join(tib,
                   dplyr::select(sngl, stop, seq, seqid),
                   by=c("stop", "seq")
                   ) -> tib

  DBI::dbWriteTable(db, tableName,
                    dplyr::select(tib, sampleid, stop, seqid),
                    append=FALSE, overwrite=TRUE,temporary=TRUE,
                    field.types=c("sampleid"="character", "stop"="integer", "seqid"="integer"))

#TODO: get ID from query. can't assign ID
#  DBI::dbExecute(db,
 #                    paste0('INSERT INTO ', vtableName , '(word, langid) VALUES (?, ?)'),
  #                   params=list(sngl$seq, sngl$stop) )



  DBI::dbClearResult(
    DBI::dbSendQuery(db,
                  paste0(
                    "CREATE INDEX stopindex ON " ,
                    tableName ,
                    " (stop, seqid)" )
                  )
  )

  return(db)
}

#' Creates temp table of *sample* haplotypes
#'
#' this creates a temporary table called refamps
#' The haplotypes table converts the difference encodings into strings (after applying the mask)
#'
#' @param db (database handle)
#' @param stop (stop coordinate of amplicon)
#' @param sequence (the DNA sequence of this individual for this amplicon)
#' @param tableName (the name of the table to be made)
makeSampleTable <- function(db, stopCoord, seq, tableName="sampleamps") {

  if (DBI::dbExistsTable(db, tableName)) {
    DBI::dbRemoveTable(db, tableName, temporary=TRUE)
  }

  tib <- tibble::tibble(stop=as.integer(stopCoord), seq=seq)
  DBI::dbWriteTable(db, tableName, tib,  append=FALSE, overwrite=TRUE,temporary=TRUE,
                    field.types=c("stop"="integer", "seq"="character"))

  return(db)
}

#' this loads in amplicon data (i.e., the coordinates of the amps)
#' into the MITO db.
#' If the alignments were to a modified reference (duplicating bits to accommodate circular alignment)
#' this routine converts those back to the original units.
#' thus the last amplicon may "wrap around" (eg., 16541-80 is a valid locus)
#'
#' @param db (the database connection)
#' @param ampFile (a BED file of amplicons; columns 2 and 3 are used only)
#' @param kitName (the NAME of the kit that this corresponds to (e.g., PrecisionID) )
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
   amps.df <- tibble::tibble(kit=kitName,
                             start= as.integer(dplyr::pull(amps,2)),
                             stop = as.integer(dplyr::pull(amps,3))) # coerce into a  tibble
   DBI::dbWriteTable(db,
                "amps", amps.df, append=append, overwrite=overwrite,
                field.types=c("kit"="character", "start"="integer","stop"="integer")
                )
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
#' @examples
#' # amps <- getRefAmpSeqs(db, kit="PrecisionID")
#' @export
getRefAmpSeqs <- function(db, kit=default_kit) {
  #TODO: Take in start/stop coordinates as well as a kit.
  loc <- DBI::dbGetQuery(db,
                    sprintf(
                    "SELECT start, stop
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





#' get all seqdiffs for some population
#'
#' @param db the database handle
#' @param pop the population requested (defaults to all). vectorized populations okay
#' @param ignoreIndels strips out indel events (default TRUE)
#' @param kit nameof amplicon sequencing kit...
#' @export
getSeqdiffs <- function(db, pop="%", ignoreIndels=FALSE, kit=default_kit) {

# optionally we want to filter out indels...
  indelString <- ""
  if (ignoreIndels) {
    indelString <- " AND event = 'X'"
  }
  genomeLength <- DBI::dbGetQuery(db, "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]

  diffs <- DBI::dbGetQuery(db,
                            paste0(
                              "SELECT populations.sampleid, position, event, basecall, amps.start, amps.stop ",
                              "FROM   full_mito_diffseqs, populations, amps ",
                              "WHERE  populations.sampleid  = full_mito_diffseqs.sampleid ",
                              "AND amps.kit == '" , kit , "' ",
                              "AND (position > amps.start AND position <= amps.stop OR ",
                                   "position + " , genomeLength , " >amps.start AND position +" , genomeLength, " <= amps.stop)" ,
                              "AND (",
                                      "pop LIKE '",
                                      stringr::str_c(pop, collapse="' OR pop LIKE '"),
                                      "')" ,
                              indelString
                             )
  )

  return(diffs)
}

#' Length of mito-genome (reference)
#'
#' This is the non-circularized (concatenated)
#' mtgenome length (16569 unless some real funny-business is happening)
#'
#' @param db the database handle
#'
#' @export
getMtgenomeLength <- function(db) {
  DBI::dbGetQuery(db,
                  "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]
}

#' Mito-genome (reference) sequence
#'
#' Returns the self-concatenated mitochondrial reference genome
#' Self-concatenation allows us to handle the circular alignment problem
#' (eg, reads, amplicons that span the canonical "start" of the mt-genome)
#' and coordinates from a circular alignment
#'
#' @param db the database handle
#' @export
getMtgenomeSequence <- function(db) {
  stringr::str_dup(
    DBI::dbGetQuery(db,
                    "SELECT sequence FROM mtgenome LIMIT 1")[[1]],
    2)
}

#' generates amplicon sequence data from difference encodings

#' @useDynLib MMDIT
#' @importFrom magrittr %>%
#' @param db the database handle
#' @param pop the population requested (defaults to all). vectorized populations okay
#' @param ignoreIndels strips out indel events (default TRUE)
#' @param kit the amplicon kit
#' @param blk blacklist of sites to filter out!
#' @export
getAmps <- function(db, pop='%', ignoreIndels=FALSE, kit=default_kit, blk=c()) {

  mtgenomeLen <- DBI::dbGetQuery(db,
                                 "SELECT seqlen FROM mtgenome LIMIT 1")[[1]]


  diffs <- dplyr::filter(
      getSeqdiffs(db, pop, ignoreIndels),
      ! position %in% blk,
      ! position %in% (blk+mtgenomeLen) )

  refamps <-   getRefAmpSeqs(db, kit)
  lookupcodes <- getSeqdiffCodes()

  # for each individual / amplicon (with a difference to the rCRS)
  #
  dplyr::group_by(diffs,
      sampleid,
      start,
      stop
    ) %>%
    dplyr::inner_join( # grab the reference sequence associated with that amplicon
      refamps,
      by=c("start", "stop")
    ) %>%
    dplyr::mutate( # deal with circular coordinates, compute the coordinate of the SNP
      # within the amplicon
      position=as.integer(position),
      relpos=as.integer(ifelse(
        position > start,
        position - start,
        position + mtgenomeLen - start))
    ) %>%
    dplyr::arrange(position) %>%
    summarize( # and interpolate; make the string based on the string differences
      sequence=seqdiffs2seq(
        Seq[[1]],
        relpos,
        lookupcodes[ event ],
        basecall
        )
      ) -> foo
# at this point, foo only has the sequence differences.
# any amp that is == to the rCRS is not present.
# we add a no-op sequence difference for individuals that == the rCRS
# (to ensure the individual isn't lost with the join)


  base::expand.grid( # first, make a data frame with all pairs of individuals and stop coordinates
    sampleid=unique(diffs$sampleid),
    stop=refamps$stop, stringsAsFactors = FALSE) %>%
    dplyr::left_join(refamps, # and add in the start and reference Seq
              by=c("stop")) %>%
    dplyr::left_join(foo, # and left-join with the sequence differences
              by=c("sampleid", "start", "stop")) %>%
    # from here, if there was a difference encoding for the amp, that amps sequence is 'sequence'
    # otherwise it's NA (no differences to reference)
    dplyr::mutate(sequence=ifelse(is.na(sequence), Seq, sequence)) %>% # change NAs to reference
    dplyr::select(sampleid, start, stop, sequence) -> amps # and reorder the columns


  return(amps)
}

#' finds nearest neighbors
#'
#' @param alldist the output from getAllDistances
#' @param knn the number of nearest neighbors to find (ties broken arbitrarily)
#' @param maxdist if >-1, returns all neighbors within distance <= maxdist (ignores knn argument)
#' @param useIdentityFunction counts the number of amplicons that differ, not the sum of the distances
getNN <- function(alldist, knn=20, maxdist=-1, useIdentityFunction=FALSE) {
   dplyr::group_by(alldist, sampleid) %>%
    dplyr::summarize(DistSum=sum(dist), SumMismatch=sum(dist>0)) %>%
    dplyr::arrange(DistSum) %>%
    dplyr::ungroup() -> bydist
   if (maxdist < 0) {
    return(  head(bydist, n=knn) )
   }
   if (useIdentityFunction) {
    return(  dplyr::filter(bydist, SumMismatch <= maxdist) )
   }
  return(  dplyr::filter(bydist, DistSum <= maxdist) )
}

#' solves the All 1NN
#'
#' (much) more efficient solutions exist (cover trees)
#' O(n log n) for n individuals and a constant genome length
#' the solution provided is O(n^2) by the same count.
#'
#' This looks at each individual in amps
#' and finds the 2NN of that individual and the associated distances
#' and the selects the 2nd individual. in theory this could give you
#' the same individual as you queried (ties are broken arbitrarily).
#' If you care about that... write your own :)
#' I just want the distances!
#'
#' @param amps the output from getAmps
#' @param ignoreIndels (whether or not indels were used to constitute the strings...)
naiveGetAllNN <- function(amps, ignoreIndels=FALSE) {
  dplyr::group_by(amps, sampleid) %>%
    dplyr::summarize(NNDist=
                list(
                  getNN(
                    getAllDistances(amps, sequence,stop, ignoreIndels=ignoreIndels),
                    knn=2)[2,]

                )
    )

}

#' solves the All 1NN
#'
#' (much) more efficient solutions exist (cover trees)
#' O(n log n) for n individuals and a constant genome length
#' the solution provided is O(n^2) by the same count.
#' This version of the code finds the NN for each amplicon (not for the whole mito)
#'
#' This looks at each individual in amps
#' and finds the 2NN of that individual and the associated distances
#' and the selects the 2nd individual. in theory this could give you
#' the same individual as you queried (ties are broken arbitrarily).
#' If you care about that... write your own :)
#' I just want the distances!
#'
#' @param amps the output from getAmps
#' @param ignoreIndels (whether or not indels were used to constitute the strings...)
naiveGetAllNNBySampAndStop <- function(amps, ignoreIndels=FALSE) {
  amps %>%
  dplyr::group_by(sampleid, stop) %>%
    dplyr::summarize(NNDist=
                       list(
                           getAllDistances(dplyr::filter(amps, stop==stop), sequence,stop, ignoreIndels=ignoreIndels)[2,]
                        )
    ) %>% tidyr::unnest()


}


writeOneNNTib <- function(amps) {
    amps %>% naiveGetAllNN(TRUE) -> all1NNTib
  readr::write_tsv( tidyr::unnest(all1NNTib), "All1NN.tsv")
  amps %>% naiveGetAllNNBySampAndStop(TRUE) -> all1NNByAmp
  readr::write_tsv(all1NNByAmp, "All1NN_byAmp.tsv")


}

#' generates amplicon sequence data from difference encodings

#' @param amps the output from getAmps
#' @param seqs query sequences
#' @param stops the stop coordinate of said sequences
#' @param ignoreIndels boolean; must match amps
getAllDistances <- function(amps, seqs, stops, ignoreIndels=FALSE) {

  tmp <- tibble::tibble(QuerySeq=seqs, stop=stops)
  dplyr::left_join(tmp,
                   amps,
                   by=c("stop")) ->cmbnd

  type <- "lv" # levenstein
  if (ignoreIndels) { # choose the appropriate distance function...
    type <- "hamming"
  }
  cmbnd$dist <- stringdist::stringdist(cmbnd$QuerySeq, cmbnd$sequence, method=type)
  return(dplyr::arrange(cmbnd, stop, dist))
}

#' gets population labels
#'
#' This takes in the database handle
#' and returns a the unique populations
#'
#' @export
getPops <- function(db) {
  DBI::dbGetQuery(db, "SELECT DISTINCT(pop) FROM populations")
}

#' reads empop; adds haplotypes to database

#' This function takes an empop file containing at least four (mandatory) columns
#' (refer to https://empop.online/downloads for the emp file format)
#' and returns a tibble with two columns - SampleID and Variant
#'
#' @importFrom magrittr %>%
#' @param db database handle to MMDIT
#' @param empopFile a file of empop sequences...
#'
addKnownHaplotypes <- function(db, empopFile) {
   hmt <- Empop2variant( empopFile, s=3 )
   hmt <- tidyr::separate(hmt,
     SampleID, c("Population"), remove=F, extra='drop')

   dplyr::bind_cols(
     hmt,
     Variant2snp(hmt$Variant)) -> hmt
   dplyr::mutate(hmt,
          Type=
            ifelse(Type=="Substitution", "X",
                       substr(Type, 1, 1)) # insertion -> I, deletion to D
   ) -> hmt

  samps <-  DBI::dbGetQuery(db,
                               "SELECT DISTINCT sampleid FROM populations")
  if (samps[1,] %in% hmt$SampleID) {
    stop("Duplicate sample IDs found!")
    return(0)
  }

  dplyr::group_by(hmt, SampleID, Population) %>%
    dplyr::summarize(whence=empopFile) %>%
    dplyr::ungroup() %>%
    dplyr::select(sampleid=SampleID, whence, pop=Population) -> samps2add

  # add the individuals to the DB
  DBI::dbWriteTable(db, "populations", samps2add, append=TRUE)
  # and add the difference encodings

  dplyr::select(hmt, sampleid=SampleID,
              position=Pos,
              event=Type,
              basecall=Allele) -> toAdd

  # add the individuals to the DB
  DBI::dbWriteTable(db, "full_mito_diffseqs", toAdd, append=TRUE)


  return(1)
}



