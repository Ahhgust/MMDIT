#load packages
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))

#' This function takes three vectors Position, Allele and Type from the SNPs
#' and returns a vector of variants (e.g. "73G", "73+G", 73-" )

#' @param Pos A vector of genomic positions for the SNPs (type numeric)
#' @param Allele A vector of nucleotide bases present in the SNPs (character strings). Deletions are represented by ""
#' @param Type A vector of the type of mutation represented by the SNP. Possible values include "Substitution", "Insertion"
#' and "Deletion" (charachter strings)
#'
#' @examples
#' # Snp2variant("73", "G", "Substitution")
#' # Snp2variant(c("73", "95", "146"), c("G", "", "C"), c("Substitution", "Deletion", "Insertion"))
#'
Snp2variant<-function(Pos, Allele, Type){
  tib<-dplyr::tibble(Pos=as.integer(Pos), Allele=Allele, Type=Type)
  Variant<-dplyr::case_when( # create and keep only one vector - Variant
    Type == "Substitution" ~ paste(Pos, Allele, sep=""), # if the Type is Substitution, concatenate Pos and Allele
    Type == "Insertion" ~ paste(Pos, Allele, sep = "+"), # if the Type is Insertion, concatenate Pos and Allele
    #seperated by +
    Type == "Deletion" ~ paste(Pos, "-", sep = ""), # if the Type is Deletion, concatenate Pos and -
    TRUE               ~ "?"
    )
}


#' Makes an "empop" dataframe
#'
#'
#' This function takes two vectors SampleID and Variant
#' and returns a tibble in the empop format.
#' Each row in the tibble will correspond to 1 individual (SampleID)
#' and their empop string (tab separated) will be the value in the second column
#'
#' @importFrom magrittr %>%
#' @param SampleID A vector of sample IDs (character strings)
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#'
#' @examples
#'
#' # write_mbop("NA12871", c("73G", "95-"))
#' # write_mbop(c("NA12871", "NA12871", "NA12872"), c("73G", "95-", "73G"))
#'
write_mbop<-function(SampleID, Variant){
  tib<-tibble::tibble(SampleID=SampleID, Variant=Variant)
  Empop<-tib%>%dplyr::group_by(SampleID)%>%dplyr::summarise(empopstring = paste(unique(Variant),collapse = "\t"))
}

#' reads empop

#' This function takes an empop file containing at least four (mandatory) columns
#' (refer to https://empop.online/downloads for the emp file format)
#' and returns a tibble with two columns - SampleID and Variant
#'
#' @importFrom magrittr %>%
#' @param EMPOP An EMPOP file (tab seperated)
#' @param s a numeric argument that tells the function how many rows of data to skip while reading in EMPOP file
#' @param ncol2skip number of columns to skip (starting from left) in the empop file default=3
#' @param guess_max a number; passed to read_delim. Helps with fast reading of files...
#' If no value is provided, the first line is ommited by default
#'
#' @examples
#' #   Empop2variant("EMPOP.emp")
#' #   Empop2variant("EMPOP.emp", s = 3)
Empop2variant<-function(empopFile, s = 1, ncol2skip=3, guess_max=100){
  # A major bug was found with this approach causing many SNPs to be dropped.
  # specifically, those whose column index exceeded that of the 1st individual
  #EMPOP<- readr::read_delim(empopFile, delim = "\t", skip = s, col_names = FALSE, guess_max=guess_max) %>% # read the empop files in
   # dplyr::rename(SampleID = X1, Haplogroup = X2, Frequencies = X3)
#  LongF<- tidyr::gather(EMPOP,"Key","Variant",-SampleID, -Haplogroup, -Frequencies, factor_key = FALSE)
#  if(ncol(LongF)==3) {
 #   dplyr::mutate(LongF,Variant = NA) %>% dplyr::select(-Haplogroup, -Frequencies)-> LongF
  #}
  #LongF%>%dplyr::select(SampleID, Variant) -> LongF

  # ; is not part of the empop file format and it is disallowed
  e <- readr::read_delim(empopFile, skip=s, delim=";", guess_max=guess_max, col_names=FALSE, col_types=cols())
  colnames(e)[1] <- "Variant"
  e <- tidyr::separate(e, Variant, "SampleID", extra='drop', remove=F, sep="\t")
  e <- tidyr::separate_rows(e, Variant, sep="\t") %>%
    dplyr::group_by(SampleID) %>%
    dplyr::filter(row_number() > ncol2skip) %>%
    dplyr::ungroup() %>%
    select(SampleID, Variant)
  # empty strings arise as the variant in the case that the individual == rCRS
  # 73A is what the rCRS has, so I'm going to add it in, otherwise, the downstream analyses generate a lot of NAs...
  # e.g., 73A is what the rCRS has...
  e <- mutate(e,
              Variant=ifelse(Variant=="", "73A", Variant) )

  return(e)
}

#' Converts 2 snp...
#'
#' This function takes a vector of variants (e.g. "73G", "73+G", 73-" ) and returns a tibble
#' with 3 columns : Pos, Allele and Type that contain information on the genomic position,
#' nucleotide base and type of mutation (i.e. "Substitution", "Insertion" and "Deletion") respectively for each variant
#'
#' @importFrom magrittr %>%
#'
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#'
#' @examples
#' # Variant2snp("73G")
#' # Variant2snp(c("73G", "95-", "146+C"))
#'
Variant2snp<-function(Variant){
  df<-(stringr::str_replace(Variant,"\\+", "\\."))# replace all occurences of "+" with a "."
  df1<-stringr::str_replace_all(df, "^[^0-9]","")%>% tibble::enframe()# replace any non numeric element at the start of the string with an empty string and make the result into a dataframe
  df2<-dplyr::rename(df1,Var = "value")%>% dplyr::select(-(1))# rename column to "Var" to be able to manipulate it downstream
  df3<-df2%>%
    tidyr::extract(Var, into=c("neg", "Pos", "Ins", "Allele"), "^([A-Z-])?(\\d+)(\\.\\d*)?([^.]*)") %>% # split Var into four columns
    dplyr::mutate(
           Pos=as.integer(Pos), # added by AW
           Len=ifelse(Ins==".", 1, as.integer(sub(".", "", Ins))), # create a column "Len"
           Allele = dplyr::case_when(                                     # populate the columen "Allele"...
             base::grepl("^del|-", Allele, ignore.case = T)~ "",         # with an empty string if vector contains "del" or "-"
             is.na(Len) | Len==1 ~ Allele,                         # with existing string if the Len is 1 or NA
             !is.na(Len) | Len>1 ~ str_dup(Allele, Len),           # with existing string mutliplied by the number in column Len
             TRUE ~ "?"                                            # with a "?" if none of these cases occur
           ),
           Type = dplyr::case_when(                                                         # create a column "Type" and populate it...
             base::grepl("\\.\\d|\\.", Ins)~ "Insertion",                                  # with "Insertion" if the column "Ins" contains a "." or digit
             base::grepl("^[ACGTRYSWKMBDHVN]$", Allele, ignore.case = TRUE)~"Substitution",# with "Substitution" if the column "Allele" contains a nucleotide or the IUPAC code
             Allele == "" ~ "Deletion",                                              # with "Deletion" if the column "Allele" has an empty string
             TRUE ~ "?"                                                              # with a "?" if none of these cases occur
           )
    ) %>%
    dplyr::select(-(neg), -(Ins), -(Len))                                # drop all unnecessary columns
}

#' Makes a giant empop file...
#'
#' This function takes a lookup-table file and a path to EMPOP files
#' and returns a bigger EMPOP file (i.e. all EMPOP files in the path bound together)
#' Note: This function is not for general use. It was constructed specifically
#' to construct haplotypes from HmtDB (https://www.hmtdb.uniba.it/) using custom up-stream processing.
#'
#' @importFrom magrittr %>%
#'
#'
#' @param LUT A lookup-table file (tab seperated)
#' @param Pa The path to folder where EMPOP files are stored (charachter string)
#'
#'@examples
#'
#' # Hmtdb2Empop(LUT, Pa = R.home())
#'
#'
Hmtdb2Empop<-function(LUT,Pa) {

  #function to add the amp column to the longfile2 format
  add_amp<-function(e) {
    tib<-Empop2variant(e) # run Empop2Variant
    tib$Amp<-stringr::str_extract(e, "[\\d]+") # add a column called "Amp"
    #and populate with only the digits extracted from filenames
    return(tib) #return the dataframe
  }

  l<-base::list.files(Pa, pattern = ".emp")# read in the name of empop files in current folder
  empop<-base::lapply(l, add_amp)%>% #run through the list and for each file in the list add an "Amp" column
    dplyr::bind_rows()%>% tidyr::drop_na()#  join all elements of the list into a dataframe and drop missing values


  t<-readr::read_tsv(LUT, col_names = FALSE )# read lookup table into R
  t%>%dplyr::select_if(~!(all(is.na(.)) | all(. == "")))-> t # remove empty columns
  t<-dplyr::rename(t, "Amp" = X5)#change the column name to "Amp" to be able to join with empop file later
  t$Amp<-as.character(t$Amp)# coerce number to charachter to be able to join with empop file later

  empop1<- dplyr::inner_join(t, empop, by = "Amp" ) # join longfile2 with lookuptable
  dplyr::anti_join(t, empop1, by = "X6")%>%dplyr::pull(X6) -> ids # pull out all individuals that have all amps same as rCRS

  empop2<-dplyr::mutate(empop1, Sort =(stringr::str_extract(Variant, "\\d+"))) # create a column "Sort" and populate it with just the numbers
  #from the Variant column
  empop2$Sort<-as.numeric(empop2$Sort)#convert "Sort" into type numeric
  empop3<-empop2%>%dplyr::arrange(Sort) # sort "Sort" in ascending order
  empop4<-dplyr::select(empop3, -(X1:X4), -(Amp), -(SampleID), -(Sort)) # delete unwanted columns
  #empop5<-rename(empop4, SampleID = ncol(empop4))# rename the last column to "SampleID"
  empop5<-dplyr::rename(empop4, SampleID = X6)# rename the last column to "SampleID"
  Empop<-empop5%>%dplyr::group_by(SampleID)%>%dplyr::summarise(empopstring =base::paste(unique(Variant),collapse = "\t"))#go from long to wide format
  Empop<-dplyr::select(Empop, SampleID, dplyr::everything())%>% # move "SampleID" to the beginging of dataframe
    dplyr::mutate(Haplogroup = "", Frequencies = "") %>% # add the two other needed columns to dataframe
    dplyr::select(SampleID, Haplogroup, Frequencies, dplyr::everything())# change order of the columns
  Empop <-tibble::add_row(Empop, .before = 1)%>% # add the first row for title of study and author details
    tibble::add_row(.before = 2) %>% # add the second row for geo background
    tibble::add_row(.before = 3, SampleID = "!#" ) # add third row for sequence range

  #create a tibble of indvidials with same amps as rCRS
  tibble::tibble(SampleID = unique(ids),
         Haplogroup = "",
         Frequencies = "",
         empopstring = "")-> foo

  Empop <- dplyr::bind_rows(Empop, foo)# bind the above tibble to the final EMPOP file

 utils::write.table(Empop, "Empop.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE, na="") # write tab-delimited file to file
#without column or row names
}
