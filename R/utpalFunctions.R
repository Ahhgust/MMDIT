#load packages
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(tibble))


#' @title
#' A helper function to runDeploid()that creates panel and plaf dataframes
#'
#' @description
#' This function takes a long dataframe (containing differences in coding compared to the rCRS,
#' in a select set of individuals) and returns panel and plaf dataframes
#'
#' @param long dataframe with information on sample ID, position and allele for a set of individuals
#'
#' @importFrom magrittr %>%
#' @export
#' @examples
#' createPlafPan(long)
createPlafPan<-function(EmpopLong) {
  Wide<-tidyr::spread(EmpopLong, Sample.name, Nuc)#go from long format to wide again
  Wide$POS<-as.numeric(Wide$POS)#coerce position to number to be able to sort
  Wide<-dplyr::arrange(Wide, POS)#sort by position

  # final panel formating
  Wide[,2:ncol(Wide)]<-base::ifelse(is.na(Wide[ ,-1])==TRUE,0,1)#replace NA and AGTC with 0 and  1 respectively
  panelFile<-Wide %>%dplyr::mutate(CHROM = "CHRM")#rename object, add mandatory CHRM column
  panelFile<-panelFile[,c(ncol(panelFile),1:ncol(panelFile)-1)]#and reorder


  EmpopLong %>%dplyr::group_by(POS)%>%
    dplyr::summarise(PLAC=length(Nuc), #to get po level allele count or numerator
                     u=length(base::unique(Nuc))) %>% # pick out the number of unique alleles
    dplyr::filter(u <=2) ->Numerator # filter out multi-allelic positions
  Numerator$POS<-as.numeric(Numerator$POS)# convert to numeric to allow sorting
  plafFile<-dplyr::arrange(Numerator, POS) # sort by position
  #add plaf column with correction
  plafFile<-dplyr::mutate(plafFile,PLAF = base::ifelse(PLAC==0|1,
                                                       1-(PLAC/(base::length(base::unique(EmpopLong$Sample.name))+2)), # correction formula
                                                       PLAC/base::length(base::unique(EmpopLong$Sample.name))))
  #final PLAF formating
  plafFile[-c(2,3)] ->plafFile # remove unnecesary columns
  plafFile<-plafFile %>%dplyr::mutate(CHROM = "CHRM")#rename object, add mandatory CHRM column
  plafFile<-plafFile[,c(base::ncol(plafFile),1:base::ncol(plafFile)-1)] #and reorder
  return(base::list(plafFile, panelFile))# return both panel and plaf dataframes as a list

}

#' @title
#' A helper function to runDeploid()that harmonizes the panel
#'
#' @description
#' This function takes a panel (createPlafPan object),
#' and a vector q of sites from am AltRef (modified createAltRef object)
#' and matches or "harmonizes" the number of sites between both dataframes
#'
#' @importFrom magrittr %>%
#' @export
#' @param panel a dataframe representing a reference panel for a set of individuals (createPlafPan object)
#' @param q vector of sites from an AltRef dataframe with reference and alternate allele counts (modified createAltRef object)
#'
#' @examples
#' harmonize(plaf,q)

#Function for harmonizing the panel so that it has the same number of sites as
#other dEploid input
harmonizePa<-function(panel,q) {

  pan<-dplyr::full_join(panel,q)%>%dplyr::arrange(POS)# join panel and POS from alt ref keeping all rows
  pan<-dplyr::mutate(pan, CHROM=ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  pan[base::is.na(pan)]<-0# convert all remaing NAs in the table to 0s
  CVPa<-pan # save final panel dataframe to enviroment

  return(CVPa) # return results of function
}

#' @title
#' A helper function to runDeploid()that harmonizes the plaf
#'
#' @description
#' This function takes a plaf (createPlafPan object),
#' and a vector q of sites from am AltRef (modified createAltRef object)
#' and matches or "harmonizes" the number of sites between both dataframes
#'
#' @importFrom magrittr %>%
#'
#' @param plaf a dataframe of population level allele frequencies (createPlanPlaf object)
#' @param q vector of sites from an AltRef dataframe with reference and alternate allele counts (modified createAltRef object)
#' @export
#' @examples
#' harmonize(plaf,q)
#Function for harmonizing the plaf so that it has the same number of sites as
#other dEploid input
harmonizePl<-function(Plaf,q) {
  plaf<-dplyr::full_join(Plaf,q)%>% dplyr::arrange(POS)# join plaf and POS from alt ref keeping all rows

  plaf<-dplyr::mutate(plaf, CHROM=base::ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  plaf[base::is.na(plaf)]<-1e-6# convert all remaing NAs in the table to 1s
  CVPl<-plaf # save final plaf table to enviroment

  return(CVPl)# return results of function
}


#' @title
#' A helper function to runDeploid()that harmonizes the AltRef
#'
#' @description
#' This function takes an AltRef (modified createAltRef object),
#' and a vector p of sites from a plaf/pan (createPlafPan object)
#' and matches or "harmonizes" the number of sites between both dataframes
#'
#' @importFrom magrittr %>%
#'
#' @param AltRef a dataframe of reference and alternate allele counts (modified createAltRef object)
#' @param p vector of sites from a plaf/pan dataframe (createPlafPan object)
#'
#' @examples
#' harmonize(AltRef,p)
#' @export
#Function for harmonizing the AltRef so that it has the same number of sites as
#other dEploid input
harmonizeAR<-function(AltRef, p){

  readCount <- AltRef$AltCount[1] + AltRef$RefCount[1]# create object representing total count

  AR<-dplyr::full_join(AltRef, p, by="POS")%>% dplyr::arrange(POS)# join alt ref tabel and POS from panel keeping all rows
  AR<-dplyr::mutate(AR, CHROM= base::ifelse(is.na(CHROM), "CHRM", "CHRM"))# convert all NAs in the CHROM column to CHRM
  AR$AltCount%>%tidyr::replace_na(0)->AR$AltCount # convert all NAs in ALT column to 0s
  AR$RefCount%>%tidyr::replace_na(readCount)->AR$RefCount # convert all NAs in REF column to 100s

  return(AR) # return results of function
}

#' @title
#' Get variants from SNP data
#'
#' @description
#' This function takes three vectors Position, Allele and Type from the SNPs
#' and returns a vector of variants (e.g. "73G", "73+G", 73-" )

#' @param Pos A vector of genomic positions for the SNPs (type numeric)
#' @param Allele A vector of nucleotide bases present in the SNPs (character strings). Deletions are represented by ""
#' @param Type A vector of the type of mutation represented by the SNP. Possible values include "Substitution", "Insertion"
#' and "Deletion" (charachter strings)
#'
#' @examples
#' Snp2variant("73", "G", "Substitution")
#' Snp2variant(c("73", "95", "146"), c("G", "", "C"), c("Substitution", "Deletion", "Insertion"))
#' @export
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


#' @title
#' Write inforatuon on sample ID and variants in the emp file format
#'
#' @description
#' This function takes two vectors SampleID and Variant
#' and returns a tibble in the emp file format.
#' Each row in the tibble will correspond to 1 individual (SampleID)
#' and their empop string (tab separated) will be the value in the second column
#'
#' @importFrom magrittr %>%
#'
#' @param SampleID A vector of sample IDs (character strings)
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#'
#' @examples
#' write_mbop("NA12871", c("73G", "95-"))
#' write_mbop(c("NA12871", "NA12871", "NA12872"), c("73G", "95-", "73G"))
#' @export
write_mbop<-function(SampleID, Variant){
  tib<-tibble::tibble(SampleID=SampleID, Variant=Variant)
  Empop<-tib%>%dplyr::group_by(SampleID)%>%dplyr::summarise(empopstring = paste(unique(Variant),collapse = "\t"))
}

#' @title
#' Get variants from an EMPOP file
#'
#' @description
#' This function takes an empop file containing at least four (mandatory) columns
#' (refer to https://empop.online/downloads for the emp file format)
#' and returns a tibble with two columns - SampleID and Variant
#' All individuals in the empop file must be single-source, with *no* heteroplasmies
#' We recommend taking the major allele in the case of heteroplasmy.
#' @export
#' @importFrom magrittr %>%
#' @param empopFile An EMPOP file (tab seperated)
#' @param s a numeric argument that tells the function how many rows of data to skip while reading in EMPOP file
#' @param ncol2skip number of columns to skip (starting from left) in the empop file default=3
#' @param guess_max a number; passed to read_delim. Helps with fast reading of files...
#' If no value is provided, the first line is ommited by default
#'
#' @examples
#' #   Empop2variant("EMPOP.emp")
#' #   Empop2variant("EMPOP.emp", s = 3)
Empop2variant<-function(empopFile, s = 1, ncol2skip=3, guess_max=100){

  e <- readr::read_delim(empopFile, skip=s, delim="\\0", guess_max=guess_max, col_names=FALSE, col_types=readr::cols())
  colnames(e)[1] <- "Variant"
  e <- tidyr::separate(e, Variant, "SampleID", extra='drop', remove=F, sep="\t") # get the sampleiD (1st column, tab delim)
  e <- tidyr::separate_rows(e, Variant, sep="\t") %>%  # and separate out the rest of the tabs (long-format)
    dplyr::group_by(SampleID) %>%
    dplyr::filter(dplyr::row_number() > ncol2skip) %>% # removing the first 3 columns in the original encoding
    dplyr::ungroup() %>%
    dplyr::select(SampleID, Variant)
  # empty strings arise as the variant in the case that the individual == rCRS
  # 73A is what the rCRS has, so I'm going to add it in, otherwise, the downstream analyses generate a lot of NAs...
  # e.g., 73A is what the rCRS has...
  e <- dplyr::mutate(e,
              Variant=ifelse(Variant=="", "73A", Variant) )

  return(e)
}

#' @title
#' Get SNPs from variant data
#'
#' @description
#' This function takes a vector of variants (e.g. "73G", "73+G", 73-" ) and returns a tibble
#' with three columns : Pos, Allele and Type that contain information on the genomic position,
#' nucleotide base and type of mutation (i.e. "Substitution", "Insertion" and "Deletion")resecptively for each variant
#'
#' @importFrom magrittr %>%
#'
#' @param Variant A vector of variants (character strings). See Snp2variant for description
#'
#' @examples
#' # Variant2snp("73G")
#' # Variant2snp(c("73G", "95-", "146+C"))
#' @export
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

#' @title
#' Join mutliple EMPOP files into one large EMPOP file using a look-up table
#'
#' @description
#' This function takes a lookup-table file and a path to EMPOP files
#' and returns a bigger EMPOP file (i.e. all EMPOP files in the path bound together)
#'
#' @importFrom magrittr %>%
#' @export
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

#' @title
#' Wrapper function for dEploid
#'
#' @description
#' This function takes information on alternate and reference allele counts,
#' number of individuals and the recombination rate and runs
#' the dEploid function on them. It retuns the dEploid output
#' with includes estimated haplotypes, mixture proprotions and other phasing statistics.
#' Additionaly this function also returns a list of nucleotide positions that violate the infinte sites model.
#'
#' @importFrom magrittr %>%
#'
#' @param AR Table of alternate and reference allele counts (Empop2AltRef object)
#' @param long Table with information on sample ID, position, allele and type of mutation
#' @param NumMCMC A numeric argument that tells the function the number of MCMC samples (default value 800)
#' @param exportPostProb Save the posterior probabilities of the final iteration of all individuals
#' @param recomb numeric argument that gives the function the constant recombination probability (default value of 0.0)
#' @param k A numeric argument that tells the function the number of individuals in the mixture (default value 5, maximum of 5)
#' @examples
#  rundEploid(AR, long, NumMCMC=3000, exportPostProb, recomb=0.0, k=5)
#' @export
rundEploid<-function(AR, long, NumMCMC=800, exportPostProb=TRUE, recomb=0.0, k=5){

  EmpopLong<-dplyr::ungroup(long) # ungroup long dataframe

  #identifying positions with infinite sites violations
  AR1<-base::subset(AR, IsRef=="N") #remove reference alleles
  V<-dplyr::bind_rows(AR1, EmpopLong)#bind AR to long table by appending rows

  V%>%dplyr::group_by(POS)%>% # group by position/site
    dplyr::mutate(UniqueAltAll=length(unique(Nuc)))-> V # calculate the number of unique alleles for each site
  V%>%dplyr::filter(UniqueAltAll >=2)%>%dplyr::select(POS)-> ISV # select only those sites that are represented by one allele

  V%>%dplyr::filter(UniqueAltAll <=1)-> V2 #select the remaining sites to make AltRef and long dataframes


  #convert remaining sites into the dEploid AltRef format
  V3<-dplyr::distinct(V2,unique(POS), .keep_all=TRUE)
  V3[-c(1,5,6)] %>% dplyr::mutate(CHROM = "CHRM")->AR2  #parse out updated AltRef table and add CHROM column
  #AR2<-as_tibble(AR2)%>%select(-POS)
  AR2<-stats::na.omit(AR2)#keep only sites that are in the mixture
  AR3<-AR2[,c(5, 4, 1, 2, 3)]%>% dplyr::rename(POS = "unique(POS)") #reorder columns and sort Pos column
  #create columns for Alt and Ref counts
  AR4<-AR3%>%dplyr::mutate(RefCount =dplyr::if_else(IsRef=="Y", Count, 0), AltCount =dplyr::if_else(IsRef == "N", Count, 0))
  AltRef <- dplyr::select(AR4, -c(Nuc,Count, IsRef)) # remove unnecessary columns

  #convert remaining sites into long format
  dplyr::select(V2, -c(Count, IsRef, UniqueAltAll))%>% dplyr::select(Sample.name, POS, Nuc)%>%dplyr::filter(Sample.name!= is.na(Sample.name)) -> EmpopLong #parse out updated long table

  #Create panel and plaf dataframes
  Plafpanel<-createPlafPan(EmpopLong)# run function that creates panel and plaf tables without candidate haps
  Pan<-Plafpanel[[2]] # split Pan from Plafpanel
  Plaf<-Plafpanel[[1]] # split Plaf from Plafpanel


  # Panel<-rename(Panel, CHROM="CHROM.x")
  # Panel<-select(Panel, -CHROM.y, -RefCount, -AltCount, -AltAllele)

  #Create vectors of positions for the harmonize functions
  p<-dplyr::tibble(POS=as.numeric(Pan$POS))# create a table to accomodate just the POS column from panel
  q<-dplyr::tibble(POS=as.numeric(AltRef$POS))# create a table to accomodate just t he POS column from AltRef   #

  #harmonize Panel, Plaf and AltRef so that they all have same number of sites
  Pan<-harmonizePa(Pan, q) # harmonize panel
  Plaf<-harmonizePl(Plaf, q) # harmonize plaf
  AltRef1<-harmonizeAR(AltRef, p) # harmonize AltRef

  #write dataframes to file
  readr::write_tsv(Pan,"Pan.tsv")#export data.frame as panelfile.tsv
  readr::write_tsv(Plaf, "Plaf.tsv") #export data.frame as panelfile.tsv
  readr::write_tsv(AltRef1[,c("CHROM", "POS", "AltCount")],"Alt.tsv" )#write position and count of alternate alleles to tsv
  readr::write_tsv(AltRef1[,c("CHROM", "POS", "RefCount")], "Ref.tsv")#write postion and count of reference alleles to tsv

  #assigning objects to dEploid input files
  Pan<-"Pan.tsv" #assinging the panel file
  Plaf<-"Plaf.tsv" #assigning the plaf file
  Alt<-"Alt.tsv" #assigning the alt file
  Ref<-"Ref.tsv" # assinging the ref file



  #core dEploid function
  Com<-(base::paste("-ref", Ref, "-alt", Alt, "-plaf", Plaf, "-panel", Pan, "-nSample", NumMCMC,
                    ifelse(exportPostProb==TRUE, "-exportPostProb", ""), "-recomb", recomb, "-k", k, sep = " "))
  dEploid.run<-dEploid(Com)


  base::colnames(dEploid.run$Haps)<-AltRef1$POS # assigne column names to the dEploid matrix to associate sites to alleles

  return(base::append(dEploid.run, unique(ISV))) #return modified dEploid object with a list of positions that violate infinite sites model

}




#' @title
#' Get proportion of mixtures from dEploid output
#'
#' @description
#' This function takes the output of dEploid and returns mixture proportions for major and minor contributors
#'
#'
#' @param dEploid.run dEploid object
#'
#'@examples
#'getMixProps(dEploid.run)
#' @export
getMixProps<-function(dEploid.run){
  MixProps<-utils::tail(dEploid.run$Proportions, n=1)%>% base::sort(decreasing = TRUE)
}

#' @title
#' Get estimated haplotypes based on the dEploid output
#'
#' @description
#' This function takes the output of dEploid and an Altref file and returns estimated haplotypes for the major
#' and minor contributors and their genomic positions
#'
#' @importFrom magrittr %>%
#'
#' @param dEploid.run dEploid output
#' @param MixProps Mixture proportions (getMixProps object)
#' @param AR Table of alternate and reference allele counts (tab-delimited plain text file, Empop2AltRef object)
#' @export
getdEploidHaps <-function(dEploid.run, AR) {
  Haps<-dplyr::as_tibble(dEploid.run$Haps)%>% tibble::rownames_to_column("Samples") # convert Haps from matrix to tibble
  tidyr::gather(Haps,"Key", "Value", -Samples)-> Long # go from wide to long

  Long[order(Long$Samples), c(1,2,3)]->HapsTab # reorder columns
  HapsTab1<- dplyr::filter(HapsTab, Value!=0)%>% dplyr::rename("POS"=Key) # remove all rows that have a reference allele (i.e. 0)
  HapsTab1$POS<-as.numeric(HapsTab1$POS) # convert the the vector POS to numeric to be able to join
  HapsTab2<-dplyr::left_join(HapsTab1, AR, by = "POS")%>% dplyr::select(Samples, POS, Nuc) # join dEploid haps with AltRef data to
}

