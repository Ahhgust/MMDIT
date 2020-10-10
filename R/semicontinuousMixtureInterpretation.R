
#' Estimates Fst and quantiles of Fst
#'
#' This function computes Buckleton's estimator of Fst
#' (theta)
#' and it computes an upper-bound. It's a friendly
#' wrapper for fst_buckleton
#'
#' @importFrom magrittr %>%
#'
#' @param alleles strings ; alleles (haplotypes);  1 row per chromosome sampled
#' @param populations strings ; same length as alleles (parallel array); population label associated with said haplotype
#' @param quantile ; [0,1] ; what quantile in the Fst distribution would you like?
#' @param nJack ; non-negative integer ; number of leave-one-out jackknifes taken
#' @param approximate ; boolean ; whether or not allele frequencies (approximate=TRUE) or site hetero/homozygosities (FALSE) should be used
#'
#' @export
estimateTheta <- function(alleles, populations, quantile=0.95, nJack=0, approximate=FALSE) {


  if (nJack==0) {
    # take the point-estimate ('mean')
    fst <- fst_buckleton(alleles, populations, nJack, approximate)
  } else {
    fst <- fst_buckleton(alleles, populations, nJack, approximate)
    # index 1 in fst is the point estimate (non-jackknifed)
    # remaining indexes are jackknifed. Take an upper-bound
    fst <- quantile(fst[-1], probs=quantile[[1]])
  }

  return(fst)
}


clopperHelper <- function(s, tot, quantile) {
  PropCIs::exactci(s, tot, quantile)$conf.int[2]
}

#' Estimates a naive (frequency-based) log10-likelihood of a set of observed alleles
#'
#'
#' This computes the likelihood of a vector of alleles (allelesObserved)
#' given a population of alleles (1 allele per individual, or with weights (populationWeights))
#'
#' The likelihood can be computed two ways:
#' Using the allele frequency correction of Clopper and Pearson (upper-bound taken from correctionQuantile)
#' OR
#' taking a simple point estimate (frequency).
#'
#' This computes: log10 ( n! * f(a) * f(b), ... ) ) for
#' n alleles of type: a, b, c, with frequencies f(a), f(b), ...
#' Two alleles (p,q) this gives the standard formulation of 2pq
#'
#' f(a) is the trickier part: with the simple point estimate it's computed
#' using the allele frequency (union of observed and population alleles)
#'
#' of quantile is sought, then the upper bound from the binomial is used
#' eg, if you have an allele frequency of 1/100, and upper-bound is taken
#' with the forensic standard being the upper-portion of the 95 CI
#'
#' To make for less redundancy, allele counts can be used; eg, you can have 10 A alleles and 10 B alleles
#' You can either give the function 10 As and 10 Bs (allelesInPopulation) or you can give it:
#' allelesInPopulation=c("A", "B"), populationCounts=c(10, 10)
#'
#' e.g., these are equivalent:
#' 10**estimateLog10LikelihoodNaive(c("A", "B"), rep(LETTERS[1:10], 10), correctionQuantile = 0.95
#' 10**estimateLog10LikelihoodNaive(c("A", "B"), LETTERS[1:10], populationCounts=rep(10, 10), correctionQuantile = 0.95)
#'
#' @param allelesThatExplain Vector of alleles that explain the mixture (treated as categorical variables)
#' @param allelesInPopulation Vector of alleles in the population. When populationCounts is NULL, each allele is assumed to have a weight of 1
#' @param populationCounts vector of integers ; the number of times each allele (allelesInPopulation) is found; default: all weights are 1
#' @param correctionQuantile real number; the quantile used in the Clopper and Pearson correction (binomial). If <= 0, the frequency estimate is used.
#'
#' @export
estimateLog10LikelihoodNaive <- function(allelesThatExplain, allelesInPopulation, populationCounts=NULL, correctionQuantile=0.95) {


  if (is.null(populationCounts)) {
    populationCounts <- rep(1, length(allelesInPopulation)) # if no weights, all alleles are assumed to have weight 1
  } else if (length(populationCounts) != length(allelesInPopulation)) {
    stop("The populationWeights needs to be the same length as the allelesInPopulation vector")
  } else if (any(populationCounts<1)) {
    stop("I need population counts, not population frequencies!")
  }


  nTot <- sum(populationCounts) + length(allelesThatExplain)

  # compute the allele proportion (numerator)
  # permit repeated alleles with weights
  tibble::tibble(Alleles=allelesInPopulation,
                 Weight=populationCounts
  ) %>%
    dplyr::group_by(Alleles) %>%
    dplyr::summarize(PopCount=sum(Weight), .groups='keep') %>%
    dplyr::ungroup() ->alleleCounts


  tibble::tibble(Alleles=allelesThatExplain) %>%
    dplyr::group_by(Alleles) %>%
    dplyr::summarize(ObservedCount=dplyr::n(), .groups='keep') %>%
    dplyr::left_join(
      alleleCounts,
      by="Alleles") %>%
    dplyr::mutate( #alleles not found in the population but in the same have a count of 0...
      PopCount= ifelse(is.na(PopCount), 0, PopCount),
      TotCount=PopCount+ObservedCount) -> summaries


  if (correctionQuantile<= 0) {
    # definitional: TotCount>0
    freqs <-   summaries$TotCount/nTot
  } else{
    # compute upper-bounds (Clopper and Pearson style)
    # arguments can be seen for using PopCount
    sapply( summaries$PopCount, clopperHelper, nTot,  correctionQuantile) -> freqs
  }

  # avoid underflow (take logs); return the log-likelihood
  # this is the log-equivalent form of equation 2 in Ge et al. 2010
  # that is: log10 ( n! * freq(a) * freq(b), ... ) ) for
  # n alleles of type: a, b, c...
  # eg, for two alleles p and q, this gives 2pq
  sum(
    log10(
      freqs
    )
  ) + log10(
    factorial(length(allelesThatExplain))
  )
}

#' Theta-corrected likelihood computation
#'
#' This computes the theta corrected likelihood of a set of alleles hypothesized (unknowns)
#' given a set of observed haplotypes (allelesObserved) as well as a sample of alleles from the relevant population
#' and an estimate of co-ancestry.
#'
#' This estimator is equation 3 in Ge et al. 2010 (doi.org/10.1016/j.legalmed.2010.02.003)
#'
#' at a high level this computes the likelihood of a set of alleles hypothesized to explain some mixture
#' Some are there by the hypothesis (a victim -> contributing allelesObserved)
#' and others come from a database (unknowns)
#' as well as a sample of alleles from the relevant population (allelesInPopulation)
#'
#' This estimates the allele frequencies (point estimate from known + unknown + population combined)
#' and weighs the likelihood of the alleles observed based on whether or not the unknowns and the knowns are identical
#' (x-vector).
#' A theta correction is also employed, where theta is an estimate of Fst (fst_buckleton)
#'
#' It takes in a set of allelesObserved to be in the mixture (according to the hypothesis)
#'
#' @param allelesObserved vector of alleles observed (empirically)  (not database, may be empty)
#' @param allelesThatExplain vector of alleles that explain the mixture (unknown haplotypes only)
#' @param theta theta-correction (Fst, taken from fst_buckleton)
#' @param allelesInPopulation a vector of alleles sampled from the population (also from database, but may in practice come from a different sub-population)
#' @param populationCounts a vector of integer weights (parallel to allelesInPopulation)
#'
#' @export
estimateLog10LikelihoodTheta <- function(allelesObserved, allelesThatExplain,theta, allelesInPopulation, populationCounts=NULL) {

  nUnknown <- length(allelesThatExplain)
  if (nUnknown < 1) {
    stop("I need at least 1 haplotype that explains the evidence")
  }

  #not the prettiest answer, but this expands the haplotypes
  if (! is.null(populationCounts)) {
    allelesInPopulation <- rep(allelesInPopulation, populationCounts)
  }

  # estimate the allele frequencies...
  # all alleles in one pot.
  alleles <- c(allelesObserved, allelesThatExplain, allelesInPopulation, recursive=TRUE)
  nAlleles <- length(alleles)

  tibble::tibble(Alleles=alleles) %>%
    dplyr::group_by(Alleles) %>%
    dplyr::summarize(freq=dplyr::n()/nAlleles, .groups='keep') %>%
    dplyr::ungroup() -> alleleFreqs

  tibble::tibble(Alleles=allelesThatExplain) %>%
    dplyr::left_join(alleleFreqs,
                     by="Alleles") %>%
    dplyr::pull(freq) -> Hi

  xvec <- as.integer( allelesThatExplain %in% allelesObserved)

  ivec <- 1:nUnknown
  k <- length(allelesObserved)

  # xi(theta) + (1-theta)*Pr(Hi)
  #    /
  # 1 + (i + k -2)theta

  ( (xvec*theta) + (1-theta)*Hi ) /

    (1+ (ivec + k -2)*(theta)) -> thetaSeries


  # sum of logs == log of products
  sum( log10(thetaSeries) ) + log10(factorial(nUnknown))
}




doLR <- function(haplotypes, populations, observed, variantGraph, numIndividuals, theta=NULL) {

  tibble::tibble(hapsDatabase=haplotypes,
                 Pops=populations) -> tib

  dplyr::group_by(tib, .data$Pops) %>%
    dplyr::mutate(
      trav=traverseSequencesGraph(variantGraph,
                                  unique( c(.data$hapsDatabase, observed, recursive=TRUE) ),
                                  0),

   #   explainingIndividuals=findExplainingIndividuals(variantGraph,
    #                                                  trav,
     #                                                 numIndividuals, 0),.groups='keep'
      ) -> mixInterpret


}


testLR <- function() {

  vgraph <- makeVariantGraph("ACATGA", c(0,0), c(4,4), c("","ACAT"))

  haplotypes <-  c("AGATGA","ACATGA","GA","GGATGA")
  populations <- rep("A", 4)

  doLR(haplotypes, populations, c(), vgraph, 2) -> foo

}



#' 2 person mixture simulations
#'
#' This simulates 2-person mixtures from the populations (pops) specified
#' and it creates a data frame with summary statistics on said mixtures.
#' RMNE type statistics are returned, as well as RMP/LR statistics (or the basis thereof)
#'
#' @param db database handle
#' @param pops population groups to use.
#' @param seed sets the seed in the random number generator
#' @param nMixes the number of simulations to do
#'
#'
#' @export
twopersonMix <- function(db, pops=c("AM", "EU"), seed=1,   nMixes=1000) {

  getMitoGenomes(db, pop=pops) -> genomes
  genomes$sampleid <- as.character(genomes$sampleid)


  genomes %>%
    dplyr::group_by(sequence) %>%
    dplyr::count() %>% dplyr::ungroup() -> genCount

  # decorate each genome with it's rarity (count)
  genomes %>%
    dplyr::left_join(genCount,
              by="sequence") -> genomes

  diffs <- getSeqdiffs(db, pop=pops, getPopulation=TRUE)

  diffs$event <- factor(diffs$event, levels=c("X", "D", "I"))

  getMtgenomeSequence(db, double=FALSE) -> rcrs

  peeps <- unique(genomes$sampleid)

  set.seed(seed)

  tibble::tibble(
    P1=sample(peeps, nMixes),
    P2=sample(peeps, nMixes),
    Ndun=-1,
    NMatch=-1,
    NExplain=-1
    ) %>%
    dplyr::filter(P1 != P2) -> peepPairs

  nMixes <- nrow(peepPairs) # adjust number of rows...

  for(i in 1:nMixes) {

    pairy <- dplyr::filter(diffs, sampleid == peepPairs$P1[[i]] | sampleid == peepPairs$P2[[i]])

    posits <- sort( unique(pairy$position) )

    tibble::tibble(
      position=posits,
      RefAllele=stringr::str_sub(rcrs, posits, posits)) -> refAlleles

    dplyr::inner_join(pairy,
               refAlleles,
               by="position") %>%
      dplyr::arrange(position,event) %>%
      dplyr::group_by(position) %>%
      dplyr::summarize(NeventTypes=dplyr::n_distinct(event),
                Event=event[[1]],
                Nrow=dplyr::n(),
                Oevent= event[[ Nrow ]],
                Alleles=dplyr::case_when(
                  Nrow == 1 & Event == "I" ~ paste("", basecall[[1]], sep=","),
                  Nrow == 1 & Event == "D" ~ paste("", RefAllele[[1]], sep=","),
                  Nrow == 1 & Event == "X" ~ paste(basecall[[1]], RefAllele[[1]], sep=","),
                  Oevent == "I" & Event == "I" ~ paste( unique(basecall), collapse=","),
                  Oevent == "D" & Event == "D" ~ "",
                  Oevent == "X" & Event == "X" ~ paste( unique(basecall), collapse=","),
                  Oevent == "D" & Event == "X" ~ paste( unique(basecall), collapse=","), # still works; D is ""
                  Oevent == "X" & Event == "D" ~ paste( unique(basecall), collapse=","), # still works; D is ""
                  TRUE ~ "?"
                ),
                .groups='keep'
      ) %>%
      tidyr::separate_rows(Alleles, sep=",") %>%
      dplyr::mutate(pos0=
               ifelse(Event=="I",position, position-1)) -> foo


    if (all(foo$Alleles!="?")) {
      foo %>% dplyr::arrange(pos0, position, Alleles) -> foo
      vgraph <- makeVariantGraph(rcrs,foo$pos0, foo$position, foo$Alleles)

      travs <- traverseSequencesGraph(vgraph, genCount$sequence, 0)
      explainy <- findExplainingIndividuals(vgraph, travs, 2, 0)
      peepPairs$NExplain[[i]] <- length(explainy)

      sapply(travs, getTraversalEditDistances, simplify=TRUE) -> eds
      which( sapply(eds, function(x) length(x)>0, simplify=TRUE) ) -> whodun

      dplyr::filter(genomes, sampleid == peepPairs$P1[[i]] | sampleid == peepPairs$P2[[i]]) %>%
        dplyr::pull(n) %>% sum() -> peepPairs$NMatch[[i]]
      peepPairs$Ndun[[i]] <- length( whodun )

    }

  }

  dplyr::inner_join(peepPairs,
                    dplyr::select(genomes, sampleid, n, sequence),
             by=c("P1"="sampleid")) %>%
    dplyr::rename(P1n=n, S1=sequence
           ) %>%
    dplyr::inner_join(
      dplyr::select(genomes, sampleid, n, sequence),
               by=c("P2"="sampleid")) %>%
    dplyr::rename(P2n=n, S2=sequence) %>%
    dplyr::mutate(NIndThatExplain=P1n*P2n) -> peepPairs




  # levenshtein distance
  purrr::map2(peepPairs$S1, peepPairs$S2, function(x,y) utils::adist(x,y)[[1]]) %>% unlist() -> pairwiseDists

  peepPairs$DistBetween <- pairwiseDists

  return(peepPairs)
}


#' 3 person mixture simulations
#'
#' This simulates 3-person mixtures from the populations (pops) specified
#' and it creates a data frame with summary statistics on said mixtures.
#' RMNE type statistics are returned, as well as RMP/LR statistics (or the basis thereof)
#'
#' @param db database handle
#' @param pops population groups to use.
#' @param seed sets the seed in the random number generator
#' @param nMixes the number of simulations to do
#'
#'
#' @export
threepersonMix <- function(db, pops=c("AM", "EU"), seed=1,   nMixes=1000) {

  getMitoGenomes(db, pop=pops) -> genomes
  genomes$sampleid <- as.character(genomes$sampleid)


  genomes %>%
    dplyr::group_by(sequence) %>%
    dplyr::count() %>% dplyr::ungroup() -> genCount

  # decorate each genome with it's rarity (count)
  genomes %>%
    dplyr::left_join(genCount,
                     by="sequence") -> genomes

  diffs <- getSeqdiffs(db, pop=pops, getPopulation=TRUE)

  diffs$event <- factor(diffs$event, levels=c("X", "D", "I"))

  getMtgenomeSequence(db, double=FALSE) -> rcrs

  peeps <- unique(genomes$sampleid)

  set.seed(seed)

  tibble::tibble(
    P1=sample(peeps, nMixes),
    P2=sample(peeps, nMixes),
    P3=sample(peeps, nMixes),
    Ndun=-1,
    NMatch=-1,
    NExplain=-1
  ) %>%
    dplyr::filter(P1 != P2) %>%
    dplyr::filter(P1 != P3) %>%
    dplyr::filter(P2 != P3) -> peepPairs

  nMixes <- nrow(peepPairs) # adjust number of rows...

  for(i in 1:nMixes) {

    pairy <- dplyr::filter(diffs, sampleid == peepPairs$P1[[i]] | sampleid == peepPairs$P2[[i]] | sampleid == peepPairs$P3[[i]])

    posits <- sort( unique(pairy$position) )

    tibble::tibble(
      position=posits,
      RefAllele=stringr::str_sub(rcrs, posits, posits)) -> refAlleles

    dplyr::inner_join(pairy,
                      refAlleles,
                      by="position") %>%
      dplyr::filter(event=='X') %>%
      dplyr::arrange(position,event) %>%
      dplyr::group_by(position) %>%
      dplyr::summarize(NeventTypes=dplyr::n_distinct(event),
                       Event=event[[1]],
                       Nrow=dplyr::n(),
                       Oevent= event[[ Nrow ]],
                       OOindex =ifelse( Nrow==3, 2, 1),
                       OOevent= event[[ OOindex ]],
                       Alleles=dplyr::case_when(
                         # only 1 event in the 3 people; fill in the reference allele as needed
                         Nrow == 1 & Event == "I" ~ paste("", basecall[[1]], sep=","),
                         Nrow == 1 & Event == "D" ~ paste("", RefAllele[[1]], sep=","),
                         Nrow == 1 & Event == "X" ~ paste(basecall[[1]], RefAllele[[1]], sep=","),

                         # 2 events; 1 needs to be the reference. 2 indel cases
                         Nrow == 2 & Oevent == "I" & Event == "I" ~ paste( unique(c(basecall, "", recursive=TRUE)), collapse=","),
                         Nrow == 2 & Oevent == "D" & Event == "D" ~ paste("", RefAllele[[1]], sep=","),
                         # mismatch and/or deletion
                         Nrow == 2 & Oevent %in% c("X", "D") & Event %in% c("X", "D") ~ paste(  unique(c(basecall, "", recursive=TRUE)), collapse=","),

                         # 3 events at the same location; no reference alleles
                         Nrow == 3 & OOevent %in% c("X", "D") & Oevent %in% c("X", "D") & Event %in% c("X", "D") ~ paste( unique(basecall), collapse=","), # delete/mismatch case
                         Nrow == 3 & OOevent == "I" & Oevent == "I" & Event == "I" ~ paste( unique(basecall), collapse=","), # all insertion case.
                         TRUE ~ "?"
                       ),
                       .groups='keep'
      ) %>%
      tidyr::separate_rows(Alleles, sep=",") %>%
      dplyr::mutate(pos0=
                      ifelse(Event=="I",position, position-1)) -> foo


    if (all(foo$Alleles!="?")) {
      foo %>% dplyr::arrange(pos0, position, Alleles) -> foo
      vgraph <- makeVariantGraph(rcrs,foo$pos0, foo$position, foo$Alleles)

      travs <- traverseSequencesGraph(vgraph, genCount$sequence, 0)
      explainy <- findExplainingIndividuals(vgraph, travs, 3, 0)
      peepPairs$NExplain[[i]] <- length(explainy)

      sapply(travs, getTraversalEditDistances, simplify=TRUE) -> eds
      which( sapply(eds, function(x) length(x)>0, simplify=TRUE) ) -> whodun

      dplyr::filter(genomes, sampleid == peepPairs$P1[[i]] | sampleid == peepPairs$P2[[i]]| sampleid == peepPairs$P3[[i]]) %>%
        dplyr::pull(n) %>% sum() -> peepPairs$NMatch[[i]]
      peepPairs$Ndun[[i]] <- length( whodun )

     # genCount[whodun,] %>% dplyr::inner_join(genomes, by='sequence') %>% dplyr::pull(sampleid)

    }

  }

  dplyr::inner_join(peepPairs,
                    dplyr::select(genomes, sampleid, n, sequence),
                    by=c("P1"="sampleid")) %>%
    dplyr::rename(P1n=n, S1=sequence
    ) %>%
    dplyr::inner_join(
      dplyr::select(genomes, sampleid, n, sequence),
      by=c("P2"="sampleid")) %>%
    dplyr::rename(P2n=n, S2=sequence) %>%
    dplyr::inner_join(
      dplyr::select(genomes, sampleid, n, sequence),
      by=c("P3"="sampleid")) %>%
    dplyr::rename(P3n=n, S3=sequence) %>%
    dplyr::mutate(NIndThatExplain=P1n*P2n*P3n) -> peepPairs


  # levenshtein distance
  purrr::map2(peepPairs$S1, peepPairs$S2, function(x,y) utils::adist(x,y)[[1]]) %>% unlist() -> pairwiseDists12
  # levenshtein distance
  purrr::map2(peepPairs$S1, peepPairs$S3, function(x,y) utils::adist(x,y)[[1]]) %>% unlist() -> pairwiseDists13
  purrr::map2(peepPairs$S2, peepPairs$S3, function(x,y) utils::adist(x,y)[[1]]) %>% unlist() -> pairwiseDists23


  peepPairs$DistBetween12 <- pairwiseDists12
  peepPairs$DistBetween13 <- pairwiseDists13
  peepPairs$DistBetween23 <- pairwiseDists23


  return(peepPairs)
}

test <- function() {

  ref <- "AAAAA"

  # insertion of an A
  # and an A->T
  # adjacent to each other

  # consistent strings:
  strings <- c(
    "AATAC",
    "AAAA",
    "AAAAA")

  # as SNPs (A insertion, A->T SNP; both alleles found)
  ally <- tibble::tibble(
    Start = c(2,2,4,4,4),
    Stop=   c(3,3,5,5,5),
    Alleles=c("A", "T", "", "C", "A"))

  gr <- makeVariantGraph(ref, ally$Start, ally$Stop, ally$Alleles)

  travs <- traverseSequencesGraph(gr, strings, 0)

  (explainy <- findExplainingIndividuals(gr, travs, 3, 0))


}



