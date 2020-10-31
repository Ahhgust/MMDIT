
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


#' Preprocess mitochondrial genomes
#'
#' This takes the output of: getMitoGenomes
#' and indexes it. In particular, it extracts the unique haplotypes,
#' adds the haplotype ID to the genomes themselves, and it returns a list of
#' exactly 2 items.
#' The first is a dataframe with the genomes data frame
#' with an n column added (# of times the haplotype was seen) as well as a SeqID (unique ID)
#' the second is a dataframe with the haplotype sequence, n (# of times seen in database, total),
#' and the same SeqID
#' @importFrom magrittr %>%
#'
#' @param genomes a data frame from MMDIT::getMitoGenomes
#' @param knownHaps a vector of strings (full mito-genomes; from hypothesized contributors)
#' @export
preprocessMitoGenomes <- function(genomes, knownHaps=c()) {


  # make up a quasi-population whose name is "K"
  # for the "known" haplotypes
  if (length(knownHaps)>0) {
    k <- unique(knownHaps)
    tibble::tibble(
      sampleid= paste0("K", 1:length(k)),
      pop="K",
      sequence=k) %>%
      dplyr::bind_rows(genomes, .) -> genomes
  }


  # unique genomes.
  genomes %>%
    dplyr::group_by(sequence) %>%
    dplyr::count() %>%
    dplyr::ungroup() -> genCount

  genCount$SeqID <- 1:nrow(genCount) # create a unique index for each sequence...


  # decorate each genome with it's rarity (count)
  # as well as a unique integer (unique to the haplotype, not the individual)
  genomes %>%
    dplyr::left_join(genCount,
                     by="sequence") %>% dplyr::ungroup() -> genomes

  return(list(genomes, genCount))
}


#' Semicontinuous LRs
#'
#' This is an omnibus wrapper for semicontinuous likelihood estimation.
#' It implements the method of:
#' Ge, Jianye, Bruce Budowle, and Ranajit Chakraborty. "Comments on" Interpreting Y chromosome STR haplotype mixture"." Legal Medicine 13.1 (2011): 52-53.
#' as applied to variant graphs (citation coming)
#'
#' The short of it, this creates a variant graph (makeVariantGraph, using pos0, pos1 and alleles)
#' and it takes genomes from the database (genomes, which is stratified by population, genCount is every unique haplotype, regardless of population)
#' and it appends a possibly empty set of known haplotypes (knownHaps) to the set of every unique database-derived haplotype
#'
#' Then every way of explaining the mixture is computed (at the level of every known haplotype).
#' The procedure is equivalent to (in the case of 2-person mixtures), taking every pair of haplotypes and computing the
#' fraction of haplotypes that explain the mixture. To make things conservative the method of Clopper and Pearson (1934) is used
#' to take the ratio (number that explain / number considered) and map that into a conservative estimate of that ratio.
#'
#' The likelihood is estimated for every population, and for every subset of knowns possible. e.g., if 1 known haplotype is
#' given, then the likelihood of both the 1 known and 0 knowns is considered.
#' If 2 knowns are hypothesized, then the lr for both knowns, the first known, the second known (individually) and 0 knowns is computed.
#'
#' The RMNE is also computed; that is, it is the number of haplotypes that explain the mixture (divided by the total, adjusted
#' by Clopper and Pearson).
#'
#'
#'
#'
#' @importFrom magrittr %>%
#'
#' @param genomes the first data frame from MMDIT::preprocessMitoGenomes
#' @param genCount the second data frame from  MMDIT::preprocessMitoGenomes
#' @param rcrs character string. the mitochondrial genome sequence (whole thing)
#' @param pos0 0-based coordinate of alleles
#' @param pos1 1-based coordinate of alleles
#' @param alleles the alleles present in the interval specified
#' @param knownHaps a vector of haplotypes hypothesized to be in the mixture
#' @param nInMix integer; the number of distinct haploid sequences present in the mixture
#' @param clopperQuantile the upper-bound confidence interval as per Clopper and Pearson
#' @param tolerance should be 0. this permits fuzzy matching between the haplotypes and the mixture. 0 == no fuzz
#'
#' @export
semicontinuousWrapper <- function(genomes, genCount, rcrs, pos0, pos1, alleles, knownHaps=c(), nInMix=2, clopperQuantile=0.95, tolerance=0) {

  # combine database alleles with known alleles (if they exist)
  allSeqs <- c(genCount$sequence, knownHaps, recursive=TRUE)

  vgraph <- makeVariantGraph(rcrs,pos0, pos1, alleles)

  travs <- traverseSequencesGraph(vgraph, allSeqs, tolerance)

  # 0-based indexing (b/c in C++); gives index in genCount$sequence
  explainy <- findExplainingIndividuals(vgraph, travs, nInMix, tolerance)

  # 1-based indexing (b/c in R)
  # this gives a vector of lists; the inner lists enumerate every valid traversal
  sapply(travs, getTraversalEditDistances, simplify=TRUE) -> eds
  # and we only care (for the RMNE) which individuals can traverse the graph.
  # in principle length(x) should be 0 or 1, but if the user is sloppy with deletion encodings
  # there can be multiple paths
  which( sapply(eds, function(x) length(x)>0, simplify=TRUE) ) -> whodun

  # population totals
  totalsByPop <- genomes %>%
    dplyr::group_by(pop) %>%
    dplyr::summarize(Tot=dplyr::n(), .groups='keep') %>%
    dplyr::ungroup()

  if (length(explainy)==0) {

    # Honestly what the common case is.
    # you have some mixture and it contains some rare haplotype
      tidyr::expand_grid(pop=totalsByPop$pop,
                         NKnown=0:length(knownHaps)) %>%
      dplyr::left_join(totalsByPop, by='pop') %>%
      dplyr::mutate(WhichKnown="Not explainable",
                    NCombinationsThatExplain=0,
                    NExplain  = 0,
                    Divisor = choose(Tot, nInMix-NKnown)) %>%
      dplyr::select(
        pop, WhichKnown, NCombinationsThatExplain, NKnown, NExplain, Divisor) %>%
      dplyr::arrange(pop, NKnown) -> likelihoods


  } else {
    # convert the explaining haplotypes as indexes into long form
    # and wrt to Sequence IDs
    tibble::tibble(SeqID=explainy, Combo=1:length(explainy)) %>%
      tidyr::unnest(cols='SeqID') %>%
      dplyr::mutate(SeqID = SeqID + 1L) -> explainLong



    #sequence IDs have no population affiliation. This creates a lookup table
    # that let's us query the number of times each Sequence ID was fond in each population
    # once for each Combo (combination of haplotypes that explain the mixture),
    # and then applies the LUT to see how often each sequence ID is associated with some population
    tidyr::expand_grid(Combo=1:max(explainLong$Combo), pop=unique(genomes$pop)) %>%
      dplyr::left_join(explainLong, by='Combo') %>%
      dplyr::left_join(
            dplyr::select(genomes, pop, SeqID, sequence),
            by=c("pop", "SeqID")) %>%
          dplyr::group_by(Combo, pop, SeqID) %>%
             dplyr::summarize(Count=sum( ! is.na(sequence)), # the count is the number of times that haplotype is observed in the DATABASE
                              Known=SeqID[[1]]> nrow(genCount), # as a grouping variable there is 1 value; if it's indexed AFTER the database, it's a known
                              sequence=sequence[[1]],
                              .groups='keep'
                              ) %>%
      dplyr::ungroup() -> explainLongByPop


   dplyr::left_join(explainLongByPop,
              totalsByPop,
              by="pop") %>%
    dplyr::group_by(Combo, pop) %>%
        dplyr::summarize(NKnown=sum(Known==TRUE),
                       WhichKnown=paste( SeqID[Known] - nrow(genCount), collapse=","),
                       NExplainTuple= prod(Count[ !Known ]), # product of the allele COUNTS; restricted to those that are unknowns (from database)
                       # each mixture may be explainable by more than 1 combination of haplotypes
                       # this computes, within a particular combination, the number of combinations of individuals that can explain the
                       # haplotypes.
                       # This is computed for each set of tuples that can explain the mixture (1 row for each)
                       Divisor=choose(Tot[[1]], nInMix-NKnown[[1]]), # number of haplotype-tuples considered (from database)
                       .groups='keep'
                       ) %>%
     dplyr::ungroup() %>%
     dplyr::group_by(pop, WhichKnown) %>% # for some combination of population and hypothesized set of known contributors
       dplyr::summarize(
         NCombinationsThatExplain=dplyr::n(), # the number of distinct tuples that can explain the mixture
         NKnown=NKnown[[1]], # the number of known haplotypes in the hypothesis
         # this sums across rows (commented above)
         # the sum of Nexplain is the total number of tuples that can explain the mixture
         # the divisor is your n choose k formulation (the number of unordered pairs possible in the case of 2-person mixtures)
         NExplain=sum(NExplainTuple), #
         Divisor=Divisor[[1]],
         .groups='keep')%>%
     dplyr::ungroup() %>%
       # and give some pretty sorting...
     dplyr::arrange(pop, NKnown) -> likelihoods

  }
 # use the method of Ge, Budowle and Charkaborty
 # and compute the log-likelihood, with a clopper and pearson upper bound
 likelihoods$LogLikelihood_GBC <- log10(
   purrr::map2_dbl(likelihoods$NExplain, likelihoods$Divisor, clopperHelper, clopperQuantile)
  )

 likelihoods$NInMix <- nInMix

# basic RMNE (the substrates of which)
# how many haplotypes that are consistent with the mixture
# are found in each population.
# the unique haplotypes (by index, aka SeqID) are in whodun
  tidyr::expand_grid(
     pop=unique(genomes$pop),
     SeqID=whodun# note no +1; which returns 1-based indexing

     ) %>%
   dplyr::left_join(
     dplyr::select(genomes, pop, SeqID, sequence),
     by=c("pop", "SeqID")
     ) %>%
   dplyr::group_by(pop) %>%
   dplyr::summarize(Count=sum( ! is.na(sequence)),
                    CountExcludingKnown=sum( ! is.na(sequence) & !(sequence %in% knownHaps) ),
                    .groups='keep'
                    ) %>%
   dplyr::ungroup() %>%
     dplyr::left_join(totalsByPop,
                      by="pop")  -> rmneLongByPop


  if (length(knownHaps)) {
    # of the known haploytpes-- which are consistent with the mixture?
    # recalling that whodun is a vector of integers which are indexes in genCount
    # known haplotypes are added as subsequent integers (e.g., 2 more with values n+1 and n+2 in the case of
    # a 2-person mixture and n distinct haplotypes in genCount)
    knownsThatFit <- whodun[ whodun > nrow(genCount) ] - nrow(genCount)
    tibble::tibble(
      pop=1:length(knownHaps),
      Count = as.integer( pop %in% knownsThatFit),
      CountExcludingKnown   =0,
      Tot = 1) %>%
      dplyr::mutate(pop=paste0("Known haplotype ", pop)) %>%
      dplyr::bind_rows( rmneLongByPop ) -> rmneLongByPop
  }

  # compute the clopper and pearson upper-bound on the proportion of alleles that are consistent with the mixture...
   rmneLongByPop$LogRMNEUB <- log10(
     purrr::map2_dbl(rmneLongByPop$Count, rmneLongByPop$Tot, clopperHelper, quantile=clopperQuantile)
    )
   rmneLongByPop$LogRMNE_ExcludeKnowns_UB <- log10(
     purrr::map2_dbl(rmneLongByPop$CountExcludingKnown, rmneLongByPop$Tot, clopperHelper, quantile=clopperQuantile)
   )
   return( list(rmneLongByPop, likelihoods))
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
twopersonMix<- function(db, pops=c("EU"), seed=1,   nMixes=1000) {

  getMitoGenomes(db, pop=pops) -> genomes
  genomes$sampleid <- as.character(genomes$sampleid)

  foo <- preprocessMitoGenomes(genomes)
  genomes <- foo[[1]]
  genCount <- foo[[2]]


  diffs <- getSeqdiffs(db, pop=pops, getPopulation=TRUE)

  diffs$event <- factor(diffs$event, levels=c("X", "D", "I"))

  getMtgenomeSequence(db, double=FALSE) -> rcrs

  peeps <- dplyr::filter(genomes, pop==pops[[1]]) %>% dplyr::pull(sampleid)

  set.seed(seed)

  tibble::tibble(
    P1=sample(peeps, nMixes),
    P2=sample(peeps, nMixes),
    RMNE=-1,
    NMatch=-1,
    NExplain=-1,
    NComboExplain=-1,
    LogLikelihood=-1
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
      # semicontinuousWrapper <- function(genomes, genCount, pos0, pos1, alleles, knownHaps=c(), nInMix=2, clopperQuantile=0.95, tolerance=0) {
      inter <-semicontinuousWrapper(genomes, genCount, rcrs, foo$pos0, foo$position, foo$Alleles, knownHaps=c(), nInMix=2, clopperQuantile = 0.95, tolerance=0)
      rmneStats <- inter[[1]]
      lrStats <- inter[[2]]
      peepPairs$RMNE[[i]] <- rmneStats$LogRMNEUB[[1]]
      peepPairs$NMatch[[i]] <- rmneStats$Count[[1]]
      peepPairs$NExplain[[i]] <- lrStats$NExplain[[1]]
      peepPairs$NComboExplain[[i]] <- lrStats$NCombinationsThatExplain[[1]]
      peepPairs$LogLikelihood[[i]] <- lrStats$LogLikelihood_GBC[[1]]
    }

  }

  return(peepPairs)
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
twopersonMixOriginal <- function(db, pops=c("EU", "AM"), seed=1,   nMixes=1000) {

  getMitoGenomes(db, pop=pops) -> genomes
  genomes$sampleid <- as.character(genomes$sampleid)

  foo <- preprocessMitoGenomes(genomes)
  genomes <- foo[[1]]
  genCount <- foo[[2]]

  fst <- estimateTheta(genomes$sequence, genomes$pop, quantile=0.99, nJack=1000)

  diffs <- getSeqdiffs(db, pop=pops, getPopulation=TRUE)

  diffs$event <- factor(diffs$event, levels=c("X", "D", "I"))

  getMtgenomeSequence(db, double=FALSE) -> rcrs

  peeps <- dplyr::filter(genomes, pop==pops[[1]]) %>% dplyr::pull(sampleid)

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

      # convert the explaining haplotypes into individual IDs
      tibble::tibble(SeqID=explainy, Combo=1:length(explainy)) %>%
        tidyr::unnest(cols='SeqID') %>%
        dplyr::mutate(SeqID = SeqID + 1L) %>% # 0- to 1-based indexing conversion
        dplyr::left_join(genomes, by="SeqID")

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

  getMitoGenomes(db, pop=pops, ignoreIndels=FALSE) -> genomes
  genomes$sampleid <- as.character(genomes$sampleid)


  genomes %>%
    dplyr::group_by(sequence) %>%
    dplyr::count() %>% dplyr::ungroup() -> genCount

  # decorate each genome with it's rarity (count)
  genomes %>%
    dplyr::left_join(genCount,
                     by="sequence") -> genomes

  diffs <- getSeqdiffs(db, pop=pops, getPopulation=TRUE, ignoreIndels=FALSE)

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
                         Nrow == 2 & Oevent %in% c("X", "D") & Event %in% c("X", "D") ~ paste(  unique(c(basecall,  RefAllele[[1]], recursive=TRUE)), collapse=","),

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


  return(dplyr::select(peepPairs, -S3))
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



