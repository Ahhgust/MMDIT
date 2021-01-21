#include <Rcpp.h>
#include <cstdlib>
#include <unordered_map>
#include <random>
using namespace Rcpp;


#define NULL_FST -1e6
// [[Rcpp::plugins(cpp11)]]

typedef std::vector<int> IntVect;

//' strings to ints
//'
//' This is a Rcpp implementation of R's
//' as.integer(as.factor(strings)) - 1
//'
//' ie, it converts strings into integers (0-based)
//' the conversion maintains the natural ordering
//' (the first string is 0, the next new string i 1, ... )
//'
//' returns number of unique elements
//'
//' @param s (vector of strings)
//' @param out (vector of integers; this is filled out by the method; same length as s)
//' @return number of distinct/unique strings
//' @export
// [[Rcpp::export]]
int
str2int(Rcpp::StringVector s, std::vector<int> &out) {

  out.reserve( s.size() ); // alloc memory if needed

  // converts a string into an int
  std::unordered_map<std::string, int> converter;

  StringVector::iterator al = s.begin();
  int i=0; // the int value of the next new string

  std::string stringy;

  for ( ; al != s.end(); ++al) {
    stringy = static_cast<std::string>(*al);
    if (! converter.count(stringy) ) {
      converter[stringy] = i;
      out.push_back(i);
      ++i;
    } else {
      out.push_back( converter[stringy] );
    }
  }

  return i;
}

// inits 2D array
//
// Pure C
// allocates a 2D array of integers
// initializes all bits to 0
//
// This is done w/ 1 memory request (calloc)
// Memory needs to be free'd (one free)
//
// returns number of unique elements
//
// adapted from: https://www.geeksforgeeks.org/dynamically-allocate-2d-array-c/
//
// @param r (number of rows requested)
// @param c (number of columns requested)
//
// @return memory on heap. 2D array
int **
init2D(int r, int c) {

  int **arr, *ptr;
  int len = sizeof(int*)*r + sizeof(int)*c*r;
  arr = (int**)calloc(len, 1); // allocate memory and set all bits to 0

  ptr = (int*)(arr+r);
  for(int i =0 ; i < r; ++i)
    arr[i] = (ptr + c*i);

  return arr;
}


// Computes FST/theta
//
// This is private function
//
// sharingTable (2D array of alleles; )
// nAlleles (number of columns in array; also number of alleles)
// nPops (number of rows in array; also number of populations)
// infinitePop (boolean, when true, allele frequency is used; when false, number of pairwise differences used)
//
// estimate of overall FST
//
double
fst_buckleton_internal(int **sharingTable,
                       int nAlleles,
                       int nPops,
                       bool infinitePop=true) {

  int i,j, k, samsies, *itr;

  // population -> total # of alleles
  std::vector<int> popTotals( nPops, 0);

  // population -> total # of alleles
  std::vector<double> popHomozygosites( nPops, 0.); // M_i in Buckleton


  // compute marginal distribution (sum of rows in sharing table). Gives sample sizes in each population
  // (this is popTotals)

  // AND within all populations, compute the homozygosities
  // ie, the probability (proportion really) of being a homozygote
  // (but this is just a diploid notion; the probability of selecting the same allele twice)
  for(i = 0; i < nPops; ++i) {
    itr = sharingTable[i];
    samsies = 0;
    for(j = 0; j < nAlleles; ++j, ++itr) {
      popTotals[i] += *itr;

      // number of individuals that are the same (un-normalize probability/proportion)
      if (infinitePop)
        samsies += (*itr)*(*itr); // this maps into a proportion
      else
        samsies += (*itr)*(*itr-1); // this is the number of pairwise comparisons that are equal (homozygous, of a population).
    }

    // normalize into a proportion
    double divisor;
    if (infinitePop)
      divisor =popTotals[i]*popTotals[i];
    else
      divisor =popTotals[i]*(popTotals[i]-1);

    popHomozygosites[i] = samsies/divisor;
  }


  // now compute the homozygosity between populations
  double betweenPopHomo=0.;

  // all pairs of populations
  for (i =0; i < nPops; ++i) {

    for (j = i + 1; j < nPops; ++j) {
      double denom = popTotals[i] * popTotals[j]; // total # of comparisons made

      int *jAlleles = sharingTable[j];
      int *iAlleles = sharingTable[i];

      samsies=0;

      // un-normalized proportion of alleles that are shared between pops i and j
      for(k = 0; k < nAlleles; ++k, ++iAlleles, ++jAlleles) {
        samsies += (*iAlleles)*(*jAlleles);
      }


      // normalized.
      betweenPopHomo += samsies/denom;

    }
  }


  // compute Mb per Buckleton
  // note this corrects by n choose 2 because I look at all unordered pairs
  // the original formulation call for all ordered pairs (twice the homozygosity)
  // original formulation would instead correct by n(n-1)
  double Mb = betweenPopHomo / ( (nPops * (nPops - 1))/2);
  if (Mb == 1.0)
    return NULL_FST;

  double Bw = 0.; // Bw from Buckleton

  for (std::vector<double>::iterator itr = popHomozygosites.begin();
       itr != popHomozygosites.end();
       ++itr
       ) {
    Bw += (*itr-Mb) /(1.0-Mb);
  }

  Bw /= nPops;

  if (Bw > 0)
    return Bw;
  // Fst artifacts can dip below 0. Let's not have that.
  return 0.;
}


//' Estimates theta
//'
//' Estimate's the overall theta/fst as per:
//' doi: 10.1016/j.fsigen.2016.03.004
//'
//' In particular, this take in a vector of strings
//' and a vector of population labels:
//' eg:
//'
//' alleles <- c("A", "A", "G", "G")
//' pops <- c("CEU", "CEU", "YRI", "YRI")
//'
//' and it estimates Buckleton's FST
//' It returns a vector of length nJack+1
//' Index 1 in the vector is the overall FST
//' subsequent indexes are the jackknife estimates.
//' To get a upper bound on FST try:
//'
//' fst_buckleton(alleles, pops, nJack=1000, approximate=FALSE) -> fsts
//' quantile(fsts[-1], 0.99)
//'
//' for a naive estimate of 99CI FST
//'
//' In general, Buckleton's estimator is perhaps a bit simple in implementation
//' e.g., it takes simple averages over population-pairs to estimate the overall FST
//' Also, it uses the number of pairwise differences to compute homozygosity/heterozygosity
//' (n*(n-1)) style.
//' The approximate option instead uses allele frequencies (which are stated as an approximation)
//' The major distinction here is that all singletons (haplotypes/allele seen once contribute nothing to
//' within-population homozygosity) when approximate is TRUE (they contribute 0 pairwise differences)
//' This makes sense if, say, all alleles/haplotype are UNIQUE (FST-> 0)
//'
//' It make less sense when sample sizes are small (like in the example)
//' This would imply that large sample sizes are needed
//'
//' @param alleles (vector of strings, 1 allele per haploid individual)
//' @param populations (vector of strings; population labels)
//' @param nJack (number of jackknifes)
//' @param approximate (boolean; treat the population as being infinite in size)
//'
//' @return Numeric vector of length nJack+1. overall FST is at index [[1]], nJack jackknife estimates follow
//'
//' @export
// [[Rcpp::export]]
NumericVector
fst_buckleton(Rcpp::StringVector alleles, Rcpp::StringVector populations, int nJack=0, bool approximate=false) {

  // for the jackknife
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0, 1.0);



  int nSeqs = alleles.size();
  if (nSeqs != populations.size()) {
    throw Rcpp::exception("Input vector sizes do not match.");
  }
  if (nSeqs==0)
    return NA_REAL;

  // encode the alleles as integers
  std::vector<int> iAlleles;
  int nAlleles = str2int(alleles, iAlleles);
  if (nAlleles==1)
    return NA_REAL;


  if (nJack > nSeqs)
    nJack = nSeqs;

  NumericVector out(1 + nJack);

  // encode the populations as integers
  std::vector<int> iPops;
  int nPops = str2int(populations, iPops);
  if (nPops < 2) {
    throw Rcpp::exception("I need at least two populations...");
  }

  // 2D array: rows are allele counts; each column is some unique allele
  // each column is some population; [i,j] gives the number of alleles of type 'j' in population i
  int **sharingTable = init2D(nPops, nAlleles);

  IntVect::iterator pops = iPops.begin();
  IntVect::iterator pals = iAlleles.begin();

  std::vector< std::pair<int, int> > jacks( nJack );
  int i=0;

  for( ; pops != iPops.end(); ++pops, ++pals) {
     ++sharingTable[ *pops ][ *pals ];

    if (nJack) {
       // use a reservoir sample. Collect some allele calls to remove.
       if (i < nJack) { // haven't collected enough yet.
         jacks[i].first = *pops;
         jacks[i].second = *pals;
       } else {
          int idx = distribution(generator)*i;
          if (idx < nJack) { // pick an index at random. if it's in range, overstrike!
            jacks[idx].first = *pops;
            jacks[idx].second = *pals;
          }
       }
    }
    ++i;
  }

  /*
  for (int i =0; i < nPops; ++i) {
    for (int j =0; j < nAlleles; ++j) {
      Rcout << sharingTable[i][j];
    }
    Rcout << "\n";
  }
*/

  out[0] = fst_buckleton_internal(sharingTable,
                                  nAlleles,
                                  nPops,
                                  !approximate);

  if (nJack==0) {
    free(sharingTable);
    return out;
  }

  // jackknife takes a bit more resources
  int stop = nJack;
  if (i < nJack)
    stop = i;


  // cache answers; alleles that we've seen before
  // this flattens the 2D array of the sharingTable to just 1 offset
  std::vector<double> sitelut( nPops * nAlleles, -1);

  // and whether or not we've sampled a singleton in some population before
  std::vector<int> singletonlut(nPops, -1);
  // for which we need to know which sites are singletons...
  std::vector<bool> singletons(nAlleles, false);

  for (i=0; i < nPops; ++i) {
    int j;
    int sum=0;
    for (j = 0; j < nAlleles; ++j) {
      sum += sharingTable[i][j];
      if (sum>1)
        break;
    }
    // allele is found once in all populations,,,
    if (sum==1)
      singletons[j]=true;

  }


  for ( i=1; i <= stop; ++i) {

    int pop = jacks[i-1].first;
    int al = jacks[i-1].second;

    int offset = pop*nAlleles + al; // the 1D offset into the 2D table. Used in the LUT
    if (sitelut[offset] >= 0) {
      out[i] = sitelut[offset];
    }

    // are we removing a singleton?
    if (singletons[al] && singletonlut[pop] >= 0.0) {
      out[i] = singletonlut[pop];
      continue;
    }

    --sharingTable[pop][al]; // withhold an allele call
    out[i] = fst_buckleton_internal(sharingTable, // recompute FST
                                    nAlleles,
                                    nPops,
                                    !approximate);
    ++sharingTable[pop][al]; // and then add it back in

    // make the return value NA if the jackknife induced a monomorphism
    if (out[i] == NULL_FST)
      out[i] = NA_REAL;

    // and update the LUTs
    //
    // in case we sample this same allele in the same population we already know the answer now!
    sitelut[offset] = out[i];

    // and did we remove a singleton? (allele found once in ALL populations)
    if (singletons[al]) // now we know the answer for removing ALL singletons in this population
      singletonlut[pop] = out[i];
  }

  free(sharingTable);
  return out;
}


