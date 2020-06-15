// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// str2int
int str2int(Rcpp::StringVector s, std::vector<int>& out);
RcppExport SEXP _MMDIT_str2int(SEXP sSEXP, SEXP outSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type out(outSEXP);
    rcpp_result_gen = Rcpp::wrap(str2int(s, out));
    return rcpp_result_gen;
END_RCPP
}
// fst_buckleton
NumericVector fst_buckleton(Rcpp::StringVector alleles, Rcpp::StringVector populations, int nJack, bool approximate);
RcppExport SEXP _MMDIT_fst_buckleton(SEXP allelesSEXP, SEXP populationsSEXP, SEXP nJackSEXP, SEXP approximateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type alleles(allelesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type populations(populationsSEXP);
    Rcpp::traits::input_parameter< int >::type nJack(nJackSEXP);
    Rcpp::traits::input_parameter< bool >::type approximate(approximateSEXP);
    rcpp_result_gen = Rcpp::wrap(fst_buckleton(alleles, populations, nJack, approximate));
    return rcpp_result_gen;
END_RCPP
}
// ksw2_gg_align
int ksw2_gg_align(std::string Tseq, std::string Qseq, Rcpp::IntegerVector opPos, Rcpp::IntegerVector ops, int sc_mch, int sc_mis, int gapo, int gape, bool extended);
RcppExport SEXP _MMDIT_ksw2_gg_align(SEXP TseqSEXP, SEXP QseqSEXP, SEXP opPosSEXP, SEXP opsSEXP, SEXP sc_mchSEXP, SEXP sc_misSEXP, SEXP gapoSEXP, SEXP gapeSEXP, SEXP extendedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type Tseq(TseqSEXP);
    Rcpp::traits::input_parameter< std::string >::type Qseq(QseqSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type opPos(opPosSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type ops(opsSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mch(sc_mchSEXP);
    Rcpp::traits::input_parameter< int >::type sc_mis(sc_misSEXP);
    Rcpp::traits::input_parameter< int >::type gapo(gapoSEXP);
    Rcpp::traits::input_parameter< int >::type gape(gapeSEXP);
    Rcpp::traits::input_parameter< bool >::type extended(extendedSEXP);
    rcpp_result_gen = Rcpp::wrap(ksw2_gg_align(Tseq, Qseq, opPos, ops, sc_mch, sc_mis, gapo, gape, extended));
    return rcpp_result_gen;
END_RCPP
}
// seqdiffs2seq
std::string seqdiffs2seq(std::string Tseq, Rcpp::IntegerVector positions, Rcpp::IntegerVector types, Rcpp::StringVector events, int initBuff);
RcppExport SEXP _MMDIT_seqdiffs2seq(SEXP TseqSEXP, SEXP positionsSEXP, SEXP typesSEXP, SEXP eventsSEXP, SEXP initBuffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type Tseq(TseqSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type positions(positionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type types(typesSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type events(eventsSEXP);
    Rcpp::traits::input_parameter< int >::type initBuff(initBuffSEXP);
    rcpp_result_gen = Rcpp::wrap(seqdiffs2seq(Tseq, positions, types, events, initBuff));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MMDIT_str2int", (DL_FUNC) &_MMDIT_str2int, 2},
    {"_MMDIT_fst_buckleton", (DL_FUNC) &_MMDIT_fst_buckleton, 4},
    {"_MMDIT_ksw2_gg_align", (DL_FUNC) &_MMDIT_ksw2_gg_align, 9},
    {"_MMDIT_seqdiffs2seq", (DL_FUNC) &_MMDIT_seqdiffs2seq, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_MMDIT(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
