// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// criterionC
NumericVector criterionC(NumericMatrix X, NumericVector y, double sigma, NumericMatrix betas, NumericVector g_index, NumericVector grps, NumericVector K);
RcppExport SEXP _SA23204165_criterionC(SEXP XSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP betasSEXP, SEXP g_indexSEXP, SEXP grpsSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type betas(betasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type g_index(g_indexSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type grps(grpsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(criterionC(X, y, sigma, betas, g_index, grps, K));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SA23204165_criterionC", (DL_FUNC) &_SA23204165_criterionC, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SA23204165(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
