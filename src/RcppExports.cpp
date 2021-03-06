// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// C_midrqLoss
double C_midrqLoss(NumericVector b, NumericMatrix G, NumericMatrix x, NumericVector yo, NumericVector offset, int type, double tau, int n, int p, int k);
RcppExport SEXP _Qtools_C_midrqLoss(SEXP bSEXP, SEXP GSEXP, SEXP xSEXP, SEXP yoSEXP, SEXP offsetSEXP, SEXP typeSEXP, SEXP tauSEXP, SEXP nSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yo(yoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(C_midrqLoss(b, G, x, yo, offset, type, tau, n, p, k));
    return rcpp_result_gen;
END_RCPP
}
// C_midrqLoss_bc
double C_midrqLoss_bc(NumericVector b, NumericMatrix G, NumericMatrix x, NumericVector yo, NumericVector offset, int type, double tau, double lambda, int n, int p, int k);
RcppExport SEXP _Qtools_C_midrqLoss_bc(SEXP bSEXP, SEXP GSEXP, SEXP xSEXP, SEXP yoSEXP, SEXP offsetSEXP, SEXP typeSEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP nSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yo(yoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(C_midrqLoss_bc(b, G, x, yo, offset, type, tau, lambda, n, p, k));
    return rcpp_result_gen;
END_RCPP
}
// C_midrqLoss_ao
double C_midrqLoss_ao(NumericVector b, NumericMatrix G, NumericMatrix x, NumericVector yo, NumericVector offset, int type, double tau, double lambda, int n, int p, int k);
RcppExport SEXP _Qtools_C_midrqLoss_ao(SEXP bSEXP, SEXP GSEXP, SEXP xSEXP, SEXP yoSEXP, SEXP offsetSEXP, SEXP typeSEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP nSEXP, SEXP pSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yo(yoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type offset(offsetSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(C_midrqLoss_ao(b, G, x, yo, offset, type, tau, lambda, n, p, k));
    return rcpp_result_gen;
END_RCPP
}
// C_rcTest
List C_rcTest(NumericMatrix x, NumericVector psi, NumericMatrix omega, int n, int p, int B);
RcppExport SEXP _Qtools_C_rcTest(SEXP xSEXP, SEXP psiSEXP, SEXP omegaSEXP, SEXP nSEXP, SEXP pSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(C_rcTest(x, psi, omega, n, p, B));
    return rcpp_result_gen;
END_RCPP
}
// C_projfun
List C_projfun(NumericMatrix x, NumericMatrix z, NumericVector sgn, int nx, int nz, int p, int ndir);
RcppExport SEXP _Qtools_C_projfun(SEXP xSEXP, SEXP zSEXP, SEXP sgnSEXP, SEXP nxSEXP, SEXP nzSEXP, SEXP pSEXP, SEXP ndirSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sgn(sgnSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type ndir(ndirSEXP);
    rcpp_result_gen = Rcpp::wrap(C_projfun(x, z, sgn, nx, nz, p, ndir));
    return rcpp_result_gen;
END_RCPP
}
// C_phifun
List C_phifun(NumericMatrix x, NumericMatrix z, int nx, int nz, int B, int ndir, int ng, NumericVector taus, IntegerVector minn, IntegerVector maxn);
RcppExport SEXP _Qtools_C_phifun(SEXP xSEXP, SEXP zSEXP, SEXP nxSEXP, SEXP nzSEXP, SEXP BSEXP, SEXP ndirSEXP, SEXP ngSEXP, SEXP tausSEXP, SEXP minnSEXP, SEXP maxnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type nz(nzSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type ndir(ndirSEXP);
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type taus(tausSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type minn(minnSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type maxn(maxnSEXP);
    rcpp_result_gen = Rcpp::wrap(C_phifun(x, z, nx, nz, B, ndir, ng, taus, minn, maxn));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Qtools_C_midrqLoss", (DL_FUNC) &_Qtools_C_midrqLoss, 10},
    {"_Qtools_C_midrqLoss_bc", (DL_FUNC) &_Qtools_C_midrqLoss_bc, 11},
    {"_Qtools_C_midrqLoss_ao", (DL_FUNC) &_Qtools_C_midrqLoss_ao, 11},
    {"_Qtools_C_rcTest", (DL_FUNC) &_Qtools_C_rcTest, 6},
    {"_Qtools_C_projfun", (DL_FUNC) &_Qtools_C_projfun, 7},
    {"_Qtools_C_phifun", (DL_FUNC) &_Qtools_C_phifun, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_Qtools(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
