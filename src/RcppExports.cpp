// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_ABC_rdirdirgamma_beta_cpp
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(const unsigned int& n_sample, const unsigned int& m_sample, const double& alpha_0, const double& beta_0, const Rcpp::NumericVector& nu_0, const Rcpp::NumericMatrix& mtx_obs, const unsigned int& reps, const double& p_norm);
RcppExport SEXP _rdirdirgamma_sample_ABC_rdirdirgamma_beta_cpp(SEXP n_sampleSEXP, SEXP m_sampleSEXP, SEXP alpha_0SEXP, SEXP beta_0SEXP, SEXP nu_0SEXP, SEXP mtx_obsSEXP, SEXP repsSEXP, SEXP p_normSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n_sample(n_sampleSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type m_sample(m_sampleSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_0(alpha_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mtx_obs(mtx_obsSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< const double& >::type p_norm(p_normSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_ABC_rdirdirgamma_beta_cpp(n_sample, m_sample, alpha_0, beta_0, nu_0, mtx_obs, reps, p_norm));
    return rcpp_result_gen;
END_RCPP
}
// compute_distances_gen_obs_cpp
Rcpp::NumericVector compute_distances_gen_obs_cpp(const Rcpp::NumericMatrix& mtx_gen, const Rcpp::NumericMatrix& mtx_obs, const double& p_norm, const bool use_optimized_summary);
RcppExport SEXP _rdirdirgamma_compute_distances_gen_obs_cpp(SEXP mtx_genSEXP, SEXP mtx_obsSEXP, SEXP p_normSEXP, SEXP use_optimized_summarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mtx_gen(mtx_genSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mtx_obs(mtx_obsSEXP);
    Rcpp::traits::input_parameter< const double& >::type p_norm(p_normSEXP);
    Rcpp::traits::input_parameter< const bool >::type use_optimized_summary(use_optimized_summarySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm, use_optimized_summary));
    return rcpp_result_gen;
END_RCPP
}
// get_number_summary_statistics
unsigned int get_number_summary_statistics(bool use_optimized_summary);
RcppExport SEXP _rdirdirgamma_get_number_summary_statistics(SEXP use_optimized_summarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type use_optimized_summary(use_optimized_summarySEXP);
    rcpp_result_gen = Rcpp::wrap(get_number_summary_statistics(use_optimized_summary));
    return rcpp_result_gen;
END_RCPP
}
// get_optimized_summary_statistics_cpp
Rcpp::NumericMatrix get_optimized_summary_statistics_cpp(const Rcpp::NumericMatrix& mtx);
RcppExport SEXP _rdirdirgamma_get_optimized_summary_statistics_cpp(SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(get_optimized_summary_statistics_cpp(mtx));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_cpp
Rcpp::NumericVector rdirichlet_cpp(const Rcpp::NumericVector& alpha, const unsigned long int seed);
RcppExport SEXP _rdirdirgamma_rdirichlet_cpp(SEXP alphaSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const unsigned long int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_cpp(alpha, seed));
    return rcpp_result_gen;
END_RCPP
}
// rdirichlet_beta_cpp
Rcpp::NumericMatrix rdirichlet_beta_cpp(const unsigned int n, Rcpp::NumericVector alpha);
RcppExport SEXP _rdirdirgamma_rdirichlet_beta_cpp(SEXP nSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirichlet_beta_cpp(n, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rdirdirgamma_cpp
RcppGSL::Matrix rdirdirgamma_cpp(const unsigned int& n, const unsigned int& m, const double& alpha_0, const double& beta_0, const Rcpp::NumericVector& nu_0, const unsigned int seed);
RcppExport SEXP _rdirdirgamma_rdirdirgamma_cpp(SEXP nSEXP, SEXP mSEXP, SEXP alpha_0SEXP, SEXP beta_0SEXP, SEXP nu_0SEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_0(alpha_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type nu_0(nu_0SEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(rdirdirgamma_cpp(n, m, alpha_0, beta_0, nu_0, seed));
    return rcpp_result_gen;
END_RCPP
}
// rdirdirgamma_beta_cpp
Rcpp::NumericMatrix rdirdirgamma_beta_cpp(const unsigned int& n, const unsigned int& m, const double& alpha_0, const double& beta_0, const Rcpp::NumericVector nu_0);
RcppExport SEXP _rdirdirgamma_rdirdirgamma_beta_cpp(SEXP nSEXP, SEXP mSEXP, SEXP alpha_0SEXP, SEXP beta_0SEXP, SEXP nu_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_0(alpha_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type beta_0(beta_0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type nu_0(nu_0SEXP);
    rcpp_result_gen = Rcpp::wrap(rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0));
    return rcpp_result_gen;
END_RCPP
}
// colsd
Rcpp::NumericVector colsd(const Rcpp::NumericMatrix& mtx);
RcppExport SEXP _rdirdirgamma_colsd(SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(colsd(mtx));
    return rcpp_result_gen;
END_RCPP
}
// norm_minkowski
double norm_minkowski(const Rcpp::NumericVector& v, const double p);
RcppExport SEXP _rdirdirgamma_norm_minkowski(SEXP vSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type v(vSEXP);
    Rcpp::traits::input_parameter< const double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_minkowski(v, p));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rdirdirgamma_sample_ABC_rdirdirgamma_beta_cpp", (DL_FUNC) &_rdirdirgamma_sample_ABC_rdirdirgamma_beta_cpp, 8},
    {"_rdirdirgamma_compute_distances_gen_obs_cpp", (DL_FUNC) &_rdirdirgamma_compute_distances_gen_obs_cpp, 4},
    {"_rdirdirgamma_get_number_summary_statistics", (DL_FUNC) &_rdirdirgamma_get_number_summary_statistics, 1},
    {"_rdirdirgamma_get_optimized_summary_statistics_cpp", (DL_FUNC) &_rdirdirgamma_get_optimized_summary_statistics_cpp, 1},
    {"_rdirdirgamma_rdirichlet_cpp", (DL_FUNC) &_rdirdirgamma_rdirichlet_cpp, 2},
    {"_rdirdirgamma_rdirichlet_beta_cpp", (DL_FUNC) &_rdirdirgamma_rdirichlet_beta_cpp, 2},
    {"_rdirdirgamma_rdirdirgamma_cpp", (DL_FUNC) &_rdirdirgamma_rdirdirgamma_cpp, 6},
    {"_rdirdirgamma_rdirdirgamma_beta_cpp", (DL_FUNC) &_rdirdirgamma_rdirdirgamma_beta_cpp, 5},
    {"_rdirdirgamma_colsd", (DL_FUNC) &_rdirdirgamma_colsd, 1},
    {"_rdirdirgamma_norm_minkowski", (DL_FUNC) &_rdirdirgamma_norm_minkowski, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_rdirdirgamma(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
