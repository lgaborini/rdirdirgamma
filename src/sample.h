#ifndef __SAMPLE_H
#define __SAMPLE_H


#include <Rcpp.h>

// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;

// Global RNG

// Allocate random number generator
// gsl_rng *r_RNG = gsl_rng_alloc(gsl_rng_mt19937);


/*

 //' Perform ABC sampling and distance calculation.
 //'
 //' Samples from Dirichlet using [rdirdirgamma_cpp()].
 //'
 //' Summary statistics between datasets:
 //'
 //' - mean
 //' - standard deviation
 //'
 //' ## RNG
 //'
 //' this function uses GSL's RNG seed, unaffected by R's RNG.
 //'
 //' @param mtx_obs the observed data matrix
 //' @param reps repetitions to average distances (default: 1)
 //' @param n_sample number of samples per source
 //' @param m_sample number of sources
 //' @param p_norm exponent of the L^p norm (can be `Inf`) (default: 2)
 //' @return a reps*2 matrix of distances between summary statistics
 //' @export
 //' @inheritParams rdirdirgamma_cpp
 // [[Rcpp::export]]
 RcppGSL::Matrix sample_ABC_rdirdirgamma_cpp(
 const unsigned int &n_sample, const unsigned int &m_sample,
 const double &alpha_0, const double &beta_0,
 const Rcpp::NumericVector &nu_0,
 const Rcpp::NumericMatrix &mtx_obs,
 const unsigned int &reps,
 const double &p_norm = 2,
 const unsigned int seed = 0
 ) {
 // too lazy
 const unsigned int n = n_sample;
 const unsigned int m = m_sample;

 const unsigned int p = nu_0.size();
 const unsigned int n_obs = mtx_obs.nrow();

 // Allocate distances between summary statistics
 Rcpp::NumericMatrix mtx_norms(reps, 2);

 // Precompute observed summary statistics
 Rcpp::NumericVector mu_obs(p);
 Rcpp::NumericVector sd_obs(p);

 mu_obs = Rcpp::colMeans(mtx_obs);

 for (unsigned int k = 0; k < p; ++k) {
 // gsl_vector_view col = gsl_matrix_const_column(mtx_obs, k);
 // sd_obs[k] = gsl_stats_mean((&col.vector)->data, 1, n_obs);
 sd_obs[k] = Rcpp::sd(mtx_obs(_, k));
 }

 Rcout << "Computed observed sd" << std::endl;

 // Generate observations
 Rcpp::NumericVector mu_gen(p);
 Rcpp::NumericVector sd_gen(p);

 Rcpp::NumericVector mu_diff(p);
 Rcpp::NumericVector sd_diff(p);

 Rcout << "Starting generating data" << std::endl;

 for (unsigned int t = 0; t < reps; ++t) {

 if(t % 1000 == 0) Rcpp::checkUserInterrupt();

 Rcout << "Allocating gen data" << std::endl;

 RcppGSL::Matrix mtx_gen(n, p);
 // Rcpp::NumericMatrix mtx_gen(n, p);

 Rcout << "Generating from Dirichlet..." << std::endl;

 mtx_gen = rdirdirgamma_cpp(n, m, alpha_0, beta_0, nu_0);

 return(mtx_gen);

 Rcout << "Computed generated mean/sd" << std::endl;

 for (unsigned int k = 0; k < p; ++k) {
 // mu_gen[k] = Rcpp::colMeans(mtx_gen);

 gsl_vector_const_view col = gsl_matrix_const_column(mtx_gen, k);
 mu_gen[k] = gsl_stats_mean((&col.vector)->data, 1, n_obs);
 sd_gen[k] = gsl_stats_sd((&col.vector)->data, 1, n_obs);
 }

 Rcout << "Computed generated mean/sd differences" << std::endl;

 mu_diff = abs(mu_gen - mu_obs);
 sd_diff = abs(sd_gen - sd_obs);

 mtx_norms(t, 0) = norm_minkowski(mu_diff, p_norm);
 mtx_norms(t, 1) = norm_minkowski(sd_diff, p_norm);

 }

 // return(mtx_norms);

 // avoid complaining
 RcppGSL::Matrix mtx_gen(n, p);
 return(mtx_gen);

 }

 */

//'
//' Procedure:
//'
//' 1. samples from Dirichlet using [rdirdirgamma_beta_cpp()].
//' 2. computes summary statistics on datasets:
//'
//' - column-wise mean
//' - column-wise standard deviation
//'
//' 3. the generated dataset is invisibly is truncated to the same amount of rows as the observed dataset.
//' 4. compute the Minkowski norms of the differences between summary statistics.
//' 5. repeat `reps` times.
//'
//' ## RNG
//'
//' this function uses R's RNG seed.
//'
//' @param mtx_obs the observed data matrix
//' @param reps number of ABC samples (default: 1)
//' @param n_sample number of samples per source
//' @param m_sample number of sources
//' @param p_norm exponent of the L^p norm (can be `Inf`) (default: `2`)
//' @export
//' @return a reps*2 matrix of distances between summary statistics
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(
      const unsigned int &n_sample, const unsigned int &m_sample,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const unsigned int &reps,
      const double &p_norm = 2
);





//' Compute distances between summary statistics.
//'
//' @param mtx_gen the generated data matrix; number of rows is free, it must have the same number of columns as `mtx_obs`
//' @param mtx_obs the observed data matrix; number of rows is free, it must have the same number of columns as `mtx_gen`
//' @param p_norm the power of the Minkowski distance (default: 2 = Euclidean)
//' @export
//' @return a length-2 vector of distances between summary statistics
// [[Rcpp::export]]
Rcpp::NumericVector compute_distances_gen_obs_cpp(
      const Rcpp::NumericMatrix &mtx_gen,
      const Rcpp::NumericMatrix &mtx_obs,
      const double &p_norm = 2
);




#endif
