#ifndef __ABC_H
#define __ABC_H

#include <RcppArmadillo.h>
// #include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;

//' Perform ABC sampling and distance calculation using the stick breaking procedure.
//'
//' Procedure:
//'
//' 1. samples from Dirichlet using [rdirdirgamma_beta_cpp()].
//' 2. computes summary statistics on datasets:
//'
//' - column-wise mean
//' - column-wise standard deviation
//'
//' or a set of column-wise quantiles
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
//' @inheritParams get_number_summary_statistics
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(
      const unsigned int &n_sample, const unsigned int &m_sample,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const unsigned int &reps,
      const double &p_norm = 2,
      const bool use_optimized_summary = false
);


//' Generate data that is accepted by ABC.
//'
//' @param mtx_obs observed data
//' @param n_sample hyperparameters that are used to generate data
//' @param m_sample hyperparameters that are used to generate data
//' @param alpha_0 hyperparameters that are used to generate data
//' @param beta_0 hyperparameters that are used to generate data
//' @param nu_0 hyperparameters that are used to generate data
//' @param summarize_eps ABC thresholds: as many as summary statistics
//' @param n_gen how many datasets are returned
//' @param max_iter how many iterations are tried
//' @param p_norm exponent of the L^p norm (can be `Inf`)
//' @return a (n x n_obs x p) array of generated data
//' @export
//' @inheritParams get_number_summary_statistics
// [[Rcpp::export]]
arma::cube generate_acceptable_data_cpp(
      const unsigned int &n_sample,
      const unsigned int &m_sample,
      const double &alpha_0,
      const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const Rcpp::NumericVector &summarize_eps,
      const unsigned int n_gen,
      const unsigned int max_iter,
      const double &p_norm,
      const bool use_optimized_summary
);

//' Perform ABC sampling using the stick breaking procedure, returning the acceptance ratio.
//'
//' @return the acceptance ratio, where 1 means that all max_iter samples were accepted.
//' @export
//' @inheritParams get_number_summary_statistics
//' @inheritParams generate_acceptable_data_cpp
// [[Rcpp::export]]
double compute_ABC_cpp(
      const unsigned int &n_sample,
      const unsigned int &m_sample,
      const double &alpha_0,
      const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const Rcpp::NumericVector &summarize_eps,
      const unsigned int max_iter,
      const double &p_norm,
      const bool use_optimized_summary
);

//' Compute distances between summary statistics.
//'
//' @param mtx_gen the generated data matrix; number of rows is free, it must have the same number of columns as `mtx_obs`
//' @param mtx_obs the observed data matrix; number of rows is free, it must have the same number of columns as `mtx_gen`
//' @param p_norm the power of the Minkowski distance (default: 2 = Euclidean)
//' @param use_optimized_summary if TRUE, use quantile matrix, else compute mean and sd vectors
//' @export
//' @return a vector of distances between summary statistics: as many entries as summary statistics
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector compute_distances_gen_obs_cpp(
      const Rcpp::NumericMatrix &mtx_gen,
      const Rcpp::NumericMatrix &mtx_obs,
      const double &p_norm = 2,
      const bool use_optimized_summary = false
);



//' Get the number of multivariate summary statistics.
//'
//' @param use_optimized_summary if TRUE, return the optimized summary statistics, else standard (mean, sd)
//' @export
//' @return an integer
// [[Rcpp::export(rng = false)]]
unsigned int get_number_summary_statistics(bool use_optimized_summary);


//' Compute optimized summary statistics.
//'
//' @export
//' @return a kxp matrix of summary statistics
//' @inheritParams get_summary_statistics_cpp
//' @inheritParams get_number_summary_statistics
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_optimized_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx);

//' Compute standard summary statistics.
//'
//' @export
//' @return a kxp matrix of summary statistics
//' @inheritParams get_summary_statistics_cpp
//' @inheritParams get_number_summary_statistics
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_standard_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx);

//' Compute summary statistics.
//'
//' @param mtx a data matrix (nxp)
//' @inheritParams get_number_summary_statistics
//' @export
//' @return a kxp matrix of summary statistics
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix get_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx, const bool use_optimized_summary);

#endif
