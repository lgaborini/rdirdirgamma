
#include <RcppArmadillo.h>
// #include <Rcpp.h>

#include <sys/time.h>

#include "rng.h"
#include "abc.h"
#include "utils.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


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

// Perform ABC sampling and distance calculation using the stick breaking procedure.
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(
      const unsigned int &n_sample, const unsigned int &m_sample,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const unsigned int &reps,
      const double &p_norm,
      const bool use_optimized_summary
) {

   const unsigned int n_obs = mtx_obs.nrow();
   const unsigned int n = n_sample;
   const unsigned int m = m_sample;
   const unsigned int p = nu_0.size();
   const unsigned int p_obs = mtx_obs.ncol();

   if (n*m < n_obs) {
      Rcpp::stop("cannot generate enough observations (needed n_obs = %i, have n_gen = %i)", n_obs, n*m);
   }

   if (p != p_obs) {
      Rcpp::stop("Error: different number of columns (nu_0: %i, observed: %i)", p, p_obs);
   }

   // Allocate distances between summary statistics
   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);
   Rcpp::NumericMatrix mtx_norms(reps, n_summary);

   if (use_optimized_summary) {

      // Quantile matrix
      Rcpp::NumericMatrix summary_obs(n_summary, p);
      Rcpp::NumericMatrix summary_gen(n_summary, p);

      summary_obs = get_optimized_summary_statistics_cpp(mtx_obs);

      for (unsigned int t = 0; t < reps; ++t) {

         if (t % 1000 == 0) Rcpp::checkUserInterrupt();

         Rcpp::NumericMatrix mtx_gen(n*m, p);

         mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

         summary_gen = get_optimized_summary_statistics_cpp(mtx_gen(Rcpp::Range(0, n_obs - 1), _));

         // Compute distances between summary statistics
         for (unsigned int i_summary = 0; i_summary < n_summary; i_summary++) {

            // Compute distances between summary statistics
            mtx_norms(t, i_summary) = norm_minkowski(summary_obs(i_summary,_) - summary_gen(i_summary,_), p_norm);
         }
      }


   } else {

      // Classic statistics: mean, sd

      // Precompute observed summary statistics
      Rcpp::NumericVector mu_obs(p);
      Rcpp::NumericVector sd_obs(p);

      mu_obs = Rcpp::colMeans(mtx_obs);
      sd_obs = colsd(mtx_obs);

      // Generate observations
      Rcpp::NumericVector mu_gen(p);
      Rcpp::NumericVector sd_gen(p);

      for (unsigned int t = 0; t < reps; ++t) {

         if (t % 1000 == 0) Rcpp::checkUserInterrupt();

         Rcpp::NumericMatrix mtx_gen(n*m, p);

         mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

         mu_gen = colMeans(mtx_gen(Rcpp::Range(0, n_obs - 1), _));
         sd_gen = colsd(mtx_gen(Rcpp::Range(0, n_obs - 1), _));

         // Compute distances between summary statistics
         mtx_norms(t, 0) = norm_minkowski(mu_gen - mu_obs, p_norm);
         mtx_norms(t, 1) = norm_minkowski(sd_gen - sd_obs, p_norm);
      }

   }


   return(mtx_norms);

}





// Compute distances between summary statistics.
Rcpp::NumericVector compute_distances_gen_obs_cpp(
      const Rcpp::NumericMatrix &mtx_gen,
      const Rcpp::NumericMatrix &mtx_obs,
      const double &p_norm,
      const bool use_optimized_summary
) {

   const unsigned int p = mtx_gen.ncol();
   const unsigned int p_obs = mtx_obs.ncol();

   if (p != p_obs) {
      Rcpp::stop("Error: different number of columns (generated: %i, observed: %i)", p, p_obs);
   }

   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);

   // Allocate distances between summary statistics
   Rcpp::NumericVector vec_norms(n_summary);

   if (use_optimized_summary) {

      // Precompute observed summary statistics
      Rcpp::NumericMatrix summary_obs(n_summary, p);
      Rcpp::NumericMatrix summary_gen(n_summary, p);

      summary_obs = get_optimized_summary_statistics_cpp(mtx_obs);
      summary_gen = get_optimized_summary_statistics_cpp(mtx_gen);

      for (unsigned int i = 0; i < n_summary; i++) {

         // Compute distances between summary statistics
         vec_norms[i] = norm_minkowski(summary_obs(i,_) - summary_gen(i,_), p_norm);
      }

   } else {


      // Precompute observed summary statistics
      Rcpp::NumericVector mu_obs(p);
      Rcpp::NumericVector sd_obs(p);

      mu_obs = Rcpp::colMeans(mtx_obs);
      sd_obs = colsd(mtx_obs);

      // Generate observations
      Rcpp::NumericVector mu_gen(p);
      Rcpp::NumericVector sd_gen(p);

      mu_gen = colMeans(mtx_gen);
      sd_gen = colsd(mtx_gen);

      // Compute distances between summary statistics
      vec_norms[0] = norm_minkowski(mu_gen - mu_obs, p_norm);
      vec_norms[1] = norm_minkowski(sd_gen - sd_obs, p_norm);


   }

   return(vec_norms);

}


Rcpp::NumericMatrix get_optimized_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx) {

   arma::mat mtx_arma = Rcpp::as<arma::mat>(mtx);

   arma::vec P = { 0.1, 0.25, 0.50, 0.75, 0.9 };
   arma::mat S = arma::quantile(mtx_arma, P, 0);
   return(Rcpp::wrap(S));

}


unsigned int get_number_summary_statistics(bool use_optimized_summary) {
   if (use_optimized_summary) {
      return(5);
   } else {
      return(2);
   }
}

