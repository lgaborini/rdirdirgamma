
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
   Rcpp::NumericMatrix mtx_distances(reps, 2);

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

      mtx_distances(t, 0) = norm_minkowski(mu_diff, p_norm);
      mtx_distances(t, 1) = norm_minkowski(sd_diff, p_norm);

   }

   // return(mtx_distances);

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

   if (reps < 1) {
      Rcpp::stop("reps must be positive");
   }
   if (n*m < n_obs) {
      Rcpp::stop("cannot generate enough observations (needed n_obs = %i, have n_generated = m*n = %i)", n_obs, n*m);
   }

   if (p != p_obs) {
      Rcpp::stop("different number of columns (nu_0: %i, observed: %i)", p, p_obs);
   }

   // Allocate distances between summary statistics
   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);
   Rcpp::NumericMatrix mtx_distances(reps, n_summary);

   // Allocate summary statistics
   Rcpp::NumericMatrix summary_gen(n_summary, p);
   Rcpp::NumericMatrix summary_obs(n_summary, p);

   // Precompute observed summary statistics
   summary_obs = get_summary_statistics_cpp(mtx_obs, use_optimized_summary);

   for (unsigned int i_rep = 0; i_rep < reps; ++i_rep) {

      if (i_rep % 1000 == 0) Rcpp::checkUserInterrupt();

      Rcpp::NumericMatrix mtx_gen(n*m, p);

      mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

      summary_gen = get_summary_statistics_cpp(
         mtx_gen(Rcpp::Range(0, n_obs - 1), _),
         use_optimized_summary
      );

      // Compute distances between summary statistics
      for (unsigned int i_summary = 0; i_summary < n_summary; i_summary++) {

         // Compute distances between summary statistics
         mtx_distances(i_rep, i_summary) = norm_minkowski(
            summary_obs(i_summary, _) - summary_gen(i_summary, _),
            p_norm
         );
      }

   }

   return(mtx_distances);

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
      Rcpp::stop("different number of columns (generated: %i, observed: %i)", p, p_obs);
   }

   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);

   // Allocate distances between summary statistics
   Rcpp::NumericVector vec_norms(n_summary);

   // Compute observed summary statistics
   Rcpp::NumericMatrix summary_obs(n_summary, p);
   Rcpp::NumericMatrix summary_gen(n_summary, p);

   summary_obs = get_summary_statistics_cpp(mtx_obs, use_optimized_summary);
   summary_gen = get_summary_statistics_cpp(mtx_gen, use_optimized_summary);

   for (unsigned int i_summary = 0; i_summary < n_summary; i_summary++) {

      // Compute distances between summary statistics
      vec_norms[i_summary] = norm_minkowski(summary_obs(i_summary,_) - summary_gen(i_summary,_), p_norm);
   }

   return(vec_norms);

}


Rcpp::NumericMatrix get_optimized_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx) {

   arma::mat mtx_arma = Rcpp::as<arma::mat>(mtx);

   // Quantiles
   // arma::vec P = { 0.1, 0.25, 0.50, 0.75, 0.9 };
   // arma::mat S = arma::quantile(mtx_arma, P, 0);

   // mean, sd, skewness, kurtosis

   return(Rcpp::wrap(arma::join_vert(
         arma::mean(mtx_arma, 0),
         arma::stddev(mtx_arma, 0, 0),   // (n-1 denominator)
         colkurtosis(mtx_arma),
         colskewness(mtx_arma)
   )));


   // return(Rcpp::wrap(S));

}

Rcpp::NumericMatrix get_standard_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx) {

   const unsigned int p = mtx.ncol();

   Rcpp::NumericMatrix mtx_summary(2, p);

   mtx_summary(0, _) = Rcpp::colMeans(mtx);
   mtx_summary(1, _) = colsd(mtx);

   return(mtx_summary);

}

Rcpp::NumericMatrix get_summary_statistics_cpp(const Rcpp::NumericMatrix &mtx, const bool use_optimized_summary) {
   if (use_optimized_summary) {
      return(get_optimized_summary_statistics_cpp(mtx));
   }
   return(get_standard_summary_statistics_cpp(mtx));
}


unsigned int get_number_summary_statistics(bool use_optimized_summary) {
   if (use_optimized_summary) {
      return(4);
   } else {
      return(2);
   }
}




// Perform ABC sampling using the stick breaking procedure, returning acceptable samples.
arma::cube generate_acceptable_data_cpp(
      const unsigned int &n_sample,
      const unsigned int &m_sample,
      const double &alpha_0,
      const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const Rcpp::NumericVector &summarize_eps,
      const unsigned int reps,
      const unsigned int max_iter,
      const double &p_norm,
      const bool use_optimized_summary
) {

   const unsigned int n_obs = mtx_obs.nrow();
   const unsigned int n = n_sample;
   const unsigned int m = m_sample;
   const unsigned int p = nu_0.size();
   const unsigned int p_obs = mtx_obs.ncol();

   if (n*m < n_obs) {
      Rcpp::stop("cannot generate enough observations (supplied n*m: %i, needed n_obs: %i)", n*m, n_obs);
   }

   if (p != p_obs) {
      Rcpp::stop("different number of columns (supplied nu_0: %i, needed: %i)", p, p_obs);
   }

   // Allocate results
   // (row, column, slice)
   arma::cube cube_samples(reps, n_obs, p);

   // Allocate distances between summary statistics
   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);
   Rcpp::NumericVector vec_distances(n_summary);

   if (summarize_eps.size() != n_summary) {
      Rcpp::stop("summarize_eps must match the number of summary statistics (supplied: %i, needed: %i)", summarize_eps.size(), n_summary);
   }

   // Generated data
   Rcpp::NumericMatrix mtx_gen(n*m, p);

   // Summary statistics
   Rcpp::NumericMatrix summary_obs(n_summary, p);
   Rcpp::NumericMatrix summary_gen(n_summary, p);

   summary_obs = get_summary_statistics_cpp(mtx_obs, use_optimized_summary);

   for (unsigned int i_rep = 0; i_rep < reps; ++i_rep) {

      bool success = false;

      for (unsigned int t = 0; t < max_iter; ++t) {

         if (t % 1000 == 0) Rcpp::checkUserInterrupt();

         mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

         summary_gen = get_summary_statistics_cpp(mtx_gen(Rcpp::Range(0, n_obs - 1), _), use_optimized_summary);

         // Allocate distances between summary statistics
         // vec_distances = compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm, use_optimized_summary);

         for (unsigned int i_summary = 0; i_summary < n_summary; i_summary++) {

            // Compute distances between summary statistics
            vec_distances[i_summary] = norm_minkowski(
               summary_obs(i_summary,_) - summary_gen(i_summary,_),
               p_norm
            );
         }


         if (is_true(all(vec_distances < summarize_eps))) {

            success = true;
            arma::mat mtx_gen_arma(n_obs, p);
            mtx_gen_arma =  mtx_gen(Rcpp::Range(0, n_obs - 1), _);

            cube_samples.row(i_rep) = mtx_gen_arma;
            break;
         }
      }

      if (success) {
         continue;
      } else {
         Rcpp::stop("Maximum number of iterations exceeded (%d).", max_iter);
      }
   }
   return(cube_samples);

}


// Perform ABC sampling using the stick breaking procedure, returning the acceptance ratio.
Rcpp::List compute_ABC_cpp(
      const unsigned int &n_sample,
      const unsigned int &m_sample,
      const double &alpha_0,
      const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const Rcpp::NumericVector &summarize_eps,
      const unsigned int reps,
      const double &p_norm,
      const bool use_optimized_summary,
      const bool return_distances
) {

   const unsigned int n_obs = mtx_obs.nrow();
   const unsigned int n = n_sample;
   const unsigned int m = m_sample;
   const unsigned int p = nu_0.size();
   const unsigned int p_obs = mtx_obs.ncol();

   if (reps < 1) {
      Rcpp::stop("reps must be positive");
   }
   if (n*m < n_obs) {
      Rcpp::stop("cannot generate enough observations (supplied n*m: %i, needed n_obs: %i)", n*m, n_obs);
   }

   if (p != p_obs) {
      Rcpp::stop("different number of columns (supplied nu_0: %i, needed: %i)", p, p_obs);
   }

   const unsigned int n_summary = get_number_summary_statistics(use_optimized_summary);

   if (summarize_eps.size() != n_summary) {
      Rcpp::stop("summarize_eps must match the number of summary statistics (supplied: %i, needed: %i)", summarize_eps.size(), n_summary);
   }

   // Allocate results
   // (row, column, slice)
   unsigned int n_accepted = 0;

   // Allocate distances between summary statistics
   Rcpp::NumericMatrix mtx_distances(reps, n_summary);

   // Generated data
   Rcpp::NumericMatrix mtx_gen(n*m, p);

   // Summary statistics
   Rcpp::NumericMatrix summary_obs(n_summary, p);
   Rcpp::NumericMatrix summary_gen(n_summary, p);

   summary_obs = get_summary_statistics_cpp(mtx_obs, use_optimized_summary);

   for (unsigned int t = 0; t < reps; ++t) {

      if (t % 1000 == 0) Rcpp::checkUserInterrupt();

      mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

      summary_gen = get_summary_statistics_cpp(mtx_gen(Rcpp::Range(0, n_obs - 1), _), use_optimized_summary);

      // Allocate distances between summary statistics
      // vec_distances = compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm, use_optimized_summary);
      mtx_distances(t,_) = compute_distances_gen_obs_cpp(mtx_gen, mtx_obs, p_norm, use_optimized_summary);

      for (unsigned int i_summary = 0; i_summary < n_summary; i_summary++) {

         // Compute distances between summary statistics
         mtx_distances(t, i_summary) = norm_minkowski(summary_obs(i_summary, _) - summary_gen(i_summary, _), p_norm);
      }


      if (is_true(all(mtx_distances(t, _) < summarize_eps))) {
         n_accepted++;
      }
   }

   // Return

   Rcpp::List l;
   l["n_accepted"] = n_accepted;
   l["accept_ratio"] = (double) n_accepted / reps;

   if (return_distances) {
      l["d_ABC"] = mtx_distances;
   }
   return (l);
}



