// https://github.com/michealhan/bbc_model/blob/8a45e0a7d87b0396b8556d9790743d592a23f899/gibbs_run.cpp

// #include <Rcpp.h>
#include <RcppGSL.h>

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <sys/time.h>

// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;

// Global RNG

// Allocate random number generator
gsl_rng *r_RNG = gsl_rng_alloc(gsl_rng_mt19937);

unsigned long int random_seed()
{
   struct timeval tv;
   gettimeofday(&tv,0);
   return (tv.tv_sec + tv.tv_usec);
}



// Generate one sample from a Dirichlet distribution.
//
// @param alpha the Dirichlet hyperparameter
// @return a numeric vector
// @export
// [[Rcpp::export]]
Rcpp::NumericVector rdirichlet_cpp(const Rcpp::NumericVector &alpha, const unsigned long int seed = 0) {

   int n = alpha.size();
   Rcpp::NumericVector results(n);

   if (seed == 0) {
      gsl_rng_set(r_RNG, random_seed());
   } else {
      gsl_rng_set(r_RNG, seed);
   }

   gsl_ran_dirichlet(r_RNG, n, alpha.begin(), results.begin());

   // Release random number generator
   // gsl_rng_free(r);

   return(results);
}

// Generate one sample from a Dirichlet distribution using the stick breaking definition.
//
// @param n how many samples to generate
// @param alpha the Dirichlet hyperparameter, with p entries
// @return a numeric matrix, n*p
// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix rdirichlet_beta_cpp(unsigned int n, Rcpp::NumericVector alpha) {

   unsigned int p = alpha.size();
   double phi;

   Rcpp::NumericMatrix r(n, p);

   for (unsigned int i = 0; i < n; i++) {
      for (unsigned int k = 0; k < p - 1; ++k) {
         phi = R::rbeta(alpha[k], sum( alpha[Range(k + 1, p - 1)] ) );
         r(i, k) = (1 - sum(r(i, _))) * phi;
      }
      r(i, p - 1) = 1 - sum(r(i, _));
   }

   return(r);
}





// Generate a Dirichlet-Dirichlet-Gamma population.

// Generate samples from m sources and p parameters, n sample per source.
// The between-source alpha hyperparameter used to generate the source parameters is mandatory.

// @param n number of samples per source
// @param m number of sources
// @param alpha_0 between-source Gamma hyperparameter, a scalar
// @param beta_0 between-source Gamma hyperparameter, a scalar
// @param nu_0 between-source alpha hyperparameter, a numeric vector.
// @param seed seed of the RNG (if 0, a random seed is used)
// @export
//
// [[Rcpp::export]]
RcppGSL::Matrix rdirdirgamma_cpp(
      const int &n, const int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const unsigned int seed = 0
   ) {

   const int p = nu_0.size();

   RcppGSL::Vector vec_gamma(m);
   RcppGSL::Matrix M(m, p);
   RcppGSL::Matrix X(m*n, p);

   const double shape = alpha_0;
   const double scale = 1/beta_0;

   if (seed == 0) {
      gsl_rng_set(r_RNG, random_seed());
   } else {
      gsl_rng_set(r_RNG, seed);
   }

   // The Gamma part
   for (int j = 0; j < m; ++j) {
      vec_gamma[j] = gsl_ran_gamma(r_RNG, shape, scale);
   }

   // The nu part
   // M <- ridirichlet_beta(n = m, a = as.numeric(nu_0))
   for (int j = 0; j < m; ++j) {

      // RcppGSL::VectorView rowview = gsl_matrix_row(M, i);
      // gsl_ran_dirichlet(r, p, nu_0.begin(), *rowview);

      gsl_vector_view rowview_source = gsl_matrix_row(M, j);
      gsl_ran_dirichlet(r_RNG, p, nu_0.begin(), (&rowview_source.vector)->data);

      // The mu*Gamma part
      // Rescale by the concentration parameter
      gsl_blas_dscal(vec_gamma[j], &rowview_source.vector);
   }

   // The mu*Gamma part
   // M_alpha <- t(t(M) %*% diag(vec_gamma))

   // The observations
   for (int j = 0; j < m; ++j) {
      // RcppGSL::VectorView rowview = gsl_matrix_row(M, i);
      // gsl_ran_dirichlet(r, p, nu_0.begin(), *rowview);

      gsl_vector_const_view rowview_source = gsl_matrix_const_row(M, j);

      for (int i = 0; i < n; ++i) {
         gsl_vector_view rowview_obs = gsl_matrix_row(X, j*n + i);
         gsl_ran_dirichlet(r_RNG, p, (&rowview_source.vector)->data, (&rowview_obs.vector)->data);
      }
   }

   // return(
   //    Rcpp::List::create(
   //       Rcpp::Named("mtx_sources") = M,
   //       Rcpp::Named("mtx_data") = X
   //    )
   // );

   return(X);
}

// Generate a Dirichlet-Dirichlet-Gamma population.

// Generate samples from m sources and p parameters, n sample per source.
// The between-source alpha hyperparameter used to generate the source parameters is mandatory.

// @param n number of samples per source
// @param m number of sources
// @param alpha_0 between-source Gamma hyperparameter, a scalar
// @param beta_0 between-source Gamma hyperparameter, a scalar
// @param nu_0 between-source alpha hyperparameter, a numeric vector.
// @export
//
// [[Rcpp::export]]
Rcpp::NumericMatrix rdirdirgamma_beta_cpp(
      const int &n, const int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0
   ) {

   const int p = nu_0.size();

   Rcpp::NumericVector vec_gamma(m);
   Rcpp::NumericMatrix M(m, p);
   Rcpp::NumericMatrix X(m*n, p);

   // Save temporary observations for a single source
   Rcpp::NumericMatrix X_single(n, p);

   // The Gamma part
   for (int j = 0; j < m; ++j) {
      vec_gamma[j] = R::rgamma(alpha_0, beta_0);
   }

   // The nu part
   // M <- ridirichlet_beta(n = m, a = as.numeric(nu_0))

   M = rdirichlet_beta_cpp(m, nu_0);

   for (int j = 0; j < m; ++j) {

      // The mu*Gamma part
      // Rescale by the concentration parameter
      M(j, _) = M(j, _) * vec_gamma[j];
   }

   // The mu*Gamma part
   // M_alpha <- t(t(M) %*% diag(vec_gamma))

   // The observations
   for (int j = 0; j < m; ++j) {
      X_single = rdirichlet_beta_cpp(n, M(j, _));
//
//       for (int k = 0; k < p; ++k) {
//          X(Range(j, j + n - 1), k) = X_single(Range(0, n - 1), k);
//       }

      for (int i = 0; i < n; ++i) {
         X(j*m + i, _) = X_single(i, _);
      }

      // X(Range(j, j + n - 1), _) = X_single;
   }

   // return(
   //    Rcpp::List::create(
   //       Rcpp::Named("mtx_sources") = M,
   //       Rcpp::Named("mtx_data") = X
   //    )
   // );

   return(X);
}



// Perform ABC sampling and distance calculation.
//
// @export
// [[Rcpp::export]]
Rcpp::NumericVector sample_ABC_rdirdirgamma_cpp_internal(
      const int &n, const int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const int &reps,
      const int &p_norm = 2,
      const unsigned int seed = 0
) {

   int p = nu_0.size();
   int n_obs = mtx_obs.nrow();

   // Allocate distances between summary statistics
   Rcpp::NumericMatrix mtx_norms(reps, 2);

   // Precompute observed summary statistics
   Rcpp::NumericVector mu_obs(p);
   Rcpp::NumericVector sd_obs(p);

   mu_obs = Rcpp::colMeans(mtx_obs);

   for (int k = 0; k < p; ++k) {
      // gsl_vector_view col = gsl_matrix_const_column(mtx_obs, k);
      // sd_obs[k] = gsl_stats_mean((&col.vector)->data, 1, n_obs);
      sd_obs[k] = Rcpp::sd(mtx_obs(_, k));
   }

   return(sd_obs);

   // Generate observations
   Rcpp::NumericVector mu_gen(p);
   Rcpp::NumericVector sd_gen(p);

   Rcpp::NumericVector mu_diff(p);
   Rcpp::NumericVector sd_diff(p);

   for (int i = 0; i < reps; ++i) {
      RcppGSL::Matrix mtx_gen(n, p);
      // Rcpp::NumericMatrix mtx_gen(n, p);

      mtx_gen = rdirdirgamma_cpp(n, m, alpha_0, beta_0, nu_0);

      for (int k = 0; k < p; ++k) {
         // mu_gen[k] = Rcpp::colMeans(mtx_gen);

         gsl_vector_const_view col = gsl_matrix_const_column(mtx_gen, k);
         mu_gen[k] = gsl_stats_mean((&col.vector)->data, 1, n_obs);
         sd_gen[k] = gsl_stats_sd((&col.vector)->data, 1, n_obs);
      }

      mu_diff = mu_gen - mu_obs;
      sd_diff = sd_gen - sd_obs;

      // return(mtx_gen);

      // gsl_vector_const_view vec_mu_diff = gsl_vector_const_view_array(mu_diff, p);
      // gsl_vector_const_view vec_sd_diff = gsl_vector_const_view_array(sd_diff, p);
      //
      // // Compute distances between summary statistics
      // mtx_norms(i, 1) = gsl_blas_dnrm2( (&vec_mu_diff.vector)->data );
      // mtx_norms(i, 2) = gsl_blas_dnrm2( (&vec_sd_diff.vector)->data );

      mtx_norms(i, 1) = sum(pow(mu_diff, p_norm));
      mtx_norms(i, 2) = sum(pow(sd_diff, p_norm));


   }

   // return(mtx_norms);



// Perform ABC sampling and distance calculation using the stick breaking procedure.
//
// @export
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(
      const int &n, const int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const int &reps,
      const int &p_norm = 2
) {

   int p = nu_0.size();

   // Allocate distances between summary statistics
   Rcpp::NumericMatrix mtx_norms(reps, 2);

   // Precompute observed summary statistics
   Rcpp::NumericVector mu_obs(p);
   Rcpp::NumericVector sd_obs(p);

   mu_obs = Rcpp::colMeans(mtx_obs);

   for (int k = 0; k < p; ++k) {
      sd_obs[k] = Rcpp::sd(mtx_obs(_, k));
   }

   Rcout << "Computed observed sd" << std::endl;

   // Generate observations
   Rcpp::NumericVector mu_gen(p);
   Rcpp::NumericVector sd_gen(p);

   Rcpp::NumericVector mu_diff(p);
   Rcpp::NumericVector sd_diff(p);

   Rcout << "Starting generating data" << std::endl;

   for (int t = 0; t < reps; ++t) {

      Rcout << "Allocating gen data" << std::endl;

      Rcpp::NumericMatrix mtx_gen(n, p);

      Rcout << "Generating from Dirichlet..." << std::endl;

      mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

      Rcout << "Computed generated mean/sd" << std::endl;

      mu_gen = colMeans(mtx_gen);
      // sd_gen = Rcpp::sd(mtx_gen);

      for (int k = 0; k < p; ++k) {
         sd_gen[k] = sd(mtx_gen(_, k));
      }

      Rcout << "Computed generated mean/sd differences" << std::endl;

      mu_diff = mu_gen - mu_obs;
      sd_diff = sd_gen - sd_obs;

      // Compute distances between summary statistics

      mtx_norms(t, 1) = sum(pow(mu_diff, p_norm));
      mtx_norms(t, 2) = sum(pow(sd_diff, p_norm));


   }

   return(mtx_norms);

}

