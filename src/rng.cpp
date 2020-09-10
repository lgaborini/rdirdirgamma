
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
Rcpp::NumericVector rdirichlet_cpp(const Rcpp::NumericVector &alpha, const unsigned long int seed) {

   unsigned int n = alpha.size();
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

// TODO:
// check out the numerically stable Dirichlet RNG
//
// https://github.com/michealhan/bbc_model/blob/8a45e0a7d87b0396b8556d9790743d592a23f899/gibbs_run.cpp



// Generate from a Dirichlet distribution using the stick breaking definition (safer).
Rcpp::NumericMatrix rdirichlet_beta_cpp(const unsigned int n, Rcpp::NumericVector alpha) {

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





// Generate a Dirichlet-Dirichlet-Gamma population (unsafe).
RcppGSL::Matrix rdirdirgamma_cpp(
      const unsigned int &n, const unsigned int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const unsigned int seed
) {

   const unsigned int p = nu_0.size();

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
   for (unsigned int j = 0; j < m; ++j) {
      vec_gamma[j] = gsl_ran_gamma(r_RNG, shape, scale);
   }

   // The nu part
   // M <- ridirichlet_beta(n = m, a = as.numeric(nu_0))
   for (unsigned int j = 0; j < m; ++j) {

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
   for (unsigned int j = 0; j < m; ++j) {
      // RcppGSL::VectorView rowview = gsl_matrix_row(M, i);
      // gsl_ran_dirichlet(r, p, nu_0.begin(), *rowview);

      gsl_vector_const_view rowview_source = gsl_matrix_const_row(M, j);

      for (unsigned int i = 0; i < n; ++i) {
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


// Generate a Dirichlet-Dirichlet-Gamma population (safer).
Rcpp::NumericMatrix rdirdirgamma_beta_cpp(
      const unsigned int &n, const unsigned int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector nu_0
) {

   const unsigned int p = nu_0.size();

   Rcpp::NumericVector vec_gamma(m);
   Rcpp::NumericMatrix M(m, p);
   Rcpp::NumericMatrix X(m*n, p);

   const double shape = alpha_0;
   const double scale = 1/beta_0;

   // Save temporary observations for a single source
   Rcpp::NumericMatrix X_source(n, p);

   // The Gamma part
   // for (int j = 0; j < m; ++j) {
   //    vec_gamma[j] = R::rgamma(alpha_0, beta_0);
   // }
   vec_gamma = Rcpp::rgamma(m, shape, scale);

   // The nu part
   // M <- ridirichlet_beta(n = m, a = as.numeric(nu_0))

   M = rdirichlet_beta_cpp(m, nu_0);

   for (unsigned int j = 0; j < m; ++j) {

      // The mu*Gamma part
      // Rescale by the concentration parameter
      M(j, _) = M(j, _) * vec_gamma[j];
   }

   // The mu*Gamma part
   // M_alpha <- t(t(M) %*% diag(vec_gamma))

   // Generate the observations

   // j-th source
   for (unsigned int j = 0; j < m; ++j) {

      // Cannot do this
      // X(Range(j, j + n - 1), _) = X_source;

      X_source = rdirichlet_beta_cpp(n, M(j, _));

      for (unsigned int i = 0; i < n; ++i) {
         X(j*n + i, _) = X_source(i, _);
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

