
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



//' Generate one sample from a Dirichlet distribution.
//'
//' @param alpha the Dirichlet hyperparameter
//' @param seed the RNG seed: if 0 (default), generate a time-based seed
//' @return a numeric vector
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rdirichlet_cpp(const Rcpp::NumericVector &alpha, const unsigned long int seed = 0) {

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



//' Generate from a Dirichlet distribution using the stick breaking definition (safer).
//'
//' @param n how many samples to generate
//' @param alpha the Dirichlet hyperparameter, with p entries
//' @return a numeric matrix, n*p
//' @export
// [[Rcpp::export]]
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





//' Generate a Dirichlet-Dirichlet-Gamma population (unsafe).
//'
//' Generate samples from m sources and p parameters, n sample per source.
//' The between-source alpha hyperparameter used to generate the source parameters is mandatory.
//'
//' @param n number of samples per source
//' @param m number of sources
//' @param alpha_0 between-source Gamma hyperparameter, a scalar
//' @param beta_0 between-source Gamma hyperparameter, a scalar
//' @param nu_0 between-source Dirichlet hyperparameter, a numeric vector
//' @export
//' @inheritParams rdirichlet_cpp
// [[Rcpp::export]]
RcppGSL::Matrix rdirdirgamma_cpp(
      const unsigned int &n, const unsigned int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const unsigned int seed = 0
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

//' Generate a Dirichlet-Dirichlet-Gamma population (safer).
//'
//' Generate samples from m sources and p parameters, n sample per source.
//' The between-source alpha hyperparameter used to generate the source parameters is mandatory.
//'
//' @export
//' @inheritParams rdirdirgamma_cpp
// [[Rcpp::export]]
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



//' Compute column-wise sd
//'
//' @param mtx a nxp matrix
//' @return a 1xp vector
// [[Rcpp::export]]
Rcpp::NumericVector colsd(const Rcpp::NumericMatrix &mtx) {
   const unsigned int p = mtx.ncol();

   Rcpp::NumericVector sd_vec(p);

   for (unsigned int k = 0; k < p; ++k) {
      sd_vec[k] = Rcpp::sd(mtx(_, k));
   }
   return(sd_vec);
}



//' Compute the Minkowski norm of a vector
//'
//' Compute the Minkowski norm of a vector.
//' $p$ can range from 1 to infinity.
//'
//' @param v a vector
//' @param p exponent of the Minkowski norm (from 1 to Inf)
//' @return a double
// [[Rcpp::export]]
double norm_minkowski(const Rcpp::NumericVector &v, const double p = 2) {

   double r;

   if (Rcpp::traits::is_infinite<REALSXP>(p)) {
      r = Rcpp::max(v);
   } else {
      r = Rcpp::sum(Rcpp::pow(v, p));
   }
   return(r);
}


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

//' Perform ABC sampling and distance calculation using the stick breaking procedure.
//'
//' Samples from Dirichlet using [rdirdirgamma_beta_cpp()].
//'
//' Summary statistics between datasets:
//'
//' - mean
//' - standard deviation
//'
//' @param n_sample passed to [rdirdirgamma_beta_cpp()]
//' @param m_sample passed to [rdirdirgamma_beta_cpp()]
//' @export
//' @return a reps*2 matrix of distances between summary statistics
//' @inheritParams sample_ABC_rdirdirgamma_cpp
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_ABC_rdirdirgamma_beta_cpp(
      const unsigned int &n_sample, const unsigned int &m_sample,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const Rcpp::NumericMatrix &mtx_obs,
      const unsigned int &reps,
      const double &p_norm = 2
) {

   const unsigned int n = n_sample;
   const unsigned int m = m_sample;
   const unsigned int p = nu_0.size();

   // Allocate distances between summary statistics
   Rcpp::NumericMatrix mtx_norms(reps, 2);

   // Rcout << "Computing summary statistics" << std::endl;

   // Precompute observed summary statistics
   Rcpp::NumericVector mu_obs(p);
   Rcpp::NumericVector sd_obs(p);

   mu_obs = Rcpp::colMeans(mtx_obs);
   sd_obs  = colsd(mtx_obs);

   // Rcout << "Computed observed sd" << std::endl;

   // Generate observations
   Rcpp::NumericVector mu_gen(p);
   Rcpp::NumericVector sd_gen(p);

   Rcpp::NumericVector mu_diff(p);
   Rcpp::NumericVector sd_diff(p);

   // Rcout << "Starting generating data" << std::endl;

   for (unsigned int t = 0; t < reps; ++t) {

      if(t % 1000 == 0) Rcpp::checkUserInterrupt();

      // Rcout << "Allocating gen data" << std::endl;

      Rcpp::NumericMatrix mtx_gen(n, p);

      // Rcout << "Generating from Dirichlet..." << std::endl;

      mtx_gen = rdirdirgamma_beta_cpp(n, m, alpha_0, beta_0, nu_0);

      // Rcout << "Computing generated mean/sd" << std::endl;

      mu_gen = colMeans(mtx_gen);
      sd_gen = colsd(mtx_gen);

      // Compute distances between summary statistics
      mu_diff = abs(mu_gen - mu_obs);
      sd_diff = abs(sd_gen - sd_obs);

      mtx_norms(t, 0) = norm_minkowski(mu_diff, p_norm);
      mtx_norms(t, 1) = norm_minkowski(sd_diff, p_norm);
   }

   return(mtx_norms);

}





//' Compute distances between summary statistics.
//'
//' @param mtx_gen the generated data matrix
//' @param mtx_obs the observed data matrix
//' @param p_norm the power of the Minkowski distance
//' @export
//' @return a length-2 vector of distances between summary statistics
// [[Rcpp::export]]
Rcpp::NumericVector compute_distances_gen_obs_cpp(
      const Rcpp::NumericMatrix &mtx_gen,
      const Rcpp::NumericMatrix &mtx_obs,
      const double &p_norm = 2
) {

   const unsigned int p = mtx_gen.ncol();

   if (mtx_gen.ncol() != mtx_obs.ncol()) {
      Rcpp::stop("Error: different number of columns");
   }

   // Allocate distances between summary statistics
   Rcpp::NumericVector vec_norms(2);

   // Rcout << "Computing summary statistics" << std::endl;

   // Precompute observed summary statistics
   Rcpp::NumericVector mu_obs(p);
   Rcpp::NumericVector sd_obs(p);

   mu_obs = Rcpp::colMeans(mtx_obs);
   sd_obs  = colsd(mtx_obs);

   // Rcout << "Computed observed sd" << std::endl;

   // Generate observations
   Rcpp::NumericVector mu_gen(p);
   Rcpp::NumericVector sd_gen(p);

   mu_gen = colMeans(mtx_gen);
   sd_gen  = colsd(mtx_gen);

   Rcpp::NumericVector mu_diff(p);
   Rcpp::NumericVector sd_diff(p);

   mu_diff = abs(mu_gen - mu_obs);
   sd_diff = abs(sd_gen - sd_obs);

   // Compute distances between summary statistics
   vec_norms(0) = norm_minkowski(mu_diff, p_norm);
   vec_norms(1) = norm_minkowski(sd_diff, p_norm);

   return(vec_norms);

}

