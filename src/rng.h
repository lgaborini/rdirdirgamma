#ifndef __RNG_H
#define __RNG_H

// #include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>

// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;

// Setup RNG for GSL
extern gsl_rng *r_RNG;

unsigned long int random_seed();


//' Generate one sample from a Dirichlet distribution.
//'
//' ## RNG
//'
//' this function uses GSL's RNG seed, unaffected by R's RNG.
//'
//' @param alpha the Dirichlet hyperparameter
//' @param seed the RNG seed: if 0 (default), generate a time-based seed
//' @return a numeric vector
//' @export
//' @family RNG functions
// [[Rcpp::export]]
Rcpp::NumericVector rdirichlet_cpp(const Rcpp::NumericVector &alpha, const unsigned long int seed = 0);

//' Generate from a Dirichlet distribution using the stick breaking definition (safer).
//'
//' ## RNG
//'
//' this function uses R's RNG seed.
//'
//' @param n how many samples to generate
//' @param alpha the Dirichlet hyperparameter, with p entries
//' @return a numeric matrix, n*p
//' @export
//' @family RNG functions
// [[Rcpp::export]]
Rcpp::NumericMatrix rdirichlet_beta_cpp(const unsigned int n, Rcpp::NumericVector alpha);

//' Generate a Dirichlet-Dirichlet-Gamma population (unsafe).
//'
//' Generate samples from m sources and p parameters, n sample per source.
//' The between-source alpha hyperparameter used to generate the source parameters is mandatory.
//'
//' ## RNG
//'
//' this function uses GSL's RNG seed, unaffected by R's RNG.
//'
//' @param n number of samples per source
//' @param m number of sources
//' @param alpha_0 between-source Gamma hyperparameter, a scalar
//' @param beta_0 between-source Gamma hyperparameter, a scalar
//' @param nu_0 between-source Dirichlet hyperparameter, a numeric vector
//' @export
//' @return a matrix with n*m rows
//' @inheritParams rdirichlet_cpp
//' @family RNG functions
// [[Rcpp::export]]
RcppGSL::Matrix rdirdirgamma_cpp(
      const unsigned int &n, const unsigned int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector &nu_0,
      const unsigned int seed = 0
);



//' Generate a Dirichlet-Dirichlet-Gamma population (safer).
//'
//' Generate samples from m sources and p parameters, n sample per source.
//' The between-source alpha hyperparameter used to generate the source parameters is mandatory.
//'
//' ## RNG
//'
//' this function uses R's RNG seed.
//'
//' @export
//' @return a matrix with n*m rows
//' @inheritParams rdirdirgamma_cpp
//' @family RNG functions
// [[Rcpp::export]]
Rcpp::NumericMatrix rdirdirgamma_beta_cpp(
      const unsigned int &n, const unsigned int &m,
      const double &alpha_0, const double &beta_0,
      const Rcpp::NumericVector nu_0
);


#endif
