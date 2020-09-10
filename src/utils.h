#ifndef __UTILS_H
#define __UTILS_H

#include <Rcpp.h>

using namespace Rcpp;


//' Compute column-wise sd
//'
//' @param mtx a nxp matrix
//' @return a 1xp vector
// [[Rcpp::export]]
Rcpp::NumericVector colsd(const Rcpp::NumericMatrix &mtx);


//' Compute the Minkowski norm of a vector
//'
//' Compute the Minkowski norm of a vector.
//' $p$ can range from 1 to infinity.
//'
//' @param v a vector
//' @param p exponent of the Minkowski norm (from 1 to Inf)
//' @return a double
// [[Rcpp::export]]
double norm_minkowski(const Rcpp::NumericVector &v, const double p = 2);



#endif
