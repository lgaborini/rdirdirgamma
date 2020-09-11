#include "utils.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]



// Compute column-wise sd
Rcpp::NumericVector colsd(const Rcpp::NumericMatrix &mtx) {
   const unsigned int p = mtx.ncol();

   Rcpp::NumericVector sd_vec(p);

   for (unsigned int k = 0; k < p; ++k) {
      sd_vec[k] = Rcpp::sd(mtx(_, k));
   }
   return(sd_vec);
}



// Compute the Minkowski norm of a vector
double norm_minkowski(const Rcpp::NumericVector &v, const double p) {

   double r;

   if (p < 1) {
      Rcpp::stop("norm_minkowski: p cannot be < 1");
   }

   if (Rcpp::traits::is_infinite<REALSXP>(p)) {
      r = max(abs(v));
   } else {
      r = pow(sum(pow(abs(v), p)), 1/p);
   }
   return(r);
}

// Source: https://github.com/RfastOfficial/Rfast/tree/master/R
arma::rowvec colkurtosis(const arma::mat &mtx) {
   arma::rowvec vec_means = arma::mean(mtx, 0);
   arma::mat M = mtx;
   M.each_row() -= vec_means;
   M = arma::pow(M, 2);

   arma::rowvec vec_k = arma::mean(arma::pow(M, 2), 0) / arma::pow(arma::mean(M, 0), 2);
   return(vec_k);

}


// Source: https://github.com/RfastOfficial/Rfast/tree/master/R
arma::rowvec colskewness(const arma::mat &mtx) {
   const unsigned int n = mtx.n_rows;

   arma::rowvec vec_means = arma::mean(mtx, 0);
   arma::mat M = mtx;
   M.each_row() -= vec_means;

   arma::rowvec up = arma::sum(arma::pow(M, 3), 0);
   arma::rowvec down = arma::pow(
      arma::sum(arma::pow(M, 2), 0) / (n - 1),
      1.5
   );

   arma::rowvec vec_s = n * up / ((n - 1) * (n - 2) * down);
   return(vec_s);

}

