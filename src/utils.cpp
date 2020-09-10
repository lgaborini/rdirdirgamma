#include <Rcpp.h>

using namespace Rcpp;


//' Compute column-wise sd
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



