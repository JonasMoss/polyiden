#ifndef SHARED
#define SHARED

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
arma::mat cum_pi_matrix_cpp(const arma::mat &pi);

#endif
