#include "shared.hpp"

//' Make cumulative pi matrix.
//'
//' @param pi The matrix of probabilities
//' @return The cumulative distribution matrix.
// [[Rcpp::export]]

arma::mat cum_pi_matrix_cpp(const arma::mat &pi) {
  return arma::cumsum(arma::cumsum(pi, 0), 1);
}
