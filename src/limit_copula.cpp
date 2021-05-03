#include "shared.h"

// [[Rcpp::export]]
double lower_limit_cpp_point_wise(const double &u, const double &v,
                                  const arma::mat &cum_pi, const int &I, const int &J) {

  arma::uvec i_vec = arma::find(cum_pi.col(J - 1) >= u);
  arma::uvec j_vec = arma::find(cum_pi.row(I - 1) >= v);

  int i_start = i_vec(0);
  int j_start = j_vec(0);
  int i_end = i_vec(0);
  int j_end = j_vec(0);

  if (i_start == I - 1) {
    i_start = I - 2;
    i_end = I - 2;
  } else if (i_start > 0)  {
    i_start = i_start - 1;
  }

  if (j_start == J - 1) {
    j_start = J - 2;
    j_end = J - 2;
  } else if (j_start > 0) {
    j_start = j_start - 1;
  }

  double current = -1/0.0;

  for (int i = i_start; i <= i_end; i++) {
    for (int j = j_start; j <= j_end; j++) {
      double candidate = cum_pi(i, j) - std::max(cum_pi(i, J - 1) - u, 0.0) - std::max(cum_pi(I - 1, j) - v, 0.0);
      current = std::max(candidate, current);
    }
  }

  return current;

}

// [[Rcpp::export]]
double upper_limit_cpp_point_wise(const double &u, const double &v,
                                  const arma::mat &cum_pi, const int &I, const int &J) {

  arma::uvec i_vec = arma::find(cum_pi.col(J - 1) >= u);
  arma::uvec j_vec = arma::find(cum_pi.row(I - 1) >= v);

  int i_start = i_vec(0);
  int j_start = j_vec(0);
  int i_end = i_vec(0);
  int j_end = j_vec(0);

  if (i_start == I - 1) {
    i_start = I - 2;
    i_end = I - 2;
  } else if (i_start > 0)  {
    i_start = i_start - 1;
  }

  if (j_start == J - 1) {
    j_start = J - 2;
    j_end = J - 2;
  } else if (j_start > 0) {
    j_start = j_start - 1;
  }

  double current = 1/0.0;

  for (int i = i_start; i <= i_end; i++) {
    for (int j = j_start; j <= j_end; j++) {
      double candidate = cum_pi(i, j) + std::max(u - cum_pi(i, J - 1), 0.0) + std::max(v - cum_pi(I - 1, j), 0.0);
      current = std::min(candidate, current);
    }
  }

  return current;

}

// [[Rcpp::export]]
arma::vec lower_limit_cpp(const arma::vec &u, const arma::vec &v, const arma::mat &cum_pi) {

  arma::vec return_vector(u.n_elem);
  int I = cum_pi.n_rows;
  int J = cum_pi.n_cols;

  for(unsigned int index = 0; index < u.n_elem; index++) {
    double maximum = lower_limit_cpp_point_wise(u(index), v(index), cum_pi, I, J);
    return_vector(index) =  std::max(0.0, std::max(u(index) + v(index) - 1.0, maximum));
  }

  return return_vector;

}

// [[Rcpp::export]]
arma::vec upper_limit_cpp(const arma::vec &u, const arma::vec &v, const arma::mat &cum_pi) {

  arma::vec return_vector(u.n_elem);
  int I = cum_pi.n_rows;
  int J = cum_pi.n_cols;

  for(unsigned int index = 0; index < u.n_elem; index++) {
    double minimum = upper_limit_cpp_point_wise(u(index), v(index), cum_pi, I, J);
    return_vector(index) = std::min(u(index), std::min(v(index), minimum));
  }

  return return_vector;

}

