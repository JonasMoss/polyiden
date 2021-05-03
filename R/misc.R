#' Make cumulative pi matrix.
#'
#' @param pi The matrix of probabilities
#' @return The cumulative distribution matrix.
cum_pi_matrix = function(pi) apply(apply(pi, 1, cumsum), 1, cumsum)

#' Transform cumulative pi matrix to pi matrix.
#'
#' @param cumpi Cumulative pi matrix.
#' @return pi matrix.

cum_pi_to_pi = function(cumpi) {

  pi_ = cumpi * 0
  I = nrow(cumpi)
  J = ncol(cumpi)
  pi_[1, ] = c(cumpi[1, 1], diff(cumpi[1 , ]))
  pi_[ , 1] = c(cumpi[1, 1], diff(cumpi[ , 1]))

  for (i in seq(I - 1) + 1) {
    for (j in seq(J - 1) + 1) {
      pi_[i, j] = cumpi[i, j] -
        cumpi[i, j - 1] -
        cumpi[i - 1, j] +
        cumpi[i - 1, j - 1]
    }
  }

  pi_

}

#' Reverse a pi matrix
#'
#' The reverse pi matrix has the property that its maximal (minimal) correlation
#'    equals the negative of the minimal (maximal) correlation of the original pi
#'    matrix.
#' @param pi A pi matrix
#' @return The reversed matrix

reverse_pi = function(pi) {

  pi_ = pi

  for (j in seq(ncol(pi_))) {
    pi_[ , j] = rev(pi[ , j])
  }

  pi_

}


#' The limit copulas satisfying the constraints of pi
#'
#' @param pi Matrix of probabilities.
#' @return The lower limit copula or upper limit copula.
#' @name limit_copula
NULL

#' @rdname limit_copula
lower_limit = function(pi) function(u, v) c(lower_limit_cpp(u, v, pi))

#' @rdname limit_copula
upper_limit = function(pi) function(u, v) c(upper_limit_cpp(u, v, pi))
