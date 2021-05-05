#' Transform probability matrices.
#'
#' The reverse pi matrix has the property that its maximal (minimal) correlation
#'    equals the negative of the minimal (maximal) correlation of the original pi
#'    matrix.
#'
#' @param pi,cum_pi The matrix of probabilities or matrix or cumulative probability matrix
#'    of probabilities.
#' @return The cumulative probability matrix corresponding to `pi` if `pi_to_cumpi` is
#'    called; the probability matrix corresponding to `cum_pi` if `cum_pi_to_pi` is called. If
#'    `reverse_pi` is called, it reverses the `pi` matrix.
#' @examples
#'    library("polyiden")
#'    copula = function(u) copula::pCopula(u, copula::normalCopula(0.3, dim = 2))
#'    pi = polyiden::generate_pi(i = 3, j = 3, copula = copula)
#'    polyiden::pi_to_cum_pi(pi)
#'    all(polyiden::cum_pi_to_pi((pi_to_cum_pi(pi))) == pi) # TRUE
#'
#'    polyiden::polyiden(pi) # [1] -0.5699454  0.8549995
#'    -polyiden::polyiden(reverse_pi(pi)) # [1] 0.8549995 -0.5699454
#' @name transform_matrices
NULL

#' @rdname transform_matrices
#' @export
pi_to_cum_pi = function(pi) apply(apply(pi, 1, cumsum), 1, cumsum)

#' @rdname transform_matrices
#' @export
cum_pi_to_pi = function(cum_pi) {

  pi = cum_pi * 0
  I = nrow(cum_pi)
  J = ncol(cum_pi)
  pi[1, ] = c(cum_pi[1, 1], diff(cum_pi[1 , ]))
  pi[ , 1] = c(cum_pi[1, 1], diff(cum_pi[ , 1]))

  for (i in seq(I - 1) + 1) {

    for (j in seq(J - 1) + 1) {

      pi[i, j] = cum_pi[i, j] -
        cum_pi[i, j - 1] -
        cum_pi[i - 1, j] +
        cum_pi[i - 1, j - 1]

    }

  }

  pi

}

#' @rdname transform_matrices
#' @export
reverse_pi = function(pi) {

  pi_reverse = pi

  for (j in seq(ncol(pi_reverse))) {
    pi_reverse[ , j] = rev(pi[ , j])
  }

  pi_reverse

}
