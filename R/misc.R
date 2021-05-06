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
#' @template references
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

#' Generate pi matrices.
#'
#' @param i,j The number of rows and number of columns of the pi matrix.
#' @param cutoffs Method to generate the uniform cutoffs. Supports `"uniform"`
#'    and `"beta"`. If `"beta"`, pass the `shape1i`, `shape2i` (for the `x`
#'    axis) and `shape1j`, `shape2j` (for the `y` axis) arguments to
#'    `...`.
#' @param copula The copula to generate probabilities from.
#' @param ... Passed to `qbeta` if applicable.
#' @return The pi matrix.
#' @examples
#'    # Make a normal copula with correlation 0.3.
#'    norm_cop = copula::normalCopula(0.3, dim = 2)
#'    copula = function(u) copula::pCopula(u, norm_cop)
#'
#'    # Generate pi matrix for this copula.
#'    pi = polyiden::generate_pi(i = 10, j = 4, copula = copula)
#'    polycor::polychor(pi) # [1] 0.2999905
#'
#'    # Using the "beta" argument we can non-uniform cutoffs.
#'    pi = polyiden::generate_pi(i = 10, j = 9, cutoffs = "beta",
#'                            copula = copula, shape1i = 2, shape2i = 7,
#'                            shape1j = 1, shape2j = 1/2)
#'
#'    polycor::polychor(pi) # [1] 0.2999905
#' @template references
#' @export

generate_pi = function(i, j, cutoffs = c("uniform", "beta"), copula, ...) {

  cutoffs = match.arg(cutoffs)

  tau_u = seq(i) / i
  tau_v = seq(j) / j

  if (cutoffs == "beta") {

    dots = list(...)
    tau_u = stats::qbeta(tau_u, dots$shape1i, dots$shape2i)
    tau_v = stats::qbeta(tau_v, dots$shape1j, dots$shape2j)

  }

  cum_pi = copula(cbind(rep(tau_u, each = length(tau_v)), rep(tau_v, length(tau_u))))
  cum_pi_to_pi(matrix(cum_pi, nrow = i, byrow = TRUE))

}

