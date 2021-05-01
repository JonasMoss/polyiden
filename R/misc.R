#' Make cumulative pi matrix.
#'
#' @param pi The matrix of probabilities
#' @return The cumulative distribution matrix.
cum_pi_matrix = function(pi) apply(apply(pi, 1, cumsum), 1, cumsum)

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
