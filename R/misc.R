#' Make cumulative pi matrix.
#'
#' @param pi The matrix of probabilities
#' @return The cumulative distribution matrix.
cum_pi_matrix = function(pi) apply(apply(pi, 1, cumsum), 1, cumsum)
