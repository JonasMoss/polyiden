#' Generate pi matrices.
#'
#' @param i,j The number of rows and number of columns of the pi matrix.
#' @param cutoffs Method to generate the uniform cutoffs. Only `"uniform"` right
#'    now.
#' @param copula The copula to generate probabilities from.
#' @return The pi matrix.
#' @examples
#'   norm_cop = copula::normalCopula(0.3, dim = 2)
#'   copula = function(u) copula::pCopula(u, norm_cop)
#'   pi = polyiden::generate_pi(i = 10, j = 4, copula = copula)
#'   polycor::polychor(pi) # [1] 0.2999905
#' @export


generate_pi = function(i, j, cutoffs = c("uniform"), copula) {

  cutoffs = match.arg(cutoffs)
  if(cutoffs == "uniform") {

    tau_u = 1:i / i
    tau_v = 1:j / j

  }

  cumpi = outer(tau_u, tau_v, Vectorize(function(u, v) copula(c(u, v))))
  cum_pi_to_pi(cumpi)

}
