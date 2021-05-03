#' Calculate partial identification bounds for the polyserial correlation
#'
#' https://www.rdocumentation.org/packages/polycor/versions/0.7-10/topics/polyserial
#'
#' @param copula Bivariate copula.
#' @param points The points where the second argument is evaluated.
#' @param marginals Specifies the marginal distribution. Either a named alternative
#'    or a list containing the distribution functions, quantile functions, and
#'    standard deviations for each marginal. Defaults to `"unifal"`, which
#'    uses standard unifal marginals.
#' @return The identification bounds.

polyserialiden = function(copula, points, marginals) {

  m_star = function(u, v) {

    candidates = function(u, v) {

      candidate_uv = function(point) {
        copula(cbind(u, rep(point, length(u)))) + pmax(v - copula(c(1, point)), 0)
      }

      min(sapply(points, candidate_uv))
    }

    pmin(u, v, mapply(FUN = candidates, u, v))
  }

  tau = c(0, points)

  quant1 = stats::qunif
  quant2 = stats::qunif
  dist1 = stats::punif
  dist2 = stats::punif
  dens1 = stats::dunif
  dens2 = stats::dunif
  sd1 = 1
  sd2 = 1

  upper_integrand = function(x) {
    u = x[1, ]
    v = x[2, ]
    matrix((m_star(u, v) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
  }

  upper_covariances = sapply(seq(length(points) - 1), function(i) cubature::hcubature(
      upper_integrand,
      lowerLimit = c(0, tau[i]),
      upperLimit = c(1, tau[i + 1]),
      vectorInterface = TRUE)$integral)

  sum(upper_covariances)

}


rho = -0
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
points = c(0.5, 1)
