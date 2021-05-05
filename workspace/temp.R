#' Calculate partial identification bounds for the polyserial correlation
#'
#' https://www.rdocumentation.org/packages/polycor/versions/0.7-10/topics/polyserial
#'
#' @param copula Bivariate copula.
#' @param points The points where the second argument is evaluated.
#' @param marginals Specifies the marginal distribution. Either a named alternative
#'    or a list containing the distribution functions, quantile functions, and
#'    standard deviations for each marginal. Defaults to `"normal"`, which
#'    uses standard unifal marginals.
#' @return The identification bounds.

polyserialiden = function(copula, points, marginals = c("normal")) {

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

  quant1 = stats::qnorm
  quant2 = stats::qnorm
  dist1 = stats::pnorm
  dist2 = stats::pnorm
  dens1 = stats::dnorm
  dens2 = stats::dnorm
  sd1 = 1
  sd2 = 1

  upper_integrand = function(x) {
    u = x[1, ]
    v = x[2, ]
    matrix((m_star(u, v) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
  }

  upper_covariances = sapply(seq(length(points)), function(i) cubature::hcubature(
      upper_integrand,
      lowerLimit = c(0, tau[i]),
      upperLimit = c(1, tau[i + 1]),
      vectorInterface = TRUE)$integral)

  sum(upper_covariances)

}


polyserialiden2 = function(copula, points, marginals = c("normal")) {

  m_star_ = function(u, v) {

    rep_points = rep(points, each = length(u)) # Send in points.
    rep_u = rep(u, length(points)) # send in v
    rep_v = rep(v, length(points)) # Send in u
    rep_cop = copula(cbind(rep_u, rep_points)) # Must send in this.

    objective = matrix(
      data = rep_cop + pmax(rep_v - rep_points, 0),
      nrow = length(points),
      byrow = TRUE)

    pmin(u, v, Rfast::colMins(objective, value = TRUE))

  }

  tau = c(0, points)

  quant1 = stats::qnorm
  quant2 = stats::qnorm
  dist1 = stats::pnorm
  dist2 = stats::pnorm
  dens1 = stats::dnorm
  dens2 = stats::dnorm
  sd1 = 1
  sd2 = 1

  upper_integrand = function(x) {
    u = x[1, ]
    v = x[2, ]
    matrix((m_star_(u, v) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
  }

  upper_covariances = sapply(seq(length(points)), function(i) cubature::hcubature(
    upper_integrand,
    lowerLimit = c(0, tau[i]),
    upperLimit = c(1, tau[i + 1]),
    vectorInterface = TRUE)$integral)

  sum(upper_covariances)

}

points = c(0.5, 1)
polyserialiden(copula, points)
polyserialiden2(copula, points)

set.seed(313)
u = runif(10)
v = runif(10)

microbenchmark::microbenchmark(m_star(u, v), m_star_(u, v), m_star_2(u, v))

set.seed(313)
u = runif(100)
v = runif(100)
microbenchmark::microbenchmark(m_star_(u, v), m_star_2(u, v))

rho = 0.35
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
points = c(0.25, 0.5, 0.75, 1)

polyserialiden(copula, points, "")



# Testing cpp

m_star_2 = function(u, v) {

  rep_points = rep(points, each = length(u)) # Send in points.
  rep_u = rep(u, length(points)) # send in v

  rep_cop = matrix(
    data = copula(cbind(rep_u, rep(points, each = length(u)))), nrow = length(points),
    byrow = TRUE)

  upper_limit_polyserial_cpp(u, v, rep_cop, points)
}

rep_points = rep(points, each = length(u)) # Send in points.
rep_cop = matrix(copula(cbind(rep_u, rep_points)), nrow = length(points), byrow = TRUE) # Must send in this.
upper_limit_polyserial_cpp(u, v, rep_cop, points)
upper_limit_polyserial_cpp_point_wise(u[1], u[2], rep_cop[, 1], points)
