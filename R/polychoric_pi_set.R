#' Calculate the partial identification set for the latent correlation
#'
#' @param pi Matrix of probabilities.
#' @param dist_1,dist_2 Distributions.
#' @param quant_1,quant_1 Quantile functions of the distributions.
#' @param sd1,sd2 Standard deviations of the distributions.
#' @param symmetric `TRUE` if the copula is symmetric, false otherwise.

polychoric_pi_set = function(pi, dist1, dist2, quant1, quant2, sd1, sd2, symmetric = FALSE) {
  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)
  tau_i = quant1(c(0, cum_pi[, J]))
  tau_j = quant2(c(0, cum_pi[I, ]))

  upper_integrand = function(x) {
    u = dist1(x[1, ])
    v = dist2(x[2, ])
    matrix(upper_limit_cpp(u, v, cum_pi) - u * v, nrow = 1)
  }

  lower_integrand = function(x) {
    u = dist1(x[1, ])
    v = dist2(x[2, ])
    matrix(lower_limit_cpp(u, v, cum_pi) - u * v, nrow = 1)
  }

  upper_correlations = outer(seq(nrow(pi)), seq(ncol(pi)), Vectorize(function(i, j) {
    cubature::hcubature(
      upper_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
      vectorInterface = TRUE)$integral
  }))

  lower_correlations = outer(seq(nrow(pi)), seq(ncol(pi)), Vectorize(function(i, j) {
    cubature::hcubature(
      lower_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
      vectorInterface = TRUE)$integral
  }))

  c(sum(lower_correlations), sum(upper_correlations)) / (sd1 * sd2)

}

polychoric_pi_set_bad = function(pi, dist1, dist2, quant1, quant2, sd1, sd2, symmetric = FALSE) {
  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)
  tau_i = quant1(c(0, cum_pi[, J]))
  tau_j = quant1(c(0, cum_pi[I, ]))

  m_copula = upper_limit(pi)
  w_copula = lower_limit(pi)

  upper_integrand = function(x) {
    m_copula(dist1(x[1]), dist2(x[2])) - dist1(x[1]) * dist2(x[2])
  }

  lower_integrand = function(x, y) {
    w_copula(dist1(x[1]), dist2(x[2])) - dist1(x[1]) * dist2(x[2])
  }

  upper_cors = outer(seq(nrow(pi)), seq(ncol(pi)), Vectorize(function(i, j) {
    cubature::hcubature(
      upper_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]))$integral
  }))

  lower_cors = outer(seq(nrow(pi)), seq(ncol(pi)), Vectorize(function(i, j) {
    cubature::hcubature(
      lower_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]))$integral
  }))

  c(sum(lower_cors), sum(upper_cors)) / (sd1 * sd2)

}



