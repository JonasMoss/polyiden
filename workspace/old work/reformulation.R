polyiden2 = function(pi) {

  quant1 = stats::qnorm
  quant2 = stats::qnorm
  dist1 = stats::pnorm
  dist2 = stats::pnorm
  dens1 = stats::dnorm
  dens2 = stats::dnorm
  sd1 = 1
  sd2 = 1

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)
  tau_i = c(0, cum_pi[, J])
  tau_j = c(0, cum_pi[I, ])

  upper_integrand = function(x) {
    u = x[1, ]
    v = x[2, ]
    matrix((upper_limit_cpp(u, v, cum_pi) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
  }

  lower_integrand = function(x) {
    u = x[1, ]
    v = x[2, ]
    matrix((lower_limit_cpp(u, v, cum_pi) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
  }

  upper_covariances = outer(
    X = seq(nrow(pi)),
    Y = seq(ncol(pi)),
    FUN = Vectorize(function(i, j) cubature::hcubature(
      upper_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
      vectorInterface = TRUE)$integral))

  lower_covariances = outer(
    X = seq(nrow(pi)),
    Y = seq(ncol(pi)),
    FUN = Vectorize(function(i, j) cubature::hcubature(
      lower_integrand,
      lowerLimit = c(tau_i[i], tau_j[j]),
      upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
      vectorInterface = TRUE)$integral))

  c(sum(lower_covariances), sum(upper_covariances)) / (sd1 * sd2)

}
