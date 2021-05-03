#' Calculate the partial identification set for the latent correlation
#'
#' Find the partial identification set for the latent caorrelation of a `pi`
#'    matrix when the marginals are known. If the marginals are unknown, the
#'    identification set is `(-1, 1)`.
#'
#' @param pi Matrix of probabilities.
#' @param marginals Specifies the marginal distribution. Either a named alternative
#'    or a list containing the distribution functions, quantile functions, and
#'    standard deviations for each marginal. Defaults to `"normal"`, which
#'    uses standard normal marginals.
#' @param method The method used for integration. One of `direct"` or `"substitution`. Defaults
#'    to `"substitution`.
#' @param symmetric `TRUE` if the copula is symmetric, false otherwise. Not yet
#'    implemented. Only `FALSE` is implemented.
#' @return The identification bounds.
#' @export
#' @examples
#'
#'    # Make pi matrix from data.
#'    data = psychTools::bfi
#'    data = data[!is.na(data$A4), ]
#'    data = data[!is.na(data$A5), ]
#'    pi = table(data$A4, data$A5) / sum(table(data$A4, data$A5))
#'
#'    # Find identification regions for the correlations for some marginals.
#'    polyiden::polyiden(pi, marginals = "normal") # [1] -0.2624196  0.7602734
#'    polyiden::polyiden(pi, marginals = "laplace") # [1] -0.3542015  0.7878503
#'    polyiden::polyiden(pi, marginals = "uniform") # [1] -0.1668584  0.7203827
#'    polyiden::polyiden(pi, marginals = "exponential") # [1] -0.2980032  0.8964526
#'
#'    # Uisng custom t-distributed marginals with nu = 7 degrees of freedom.
#'
#'    nu = 7
#'    tnu_marginals = list(
#'      quant1 = function(p) qt(p, nu),
#'      quant2 = function(p) qt(p, nu),
#'      dens1 = function(q) dt(q, nu),
#'      dens2 = function(q) dt(q, nu),
#'      sd1 = sqrt(nu/(nu - 2)),
#'      sd2 = sqrt(nu/(nu - 2))
#'    )
#'
#'    polyiden::polyiden(pi, marginals = tnu_marginals) # [1] -0.3225405  0.7783182
#'
polyiden = function(pi, marginals = c("normal", "uniform", "exponential", "laplace"),
                    method = c("substitution", "direct"), symmetric = FALSE) {

  method = match.arg(method)

  if (is.character(marginals)) {

    marginals = match.arg(marginals)

    if (marginals == "normal") {

      quant1 = stats::qnorm
      quant2 = stats::qnorm
      dist1 = stats::pnorm
      dist2 = stats::pnorm
      dens1 = stats::dnorm
      dens2 = stats::dnorm
      sd1 = 1
      sd2 = 1

    } else if (marginals == "uniform") {

      quant1 = stats::qunif
      quant2 = stats::qunif
      dist1 = stats::punif
      dist2 = stats::punif
      dens1 = stats::dunif
      dens2 = stats::dunif
      sd1 = sqrt(1/12)
      sd2 = sqrt(1/12)

    } else if (marginals == "exponential") {

      quant1 = stats::qexp
      quant2 = stats::qexp
      dist1 = stats::pexp
      dist2 = stats::pexp
      dens1 = stats::dexp
      dens2 = stats::dexp
      sd1 = 1
      sd2 = 1

    } else if (marginals == "laplace") {

      quant1 = extraDistr::qlaplace
      quant2 = extraDistr::qlaplace
      dist1 = extraDistr::plaplace
      dist2 = extraDistr::plaplace
      dens1 = extraDistr::dlaplace
      dens2 = extraDistr::dlaplace
      sd1 = sqrt(2)
      sd2 = sqrt(2)

    }

  } else {

    quant1 = marginals$quant1
    quant2 = marginals$quant2

    if(method == "direct") {
      dist1 = marginals$dist1
      dist2 = marginals$dist2
    }

    if(method == "substitution") {
      dens1 = marginals$dens1
      dens2 = marginals$dens2
    }

    sd1 = marginals$sd1
    sd2 = marginals$sd1

  }

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  if(symmetric) {
    upper_limit = upper_limit_symmetry_cpp
    lower_limit = lower_limit_symmetry_cpp
  } else {
    upper_limit = upper_limit_cpp
    lower_limit = lower_limit_cpp
  }

  if(method == "direct") {
    tau_i = quant1(c(0, cum_pi[, J]))
    tau_j = quant2(c(0, cum_pi[I, ]))

    upper_integrand = function(x) {
      u = dist1(x[1, ])
      v = dist2(x[2, ])
      matrix(upper_limit(u, v, cum_pi) - u * v, nrow = 1)
    }

    lower_integrand = function(x) {
      u = dist1(x[1, ])
      v = dist2(x[2, ])
      matrix(lower_limit(u, v, cum_pi) - u * v, nrow = 1)
    }

  }

  if(method == "substitution") {
    tau_i = c(0, cum_pi[, J])
    tau_j = c(0, cum_pi[I, ])

    upper_integrand = function(x) {
      u = x[1, ]
      v = x[2, ]
      matrix((upper_limit(u, v, cum_pi) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }

    lower_integrand = function(x) {
      u = x[1, ]
      v = x[2, ]
      matrix((lower_limit(u, v, cum_pi) - u * v) / (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }

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


#' Calculate partial identification bounds for the polyserial correlation
#'
#' https://www.rdocumentation.org/packages/polycor/versions/0.7-10/topics/polyserial
#'
#' @param pi Distribution function. Its first argument should be continous and its
#'    second argument categorical.
#' @param marginals Specifies the marginal distribution. Either a named alternative
#'    or a list containing the distribution functions, quantile functions, and
#'    standard deviations for each marginal. Defaults to `"normal"`, which
#'    uses standard normal marginals.
#' @return The identification bounds.
polyserialiden = function(pi, marginals) stop("not implemented")

#' Parametric latent correlations
#'
#' @param pi Matrix of probabilities.
#' @param distribution The joint distribution of the latent variables,
#'    parameterized by the correlation `"rho"`.
#' @param discrepancy Discrepancy function to minimize. Defaults to
#'    `-y * log(x))`, which corresponds to minimization of the Kullback-Leibler
#'    divergence.
#' @return The closest latent correlation as judged by the discrepancy function.
latent_correlation = function(pi, distribution = c("normal"), discrepancy = function(x, y) -y * log(x)) {

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  tau_u = cum_pi[, J]
  tau_v = cum_pi[I, ]

  if(distribution == "normal") {
    distribution = function(x, r) {
      copula::pCopula(x, copula::normalCopula(r, dim = 2))
    }
  }

  objective = function(r) {

    probs = outer(
      X = tau_u,
      Y = tau_v,
      FUN = Vectorize(function(u, v) distribution(c(u, v), r))
      )

    sum(discrepancy(cum_pi_to_pi(probs), pi))

  }

  stats::optimize(objective, interval = c(-1, 1))$minimum

}
