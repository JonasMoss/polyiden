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
#'      dist1 = function(q) pt(q, nu),
#'      dist2 = function(q) pt(q, nu),
#'      sd1 = sqrt(nu/(nu - 2)),
#'      sd2 = sqrt(nu/(nu - 2))
#'    )
#'
#'    polyiden::polyiden(pi, marginals = tnu_marginals) # [1] -0.3225405  0.7783182
#'
polyiden = function(pi, marginals = c("normal", "uniform", "exponential", "laplace"), symmetric = FALSE) {

  if (is.character(marginals)) {

    marginals = match.arg(marginals)

    if (marginals == "normal") {

      quant1 = stats::qnorm
      quant2 = stats::qnorm
      dist1 = stats::pnorm
      dist2 = stats::pnorm
      sd1 = 1
      sd2 = 1

    } else if (marginals == "uniform") {

      quant1 = stats::qunif
      quant2 = stats::qunif
      dist1 = stats::punif
      dist2 = stats::punif
      sd1 = sqrt(1/12)
      sd2 = sqrt(1/12)

    } else if (marginals == "exponential") {

      quant1 = stats::qexp
      quant2 = stats::qexp
      dist1 = stats::pexp
      dist2 = stats::pexp
      sd1 = 1
      sd2 = 1

    } else if (marginals == "laplace") {

      quant1 = extraDistr::qlaplace
      quant2 = extraDistr::qlaplace
      dist1 = extraDistr::plaplace
      dist2 = extraDistr::plaplace
      sd1 = sqrt(2)
      sd2 = sqrt(2)

    }

  } else {

    quant1 = marginals$quant1
    quant2 = marginals$quant2
    dist1 = marginals$dist1
    dist2 = marginals$dist2
    sd1 = marginals$sd1
    sd2 = marginals$sd1

  }

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
#' @param copula Copula
#' @param marginals Specifies the marginal distribution. Either a named alternative
#'    or a list containing the distribution functions, quantile functions, and
#'    standard deviations for each marginal. Defaults to `"normal"`, which
#'    uses standard normal marginals.
#' @param distribution Joint distribution. Ignored if `copula` is specified.
#' @return The latent correlation.
latent_correlation = function(pi, copula, marginals, distribution) stop("not implemented")
