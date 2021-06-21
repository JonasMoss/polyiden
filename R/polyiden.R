#' Calculate the partial identification set for the latent correlation
#'
#' Find the partial identification set for the latent caorrelation of a `pi`
#'    matrix when the marginals are known. If the marginals are unknown, the
#'    identification set is `(-1, 1)`.
#'
#' @param pi Matrix of probabilities.
#' @param marginals Specifies the marginal distribution. Either a named
#'    alternative or a list containing the distribution functions, quantile
#'    functions, and standard deviations for each marginal. (See the example.)
#'    Defaults to `"normal"`, which uses standard normal marginals.
#' @param method The method used for integration. One of `"direct"` or
#'    `"substitution"`. Defaults to `"substitution"`.
#' @param symmetric If `TRUE`, uses the approximation to the symmetric
#'    quasi-copula bounds to approximate the identification set when the
#'    copula is symmetric.
#' @param exact If `TRUE` and `symmetric = TRUE`, calculates the identification
#'    set for symmetric copulas. This may be extremely slow, and is not
#'    recommended. Usually yields results virtually equal to `symmetric = TRUE`
#'    and `exact = FALSE`.
#' @param max_rep Maximal repetitions when looking for an interior point; only
#'    used when `symmetric = TRUE` and `exact = TRUE`.
#' @param which_correlations Which correlations to return? Choose from `"lower"`,
#'    `"upper"`, and `"both"`.
#' @return The endpoints of the identification set.
#' @export
#' @examples
#'
#' # Make pi matrix from data.
#' data <- psychTools::bfi
#' data <- data[!is.na(data$A4), ]
#' data <- data[!is.na(data$A5), ]
#' pi <- table(data$A4, data$A5) / sum(table(data$A4, data$A5))
#'
#' # Find identification regions for the correlations for some marginals.
#' polyiden::polyiden(pi, marginals = "normal") # [1] -0.2624196  0.7602734
#' polyiden::polyiden(pi, marginals = "laplace") # [1] -0.3542015  0.7878503
#' polyiden::polyiden(pi, marginals = "uniform") # [1] -0.1668584  0.7203827
#' polyiden::polyiden(pi, marginals = "exponential") # [1] -0.2980032  0.8964526
#'
#' # Assuming symmetry, the identification sets become slightly smaller.
#' polyiden::polyiden(pi, marginals = "normal", symmetric = TRUE)
#' # [1] -0.1327938  0.6988629
#'
#' # Using custom t-distributed marginals with nu = 7 degrees of freedom.
#'
#' nu <- 7
#' tnu_marginals <- list(
#'   quant1 = function(p) qt(p, nu),
#'   quant2 = function(p) qt(p, nu),
#'   dens1 = function(q) dt(q, nu),
#'   dens2 = function(q) dt(q, nu),
#'   sd1 = sqrt(nu / (nu - 2)),
#'   sd2 = sqrt(nu / (nu - 2))
#' )
#'
#' polyiden::polyiden(pi, marginals = tnu_marginals) # [1] -0.3225405  0.7783182
#' @template references
polyiden <- function(
  pi,
  marginals = c("normal", "uniform", "exponential", "laplace"),
  method = c("substitution", "direct"),
  symmetric = FALSE,
  exact = FALSE,
  max_rep = 100000,
  which_correlations = c("both", "lower", "upper")) {

  method <- match.arg(method)
  which_correlatons <- match.arg(which_correlations)

  if (exact & symmetric) {
    return(c(
      symmetry_correlation(
        pi_to_cum_pi(pi),
        marginals = marginals,
        method = method,
        direction = "lower",
        max_rep = max_rep),
      symmetry_correlation(
        pi = pi_to_cum_pi(pi),
        marginals = marginals,
        method = method,
        direction = "upper",
        max_rep = max_rep)
    ))
  }

  if (is.character(marginals)) {
    marginals <- match.arg(marginals)

    if (marginals == "normal") {
      quant1 <- stats::qnorm
      quant2 <- stats::qnorm
      dist1 <- stats::pnorm
      dist2 <- stats::pnorm
      dens1 <- stats::dnorm
      dens2 <- stats::dnorm
      sd1 <- 1
      sd2 <- 1
    } else if (marginals == "uniform") {
      quant1 <- stats::qunif
      quant2 <- stats::qunif
      dist1 <- stats::punif
      dist2 <- stats::punif
      dens1 <- stats::dunif
      dens2 <- stats::dunif
      sd1 <- sqrt(1 / 12)
      sd2 <- sqrt(1 / 12)
    } else if (marginals == "exponential") {
      quant1 <- stats::qexp
      quant2 <- stats::qexp
      dist1 <- stats::pexp
      dist2 <- stats::pexp
      dens1 <- stats::dexp
      dens2 <- stats::dexp
      sd1 <- 1
      sd2 <- 1
    } else if (marginals == "laplace") {
      quant1 <- extraDistr::qlaplace
      quant2 <- extraDistr::qlaplace
      dist1 <- extraDistr::plaplace
      dist2 <- extraDistr::plaplace
      dens1 <- extraDistr::dlaplace
      dens2 <- extraDistr::dlaplace
      sd1 <- sqrt(2)
      sd2 <- sqrt(2)
    }
  } else {
    quant1 <- marginals$quant1
    quant2 <- marginals$quant2

    if (method == "direct") {
      dist1 <- marginals$dist1
      dist2 <- marginals$dist2
    }

    if (method == "substitution") {
      dens1 <- marginals$dens1
      dens2 <- marginals$dens2
    }

    sd1 <- marginals$sd1
    sd2 <- marginals$sd1
  }

  cum_pi <- pi_to_cum_pi(pi)
  k1 <- nrow(pi)
  k2 <- ncol(pi)

  if (symmetric) {
    upper_limit <- function(u, v, cum_pi) {
      pmin(upper_limit_cpp(u, v, cum_pi), upper_limit_cpp(v, u, cum_pi))
    }

    lower_limit <- function(u, v, cum_pi) {
      pmax(lower_limit_cpp(u, v, cum_pi), lower_limit_cpp(v, u, cum_pi))
    }
  } else {
    upper_limit <- upper_limit_cpp
    lower_limit <- lower_limit_cpp
  }

  if (method == "direct") {
    tau_i <- quant1(c(0, cum_pi[, k2]))
    tau_j <- quant2(c(0, cum_pi[k1, ]))

    upper_integrand <- function(x) {
      u <- dist1(x[1, ])
      v <- dist2(x[2, ])
      matrix(upper_limit(u, v, cum_pi) - u * v, nrow = 1)
    }

    lower_integrand <- function(x) {
      u <- dist1(x[1, ])
      v <- dist2(x[2, ])
      matrix(lower_limit(u, v, cum_pi) - u * v, nrow = 1)
    }
  }

  if (method == "substitution") {
    tau_i <- c(0, cum_pi[, k2])
    tau_j <- c(0, cum_pi[k1, ])

    upper_integrand <- function(x) {
      u <- x[1, ]
      v <- x[2, ]
      matrix((upper_limit(u, v, cum_pi) - u * v) /
               (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }

    lower_integrand <- function(x) {
      u <- x[1, ]
      v <- x[2, ]
      matrix((lower_limit(u, v, cum_pi) - u * v) /
               (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }
  }

  if(which_correlatons == "both" | which_correlatons == "upper") {
    upper_covariances <- outer(
      X = seq(nrow(pi)),
      Y = seq(ncol(pi)),
      FUN = Vectorize(function(i, j) {
        cubature::hcubature(
          upper_integrand,
          lowerLimit = c(tau_i[i], tau_j[j]),
          upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
          vectorInterface = TRUE
        )$integral
      })
    )
  }

  if(which_correlatons == "both" | which_correlatons == "lower") {
    lower_covariances <- outer(
      X = seq(nrow(pi)),
      Y = seq(ncol(pi)),
      FUN = Vectorize(function(i, j) {
        cubature::hcubature(
          lower_integrand,
          lowerLimit = c(tau_i[i], tau_j[j]),
          upperLimit = c(tau_i[i + 1], tau_j[j + 1]),
          vectorInterface = TRUE
        )$integral
      })
    )
  }

  if(which_correlatons == "both") {
    c(sum(lower_covariances), sum(upper_covariances)) / (sd1 * sd2)
  } else if(which_correlatons == "upper") {
    sum(upper_covariances) / (sd1 * sd2)
  } else if(which_correlatons == "lower") {
    sum(lower_covariances) / (sd1 * sd2)
  }

}

#' Calculate partial identification bounds for the polyserial correlation
#'
#' @param copula The copula.
#' @param points The points where the second argument is evaluated.
#' @param marginals Specifies the marginal distribution. Either a named
#'    alternative or a list containing the distribution functions, quantile
#'    functions, and standard deviations for each marginal. Defaults to
#'    `"normal"`, which uses standard normal marginals.
#' @param method The method used for integration. One of `direct"` or
#'    `"substitution`. Defaults to `"substitution`.
#' @export
#' @return The identification bounds.
#' @template references

polyserialiden <- function(
  copula,
  points,
  marginals = c("normal", "uniform", "exponential", "laplace"),
  method = c("substitution", "direct")) {

  method <- match.arg(method)

  if (is.character(marginals)) {
    marginals <- match.arg(marginals)

    if (marginals == "normal") {
      quant1 <- stats::qnorm
      quant2 <- stats::qnorm
      dist1 <- stats::pnorm
      dist2 <- stats::pnorm
      dens1 <- stats::dnorm
      dens2 <- stats::dnorm
      sd1 <- 1
      sd2 <- 1
    } else if (marginals == "uniform") {
      quant1 <- stats::qunif
      quant2 <- stats::qunif
      dist1 <- stats::punif
      dist2 <- stats::punif
      dens1 <- stats::dunif
      dens2 <- stats::dunif
      sd1 <- sqrt(1 / 12)
      sd2 <- sqrt(1 / 12)
    } else if (marginals == "exponential") {
      quant1 <- stats::qexp
      quant2 <- stats::qexp
      dist1 <- stats::pexp
      dist2 <- stats::pexp
      dens1 <- stats::dexp
      dens2 <- stats::dexp
      sd1 <- 1
      sd2 <- 1
    } else if (marginals == "laplace") {
      quant1 <- extraDistr::qlaplace
      quant2 <- extraDistr::qlaplace
      dist1 <- extraDistr::plaplace
      dist2 <- extraDistr::plaplace
      dens1 <- extraDistr::dlaplace
      dens2 <- extraDistr::dlaplace
      sd1 <- sqrt(2)
      sd2 <- sqrt(2)
    }
  } else {
    quant1 <- marginals$quant1
    quant2 <- marginals$quant2

    if (method == "direct") {
      dist1 <- marginals$dist1
      dist2 <- marginals$dist2
    }

    if (method == "substitution") {
      dens1 <- marginals$dens1
      dens2 <- marginals$dens2
    }

    sd1 <- marginals$sd1
    sd2 <- marginals$sd1
  }

  upper_limit <- function(u, v) {
    rep_points <- rep(points, each = length(u))
    rep_u <- rep(u, length(points))
    rep_v <- rep(v, length(points))
    rep_cop <- copula(cbind(rep_u, rep_points))

    objective <- matrix(
      data = rep_cop + pmax(rep_v - rep_points, 0),
      nrow = length(points),
      byrow = TRUE
    )

    pmin(u, v, Rfast::colMins(objective, value = TRUE))
  }

  lower_limit <- function(u, v) {
    rep_points <- rep(points, each = length(u))
    rep_u <- rep(u, length(points))
    rep_v <- rep(v, length(points))
    rep_cop <- copula(cbind(rep_u, rep_points))

    objective <- matrix(
      data = rep_cop - pmax(rep_points - rep_v, 0),
      nrow = length(points),
      byrow = TRUE
    )

    pmax(0, u + v - 1, Rfast::colMaxs(objective, value = TRUE))
  }

  if (method == "direct") {
    tau <- quant2(c(0, points, 1))
    lower <- quant1(0)
    upper <- quant1(1)

    upper_integrand <- function(x) {
      u <- dist1(x[1, ])
      v <- dist2(x[2, ])
      matrix(upper_limit(u, v) - u * v, nrow = 1)
    }

    lower_integrand <- function(x) {
      u <- dist1(x[1, ])
      v <- dist2(x[2, ])
      matrix(lower_limit(u, v) - u * v, nrow = 1)
    }
  } else if (method == "substitution") {
    tau <- c(0, points, 1)
    lower <- 0
    upper <- 1

    upper_integrand <- function(x) {
      u <- x[1, ]
      v <- x[2, ]
      matrix((upper_limit(u, v) - u * v) /
               (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }

    lower_integrand <- function(x) {
      u <- x[1, ]
      v <- x[2, ]
      matrix((lower_limit(u, v) - u * v) /
               (dens1(quant1(u)) * dens2(quant2(v))), nrow = 1)
    }
  }

  upper_covariances <- sapply(seq(length(tau) - 1), function(i) {
    cubature::hcubature(
      upper_integrand,
      lowerLimit = c(lower, tau[i]),
      upperLimit = c(upper, tau[i + 1]),
      vectorInterface = TRUE
    )$integral
  })

  lower_covariances <- sapply(seq(length(tau) - 1), function(i) {
    cubature::hcubature(
      lower_integrand,
      lowerLimit = c(lower, tau[i]),
      upperLimit = c(upper, tau[i + 1]),
      vectorInterface = TRUE
    )$integral
  })

  c(sum(lower_covariances), sum(upper_covariances)) / (sd1 * sd2)
}

#' Parametric latent correlations
#'
#' @param pi Matrix of probabilities.
#' @param distribution The joint distribution of the latent variables,
#'    parameterized by the correlation `"rho"`. Only `"normal"` supported.
#' @param discrepancy Discrepancy function to minimize. Defaults to
#'    `-y * log(x))`, which corresponds to minimization of the Kullback-Leibler
#'    divergence.
#' @return The closest latent correlation as judged by the discrepancy function.
#' @export
#' @template references

latent_correlation <- function(
  pi,
  distribution = c("normal"),
  discrepancy = function(x, y) -y * log(x)) {

  cum_pi <- pi_to_cum_pi(pi)
  k1 <- nrow(pi)
  k2 <- ncol(pi)

  tau_u <- cum_pi[, k2]
  tau_v <- cum_pi[k1, ]

  if (distribution == "normal") {
    distribution <- function(x, r) {
      copula::pCopula(x, copula::normalCopula(r, dim = 2))
    }
  }

  objective <- function(r) {
    probs <- outer(
      X = tau_u,
      Y = tau_v,
      FUN = Vectorize(function(u, v) distribution(c(u, v), r))
    )

    sum(discrepancy(cum_pi_to_pi(probs), pi))
  }

  stats::optimize(objective, interval = c(-1, 1))$minimum
}
