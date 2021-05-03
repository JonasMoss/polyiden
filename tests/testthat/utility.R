# Copied from the paper
# "PARTIAL  IDENTIFICATION  OF  LATENT  CORRELATIONS  WITHBINARY  DATA"

#Details:
#
# The function `lcorr` calculates exact bounds for the latent correlation
#    based on the 2 x 2 table `p` when both marginals are the same.
#
# Calling `lcorr` with `dist = pnorm` is equivalent to calling `tetrachoric`,
#    which calculates the bounds on the latent correlation when both marginals
#    are normal. Calling `lcorr` with `dist = punif` quivalent to
#    calling `spearman`, which calculates the bounds of Spearman's rho.
#
# Parameters:
#  p: A matrix of positive numbers summing to 1 or a two times two table.
#  dist: The marginal distribution on both axes.
#  quant: The marginal quantile function on both axes.
#  var: Variance of the marginal distribution.
#
# Examples:
# tetrachoric(smallpox) # c(-0.85  0.99)
# polycor::polychor(smallpox) # 0.59

lcorr <- function(p, dist = pnorm, quant = qnorm) {
  var = 1 # Don't change; does not work.
  p <- p / sum(p)
  dist_name <- substitute(dist)
  if (dist_name == quote(pnorm)) {
    # Use special implementation of normal here?
    tetrachoric(p)
  } else if (dist_name == quote(punif)) {
    lspearman(p)
  } else {
    c(
      lower = -lcorr_(p[, c(2, 1)], Fx = dist, Fy = dist, Qx = quant,
                      Qy = quant, var = var),
      upper = lcorr_(p, Fx = dist, Fy = dist, Qx = quant, Qy = quant, var = var)
    )
  }
}


#Function "lspearman". Computes the bounds from Proposition 2.
lspearman <- function(p) {
  p <- p / sum(p)
  c(
    lower = 6 * p[1, 1] * p[2, 2] * (p[1, 1] + p[2, 2]) - 1,
    upper = 1 - 6 * p[1, 2] * p[2, 1] * (p[1, 2] + p[2, 1])
  )
}

#Function "tetrachoric": Computes the bounds from Proposition 1,
#assuming standard normal marginals.
tetrachoric <- function(p) {
  c(
    lower = -lcorr_(p[, c(2, 1)]),
    upper = lcorr_(p)
  )
}

# Support-function "lcorr_", which numerically computes the integrals
# J_2, \ldots, J_8, as described in the online supplementary material.
#
# Parameters:
# p: A matrix of positive numbers summing to 1.
# Fx: The marginal distribution on the x axis.
# Fy: The marginal distribution on the y axis.
# Qx: The marginal quantile function on the x axis.
# Qy: The marginal quantile function on the y axis.
# var: Common variance of the marginals.

lcorr_ <- function(p, Fx = pnorm, Fy = pnorm, Qx = qnorm, Qy = qnorm, var = 1) {
  p <- p / sum(p)
  1 + var * (-J2(p, Fx, Fy, Qx, Qy) + J3(p, Fx, Fy, Qx, Qy) +
    J4(p, Fx, Fy, Qx, Qy) + J5(p, Fx, Fy, Qx, Qy) + J6(p, Fx, Fy, Qx, Qy) +
    J7(p, Fx, Fy, Qx, Qy) + J8(p, Fx, Fy, Qx, Qy))
}

J2 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1])
  a2 <- Qx(p[1, 1] + p[1, 2] + p[2, 1])
  integrand <- function(x) stats::integrate(Fy, a1, x)$value
  2 * stats::integrate(Vectorize(integrand), a1, a2)$value
}

J3 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1])
  a2 <- Qx(p[1, 1] + p[1, 2])
  b1 <- Qy(p[1, 1])
  b2 <- Qy(p[1, 1] + p[2, 1])
  C <- p[1, 1]
  C * (a2 - a1) * (b2 - b1)
}

J4 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1] + p[1, 2])
  a2 <- Qx(p[1, 1] + p[1, 2] + p[2, 1])
  b1 <- Qy(p[1, 1] + p[2, 1])
  b2 <- Qy(p[1, 1] + p[1, 2] + p[2, 1])

  stats::integrate(Fx, a1, a2)$value * (b2 - b1) +
    stats::integrate(Fy, b1, b2)$value * (a2 - a1) -
    (p[1, 1] + p[1, 2] + p[2, 1]) * (b2 - b1) * (a2 - a1)
}

J5 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1])
  a2 <- Qx(p[1, 1] + p[1, 2])
  b1 <- Qy(p[1, 1] + p[2, 1])
  b2 <- Qy(p[1, 1] + p[1, 2] + p[2, 1])
  D <- stats::integrate(function(x) Fx(x) * (a2 - x), a1, a2)$value
  D * (b2 - b1) / (a2 - a1)
}

J6 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1])
  a2 <- Qx(p[1, 1] + p[1, 2])
  b1 <- Qy(p[1, 1] + p[2, 1])
  b2 <- Qy(p[1, 1] + p[1, 2] + p[2, 1])
  integrand <- function(x) {
    stats::integrate(function(y) Fy(y) - p[2, 1],
      lower = b1,
      upper = b1 + (b2 - b1) * (x - a1) / (a2 - a1)
    )$value
  }
  stats::integrate(Vectorize(integrand), a1, a2)$value
}

J7 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1] + p[1, 2])
  a2 <- Qy(p[1, 1] + p[1, 2] + p[2, 1])
  b1 <- Qy(p[1, 1])
  b2 <- Qy(p[1, 1] + p[2, 1])
  integrand <- function(x) (b2 - b1) / (a2 - a1) * (x - a2) * (p[1, 2] - Fx(x))
  stats::integrate(integrand, a1, a2)$value
}

J8 <- function(p, Fx, Fy, Qx, Qy) {
  a1 <- Qx(p[1, 1] + p[1, 2])
  a2 <- Qy(p[1, 1] + p[1, 2] + p[2, 1])
  b1 <- Qy(p[1, 1])
  b2 <- Qy(p[1, 1] + p[2, 1])
  integrand <- function(x) {
    stats::integrate(function(y) Fy(y),
      lower = b1,
      upper = b1 + (b2 - b1) * (x - a1) / (a2 - a1)
    )$value
  }
  stats::integrate(Vectorize(integrand), a1, a2)$value
}
