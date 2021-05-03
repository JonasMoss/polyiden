#' Lower bound copula for a rectangle
#'
#' @param lower,upper Lowe and upper limits of the rectangle. Both are vectors
#'    of two elements.
#' @param copula The copula cumulative distribution function to work with.
#' @return The lower bound copula function.

ubcopula = function(lower, upper, copula) {

  function(x) {
    fn = function(par) {
      copula(par) + max(x[1] - par[1], 0) + max(x[2] - par[2], 0)
    }

    par = c((lower[1] + upper[1])/2, (lower[2] + upper[2])/2)

    value = optim(
      par = par,
      fn = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper)$value

    min(x[1], x[2], value)
  }

}

lbcopula = function(lower, upper, copula) {

  function(x) {
    fn = function(par) {
      -(copula(par) - max(par[1] - x[1], 0) - max(par[2] - x[2], 0))
    }

    par = c((lower[1] + upper[1])/2, (lower[2] + upper[2])/2)

    value = -optim(
      par = par,
      fn = fn,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper)$value

    max(0, x[1] + x[2] - 1, value)
  }

}


#' Minimal and maximal correlations over a rectangle.
#' @param lower,upper Lowe and upper limits of the rectangle. Both are vectors
#'    of two elements.
#' @param copula The copula cumulative distribution function to work with.
#' @param M Integration range.
#' @param maxEval The maximum number of function evaluations needed. Defaults
#'    to `5000`.
#' @return The minimal correlation.

minimal_correlation = function(lower, upper, copula, M = 5, maxEval = 20000) {

  cop = lbcopula(lower, upper, copula)

  f = function(x) cop(c(pnorm(x[1]), pnorm(x[2]))) - pnorm(x[1])*pnorm(x[2])

  cubature::pcubature(
    f = f,
    lowerLimit = c(-M, -M),
    upperLimit = c(M, M),
    maxEval = maxEval)$integral

}

maximal_correlation = function(lower, upper, copula, M = 5, maxEval = 20000) {

  cop = ubcopula(lower, upper, copula)

  f = function(x) cop(c(pnorm(x[1]), pnorm(x[2]))) - pnorm(x[1])*pnorm(x[2])

  cubature::pcubature(
    f = f,
    lowerLimit = c(-M, -M),
    upperLimit = c(M, M),
    maxEval = maxEval)$integral

}

# Gaussian copula with correlation 0.5.
copula = function(x) copula::pCopula(x, copula::normalCopula(param = 0.5))

p = cbind(c(1, 0.5),
          c(0.5, 1))

coordinate = c(0.5, 0.5)
rho = 0.5

epses = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.49)

maximals = sapply(
  X = epses,
  FUN = function(eps) {
    maximal_correlation(coordinate - eps, coordinate + eps, copula)
  }
)

minimals = sapply(
  X = epses,
  FUN = function(eps) {
    minimal_correlation(coordinate - eps, coordinate + eps, copula)
  }
)

# tetrachorics is the function from utility.R

maximal_rho = tetrachoric(p)[2]
minimal_rho = tetrachoric(p)[1]

plot_bounds = function() {

  pdf("bounds.pdf", height = 8, width = 8)
  par(las = 1)
  plot(c(0, epses, 0.5), c(tetrachoric(p)[2], maximals, rho), ylim = c(-1, 1),
       type = "l", pch = 20, xlab = expression(epsilon),
       ylab = "Interval bounds", cex.lab = sqrt(2))
  points(c(0, epses, 0.5), c(tetrachoric(p)[1], minimals, rho), col = "red",
         type = "l", pch = 20)
  legend("bottomright", c("Upper limit", "Lower limit"), lwd = c(1, 1),
         pch = c(20, 20), col = c("black", "red"), bty = "n")
  dev.off()

}

plot_bounds()

saveRDS(maximals, file = "maximals.Rds")
saveRDS(minimals, file = "minimals.Rds")
