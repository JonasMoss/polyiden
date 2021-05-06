#' Calculate the generalized polychoric correlation
#'
#' Computes the polychoric correlation using any discrepancy function.
#'
#' @param x A contingency table of counts or an ordered categorical variable; the
#'     latter can be numeric, logical, a factor, an ordered factor, or a
#'     character variable, but if a factor, its levels should be in proper
#'     order, and the values of a character variable are ordered alphabetically.
#' @param y If x is a variable, a second ordered categorical variable.
#' @param discrepancy The desired discrepancy function -ylog(x) correponds to the
#'     ordinary polychoric correlation.
#' @return The generalized polychoric correlation.
#'
polychoric = function(pi, discrepancy = function(x, y) (x - y)^2) {

  tau_x = qnorm(c(0, cumsum(colSums(pi))))
  tau_y = qnorm(c(0, cumsum(rowSums(pi))))

  ls = function(r) {

    probability_function = Vectorize(function(i, j) {
      c(mvtnorm::pmvnorm(
        lower = c(tau_y[i], tau_x[j]),
        upper = c(tau_y[i + 1], tau_x[j + 1]),
        corr = matrix(c(1, r, r, 1), nrow = 2)))
    })

    probabilities = outer(seq(ncol(pi)), seq(nrow(pi)), probability_function)

    sum(discrepancy(probabilities, pi))

  }

  optimize(ls, interval = c(-1, 1))$minimum

}



# Example
dat = psychTools::bfi[!is.na(psychTools::bfi$A5) & !is.na(psychTools::bfi$O5), ]
x = dat$A5
y = dat$O3
rm(dat)
pi = table(x, y) / sum(table(x, y))



polychoric(pi)
polychoric(pi, discrepancy = function(x, y) -y*log(x))
polychoric(pi, discrepancy = function(x, y) abs(x - y))
polychoric(pi, discrepancy = disc)
psych::polychoric(cbind(x, y))
polycor::polychor(x, y)

#polycor:::binBvn(r, tau_x[2:6], tau_y[2:6]) - probabilities

alphas = seq(0, 1, by = 0.01)
# BHHJ
disc = function(x, y, alpha) {
  if(alpha == 0) return(-y*log(x))
  x ^ (1 + alpha) - (1 + 1/alpha) * y * x ^alpha + 1/alpha * y ^ (1 + alpha)
}

result = sapply(alphas, function(alpha)
  polychoric(pi, discrepancy = function(x, y) disc(x, y, alpha)))

plot(alphas, result)
