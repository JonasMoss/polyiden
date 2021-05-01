data = psychTools::bfi
data = data[!is.na(data$A4), ]
data = data[!is.na(data$A5), ]
pi = table(data$A4, data$A5) / sum(table(data$A4, data$A5))

# pi = runif(4)
# pi = pi / sum(pi)
# pi = matrix(pi, nrow = 2)

polyiden(pi, marginals = "normal")
polyiden(pi, marginals = "uniform")
polyiden(pi, marginals = "exponential")
polyiden(pi, marginals = "laplace")

# Uisng custom t-distributed marginals with nu = 7 degrees of freedom.

nu = 7
tnu_marginals = list(
  quant1 = function(p) qt(p, nu),
  quant2 = function(p) qt(p, nu),
  dist1 = function(q) pt(q, nu),
  dist2 = function(q) pt(q, nu),
  sd1 = sqrt(nu/(nu - 2)),
  sd2 = sqrt(nu/(nu - 2))
)

polyiden(pi, marginals = tnu_marginals)

