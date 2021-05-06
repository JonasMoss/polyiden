rho = 0.35
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
js = c(2:9, 10 * seq(100))
pis = sapply(js, function(j) generate_pi(i = 4, j = j, copula = copula))
bounds = sapply(pis, polyiden::polyiden)

a = polyserialiden(copula, points)
b = -polyserialiden(function(u) copula::pCopula(u, copula::normalCopula(-rho, dim = 2)), points)
'
plot(js, bounds[1, ], ylim = c(-1, 1), log = "x", type = "l")
points(js, bounds[2, ], type = "l")
abline(h = a, lty = 2, col = "grey")
abline(h = b, lty = 2, col = "grey")
'
