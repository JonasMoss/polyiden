norm_cop = copula::normalCopula(0.7, dim = 2)
copula = function(u) copula::pCopula(u, norm_cop)

tau_u = 1:6 / 6
tau_v = 1:6 / 6

cumpi = outer(tau_u, tau_v, Vectorize(function(u, v) copula(c(u, v))))

pi = cum_pi_to_pi(cumpi)

polycor::polychor(pi)


pi = generate_pi(i = 10, j = 4, copula = copula)
polycor::polychor(pi)


ijs = 1:10
pis = sapply(ijs, function(ij) generate_pi(i = ij, j = ij, copula = copula))
intervals = sapply(pis, polyiden)

plot(ijs, intervals[1, ], ylim = c(-1, 1), type = "b", pch = 20)
points(ijs, intervals[2, ], type = "b", pch = 20)
abline(h = 0.7)

