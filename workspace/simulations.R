rho = -0.2
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
pi = generate_pi(i = 7, j = 15, copula = copula)
polycor::polychor(pi)
polyiden::polyiden(pi, symmetric = TRUE)
polyiden2(pi, symmetric = TRUE)
polyiden::polyiden(pi, method = "substitution")

microbenchmark::microbenchmark(polyiden::polyiden(pi, method = "direct"),
                               polyiden::polyiden(pi, method = "substitution"))
polyiden::polyiden(pi, method = "direct")
polyiden::polyiden(pi, method = "substitution")


rho = 0.7
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))

ijs = 2:30
pis = sapply(ijs, function(ij) generate_pi(i = ij, j = ij, copula = copula))
intervals_normal = sapply(pis, polyiden, marginals = "normal")
intervals_uniform = sapply(pis, polyiden, marginals = "uniform")
intervals_laplace = sapply(pis, polyiden, marginals = "laplace")

plot(ijs, intervals_normal[1, ], ylim = c(-1, 1), type = "b", pch = 20)
points(ijs, intervals_normal[2, ], type = "b", pch = 20)
points(ijs, intervals_uniform [1, ], type = "b", pch = 20, col = "blue")
points(ijs, intervals_uniform [2, ], type = "b", pch = 20, col = "blue")
points(ijs, intervals_laplace[1, ], type = "b", pch = 20, col = "red")
points(ijs, intervals_laplace[2, ], type = "b", pch = 20, col = "red")
abline(h = rho)


rho = 0.7
norm_cop = copula::normalCopula(rho, dim = 2)
copula = function(u) copula::pCopula(u, norm_cop)

ijs = c(2:9, 10 * (1 : 10))
pis = sapply(ijs, function(ij) generate_pi(i = ij, j = ij, copula = copula))
tictoc::tic()
intervals_normal = sapply(pis, polyiden, marginals = "normal")
tictoc::toc()
intervals_uniform = sapply(pis, polyiden, marginals = "uniform")
intervals_laplace = sapply(pis, polyiden, marginals = "laplace")
#intervals_exponential = sapply(pis, polyiden, marginals = "exponential")

plot(ijs, intervals_normal[1, ], ylim = c(-1, 1), type = "b", pch = 20)
points(ijs, intervals_normal[2, ], type = "b", pch = 20)
points(ijs, intervals_uniform [1, ], type = "b", pch = 20, col = "blue")
points(ijs, intervals_uniform [2, ], type = "b", pch = 20, col = "blue")
points(ijs, intervals_laplace[1, ], type = "b", pch = 20, col = "red")
points(ijs, intervals_laplace[2, ], type = "b", pch = 20, col = "red")
#points(ijs, intervals_exponential[1, ], type = "b", pch = 20, col = "green")
#points(ijs, intervals_exponential[2, ], type = "b", pch = 20, col = "green")
abline(h = rho)

