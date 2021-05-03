rho = 0.7
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
pi = generate_pi(i = 5, j = 5, copula = copula)
polycor::polychor(pi)

reverse_pi(pi)
cum_pi_matrix(pi)
cum_pi_matrix(reverse_pi(pi))



polyiden::polyiden(pi)
polyiden::polyiden(reverser(pi))





pi = generate_pi(i = 100, j = 100, copula = copula)


polyiden::polyiden(pi)
polyiden::polyiden(reverser(pi))

