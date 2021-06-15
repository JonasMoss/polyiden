pi = matrix(c(
  0.1, 0.1, 0.1,
  0.15, 0.15, 0.1,
  0.1, 0.1, 0.1
), nrow = 3, byrow = TRUE)



rho = 0.3
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))

pi = polyiden::generate_pi(
  i = I,
  j = J,
  cutoffs = "beta",
  copula = copula,
  shape1i = 2.7,
  shape2i = 1.1,
  shape1j = 3,
  shape2j = 5)

set.seed(313)
I = 3
J = 2

cum_pi = pi_to_cum_pi(pi)
u = sort(unique(c(cum_pi[I, ], cum_pi[, J])))

f = function(u, v) {
  pmin(upper_limit_cpp(u, v, cum_pi), upper_limit_cpp(v, u, cum_pi))
}

f = function(u, v) {
  pmax(lower_limit_cpp(u, v, cum_pi), lower_limit_cpp(v, u, cum_pi))
}

f = function(u, v) {
  lower_limit_cpp(u, v, cum_pi)
}

cum_pi_to_pi(outer(u, u, f))
pi
