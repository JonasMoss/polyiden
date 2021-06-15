pi = matrix(c(
  0.10, 0.30, 0.40,
  0.15, 0.50, 0.80,
  0.20, 0.60, 1.00
), nrow = 3, byrow = TRUE)


pi = matrix(c(
  0.05, 0.10, 0.17, 0.19, 0.20,
  0.10, 0.30, 0.32, 0.35, 0.40,
  0.17, 0.32, 0.40, 0.50, 0.60,
  0.19, 0.35, 0.50, 0.70, 0.80,
  0.20, 0.40, 0.60 ,0.80, 1.00
), nrow = 5, byrow = TRUE)


cum_pi = cum_pi_to_pi(pi)

f = function(u, v) {
  vals1 = upper_limit_cpp(u, v, pi)
  vals2 = upper_limit_cpp(1-u, 1-v, pi)
  pmin(vals1, vals2 + u + v - 1)
}

u = seq(0,1, by = 0.1)
cum_pi_to_pi(outer(u, u, f))


f(0.4, 0.2)
f(0.8, 0.6)
f(0.4, 0.6)
f(0.8, 0.2)



f(0.4, 0.2)
f(0.8, 0.6)
f(0.4, 0.6)
f(0.8, 0.2)



u = 0.23
v = 0.54
f(u, v)
f(1 - u, 1 - v) + u + v - 1
