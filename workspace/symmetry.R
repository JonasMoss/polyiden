pi = matrix(c(
  0.1, 0.1, 0.1,
  0.15, 0.1, 0.1,
  0.1, 0.15, 0.1
), nrow = 3, byrow = TRUE)

sum(pi)
cumpi = pi_to_cum_pi(pi)

polyiden::polyiden(pi)
polyiden::polyiden(pi, symmetric = TRUE)
polyiden::polyiden(reverse_pi(pi), symmetric = TRUE)

pi = matrix(c(
  0.4, 0.8,
  0.5, 1), nrow = 2, byrow = TRUE)

pi = cum_pi_to_pi(pi)

polyiden(pi, symmetric = TRUE)
polyiden(pi)

pi = matrix(c(
  0, 0.4, 0.5,
  0.4, 0.8, 0.8,
  0.5, 0.8,  1
), nrow = 3, byrow = TRUE)

cum_pi_to_pi(pi)

polyiden(cum_pi_to_pi(pi))




pi = matrix(c(
  0.5, 0.8,
  0.6, 1), nrow = 2, byrow = TRUE)

pi = cum_pi_to_pi(pi)

polyiden(pi, symmetric = TRUE)

allowed = Vectorize(function(x, y) {
  pi = matrix(c(
    x, 0.4, 0.5,
    0.4, y, 0.8,
    0.5, 0.8,  1
  ), nrow = 3, byrow = TRUE)
  if(all(cum_pi_to_pi(pi) >= 0)) {
    polyiden::polyiden(cum_pi_to_pi(pi))[2]
  } else {
    0
  }
})

x = seq(0.1, 0.4, by = 0.01)
y = seq(0.6, 0.7, by = 0.01)
z = outer(x, y, allowed)
contour(x, y, z)
abline(a = 0.8, b = -1, col = "red", lwd = 2)
min(z)
max(z)
