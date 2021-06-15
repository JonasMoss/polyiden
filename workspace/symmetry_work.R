
constrOptim(
  theta = (optimum1 + optimum2)/2,
  f = f,
  ci = constraints$ci,
  ui = constraints$ui,
  grad = grad)


eps = 0.00000001
g = function(n_reps = 100000) {
  for(i in 1:n_reps) {
    samp = runif(length(optimum), -eps, eps)
    if(all(const.mat %*% (optimum + samp) > const.rhs))
      return(samp)

    return(FALSE)
  }
}

const.mat %*% (optimum + runif(length(optimum), -eps, eps)) >= const.rhs



cbind(const.mat, values = const.mat %*% optimum, const.rhs)



















pi = matrix(c(
  0.1, 0.1, 0.1,
  0.15, 0.1, 0.1,
  0.1, 0.15, 0.1
), nrow = 3, byrow = TRUE)

sum(pi)
pi = pi_to_cum_pi(pi)

pi = symmetrize(pi)

pi = matrix(c(
  0.1, 0.2,
  0.15, 0.1,
  0.2, 0.25
), nrow = 3, byrow = TRUE)

sum(pi)
pi = pi_to_cum_pi(pi)


pi = matrix(runif(4), nrow = 2)
pi = pi / sum(pi)
pi = pi_to_cum_pi(pi)
symmetry_constraints(pi)

pi = symmetrize(pi)
