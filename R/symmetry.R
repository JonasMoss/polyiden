#' Symmetrize a cumulative pi matrix
#'
#' Symmetrizes a cumulative pi matrix. Assumes the marginals have no
#'    equal elements. Returns a symmetric matrix with `NA`s for elements that
#'    are not known.
#' @param pi A cumulative pi matrix.
#' @return The symmetrized cumulative pi matrix.
symmetrize = function(pi) {
  I = ncol(pi)
  J = nrow(pi)
  indices = sort(unique(c(pi[J, ], pi[, I])))
  K = length(indices)
  mat = matrix(NA, ncol = K, nrow = K)
  mat[K, ] = indices
  mat[, K] = indices

  i_index = rep(0, I - 1)
  j_index = rep(0, J - 1)

  for(i in seq(I - 1)) {
    i_index[i] = which(indices == c(pi[J, i]))
  }

  for(j in seq(J - 1)) {
    j_index[j] = which(indices == c(pi[j, I]))
  }

  for(i in seq(i_index)) {
    for(j in seq(j_index)) {
      mat[i_index[i], j_index[j]] = pi[j, i]
      mat[j_index[j], i_index[i]] = pi[j, i]
    }
  }

  mat

}


#' Generate the constraints
symmetry_constraints = function(pi) {
  pi_s = symmetrize(pi)
  pi_s[upper.tri(pi_s, diag = FALSE)] = 0
  K = nrow(pi_s)
  aug_pi_s = rbind(rep(0, K + 1), cbind(rep(0, K), pi_s))

  constraints = matrix(0, nrow = choose(K + 1, 2), ncol = choose(K + 2, 2))
  counter = 1

  Kp = K + 1

  f = function(i, j) {
    i = i - 1
    j = j - 1
    (Kp * (Kp - 1) /2) - (Kp - j) * ((Kp - j) - 1) / 2 + i + 1
  }

  for(i in 2:Kp) {
    for(j in 2:i) {
      constraints[counter, f(i, j)] = 1
      constraints[counter, f(i - 1, j - 1)] = 1
      if(i > j) {
        constraints[counter, f(i - 1, j )] = -1
        constraints[counter, f(i , j - 1)] = -1
      } else {
        constraints[counter, f(i, j - 1)] = -2
      }
      counter = counter + 1
    }
  }

  aug_pi_s[upper.tri(aug_pi_s, diag = FALSE)] = 0
  x_indices = which(is.na(aug_pi_s), arr.ind = TRUE)
  x_indices = apply(x_indices, 1, function(u) f(u[1], u[2]))
  ui = constraints[, x_indices]

  values = aug_pi_s [(!is.na(aug_pi_s) & !upper.tri(aug_pi_s, diag = FALSE))]
  ci = constraints[, setdiff(seq(ncol(constraints)), x_indices)] %*% values

  # Remove 0 rows.
  indices = !apply(ui, 1, function(x) sum(x^2) == 0)

  list(ci = -ci[indices], ui = ui[indices, ])

}


#' Reduce a symmetric matrix problem.
#'
#' The inequality constraints in the symmetric matrix problem sometimes has
#'    an empty interior. In this case `constrOptim` will not work, and we have
#'    to remove redunant variables.
#'
#' @param constraints List of constraints.
#' @return Not sure yet.
reduce_constraints = function(constraints) {
  const.rhs = constraints$ci
  const.dir = rep(">=", length(const.rhs))
  const.mat = constraints$ui
  K = ncol(const.mat)
  upper = matrix(NA, ncol = K, nrow = K)
  lower = matrix(NA, ncol = K, nrow = K)
  variable = rep(TRUE, K)

  vec = rep(0, K)
  for(i in seq(K)) {
    vec[i] = -1
    upper[i, ] = lpSolve::lp(direction = "min",
                    objective.in = vec,
                    const.mat = const.mat,
                    const.dir = const.dir ,
                    const.rhs = const.rhs)$solution
    vec[i] = 1
    lower[i, ] = lpSolve::lp(direction = "min",
                    objective.in = vec,
                    const.mat = const.mat,
                    const.dir = const.dir ,
                    const.rhs = const.rhs)$solution
    if(lower[i, i] == upper[i, i]) variable[i] = FALSE
    vec[i] = 0
  }

  candidate_data = rbind(
    upper = upper[variable, variable],
    lower = lower[variable, variable])

  weight = runif(nrow(candidate_data))

  list(
    candidate = colSums(candidate_data * weight / (sum(weight))),
    candidate_data = candidate_data,
    constants = colMeans(upper[!variable, !variable]),
    variable = variable)
}

