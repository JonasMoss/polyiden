## Testing small.

pi = matrix(c(
  0.4, 0.8,
  0.5, 1), nrow = 2, byrow = TRUE)


inner2 = function(pi, direction = c("positive", "negative")) {

  direction = match.arg(direction)

  constraints = symmetry_constraints(pi)
  reduction = reduce_constraints(constraints)

  # Make new symmetrized "pi_s" matrix with fewer NAs.
  pi_s = symmetrize(pi)
  K = nrow(pi_s)
  aug_pi_s = rbind(rep(0, K + 1), cbind(rep(0, K), pi_s))
  aug_pi_s[upper.tri(aug_pi_s, diag = FALSE)] = 0
  x_indices = which(is.na(aug_pi_s), arr.ind = TRUE)[!reduction$variable, ]
  aug_pi_s[x_indices] = reduction$constants
  pi_s = aug_pi_s[2:(K+1), 2:(K+1)]

  # Construct the new constraints.
  indices = rowSums(constraints$ui[, reduction$variable]^2) != 0
  ui = constraints$ui[indices, reduction$variable]
  ci = constraints$ci[indices] -
    constraints$ui[indices, !reduction$variable] %*% reduction$constants

  # ui %*% reduction$candidate >= ci

  # Define the functions
  s = pi_s
  x_indices = which(is.na(pi_s), arr.ind = TRUE)

  if(direction == "positive") {
    f = function(x) {
      s[x_indices] = x
      s = s + t(s) - diag(diag(s))
      -polyiden::polyiden(cum_pi_to_pi(s))[2]
    }
  } else {
    f = function(x) {
      s[x_indices] = x
      s = s + t(s) - diag(diag(s))
      polyiden::polyiden(cum_pi_to_pi(s))[1]
    }
  }

  grad = function(x) numDeriv::grad(f, x)

  constrOptim(
    theta = reduction$candidate,
    f = f,
    ci = ci,
    ui = ui,
    grad = grad)

}

#' Find inner point
#'
#' Find an inner point satisfying the constraints implied by a pi matrix.
#'
#' @param pi Pi matrix.
#' @return Inner point
inner_point = function(pi) {
  pi_s = symmetrize(pi)
  pi_s[upper.tri(pi_s, diag = FALSE)] = 0
  x_indices = which(is.na(pi_s))
  K = sum(is.na(pi_s))

  reduce_constraints(constraints)

  constraints = symmetry_constraints(pi)
  const.rhs = constraints$ci
  const.dir = rep(">=", length(const.rhs))
  const.mat = constraints$ui

  objective.in = rep(0, K)
  optimas = matrix(NA, nrow = 2 * K, ncol = K)
  for(i in seq(2 * K)) {

    if(i > K) {
      objective.in[i - K] = -1
    } else {
      objective.in[i] = 1
    }

    optimas[i, ] = lpSolve::lp(
      direction = "min",
      objective.in = objective.in,
      const.mat = const.mat,
      const.dir = const.dir,
      const.rhs = const.rhs)$solution

    objective.in = rep(0, K)

  }

  const.mat %*% optimas[1, ] > const.rhs
  const.mat %*% optimas[2, ] > const.rhs
}

f = function(x) {
  s[x_indices] = x
  s = s + t(s) - diag(diag(s))
  polyiden::polyiden(cum_pi_to_pi(s))[1]
}

grad = function(x) {
  numDeriv::grad(f, x)
}

constraints = symmetry_constraints(pi)

# A list of RHS of constraints (except the non negative)

const.rhs = constraints$ci
const.dir = rep(">=", length(const.rhs))
const.mat = constraints$ui
objective.in <- rep(1, le)

# Find the optimal solution
optimum1 <-  lp(direction = "min",
                objective.in = c(0, 0, 0, 0, 1, 0),
                const.mat = const.mat,
                const.dir = const.dir ,
                const.rhs = const.rhs)$solution










pi = matrix(c(
  0.1, 0.13, 0.1,
  0.15, 0.12, 0.1,
  0.1, 0.1, 0.1
), nrow = 3, byrow = TRUE)
pi = pi_to_cum_pi(pi)

pi = matrix(c(
  runif(9)
), nrow = 3, byrow = TRUE)
pi = pi / sum(pi)
pi = pi_to_cum_pi(pi)

pi_s = symmetrize(pi)
pi_s[upper.tri(pi_s, diag = FALSE)] = 0
x_indices = which(is.na(pi_s))
s = pi_s

f = function(x) {
  s[x_indices] = x
  s = s + t(s) - diag(diag(s))
  polyiden::polyiden(cum_pi_to_pi(s))[1]
}

grad = function(x) {
  numDeriv::grad(f, x)
}

constraints = symmetry_constraints(pi)

# A list of RHS of constraints (except the non negative)
const.rhs = constraints$ci
const.dir = rep(">=", length(const.rhs))
const.mat = constraints$ui
objective.in <- rep(1, K)

# Find the optimal solution
optimum1 <-  lp(direction = "min",
                objective.in = c(0, 0, 0, 0, 1, 0),
                const.mat = const.mat,
                const.dir = const.dir ,
                const.rhs = const.rhs)$solution

# Find the optimal solution
optimum2 <-  lp(direction = "min",
                objective.in = runif(K),
                const.mat = const.mat,
                const.dir = const.dir ,
                const.rhs = const.rhs)$solution

constraints

hrep = cbind(rep(0, length(constraints$ci)), -constraints$ci, -constraints$ui)
rcdd::redundant(hrep, representation = "H")





pi = matrix(c(
  0.1, 0.13, 0.1,
  0.15, 0.12, 0.1,
  0.1, 0.1, 0.1
), nrow = 3, byrow = TRUE)

pi = matrix(c(
  0.1, 0.13, 0.1,
  0.15, 0.12, 0.14,
  0.1, 0.1, 0.06
), nrow = 3, byrow = TRUE)


dat = na.omit(cbind(psychTools::bfi$A3, psychTools::bfi$A4))
pi = table(dat[, 1], dat[, 2]) / sum(table(dat[, 1], dat[, 2]) )
pi = pi_to_cum_pi(pi)

constraints = symmetry_constraints(pi)
reduction = reduce_constraints(constraints)

pi_s = symmetrize(pi)
K = nrow(pi_s)
aug_pi_s = rbind(rep(0, K + 1), cbind(rep(0, K), pi_s))
aug_pi_s[upper.tri(aug_pi_s, diag = FALSE)] = 0
x_indices = which(is.na(aug_pi_s), arr.ind = TRUE)[!reduction$variable, ]
aug_pi_s[x_indices] = reduction$constants
Kp = K
pi_s = aug_pi_s[2:(K+1), 2:(K+1)]
constraints = symmetry_constraints(pi_s)

g = function(i, j) {
  i = i - 1
  j = j - 1
  (Kp * (Kp - 1) /2) - (Kp - j) * ((Kp - j) - 1) / 2 + i + 1
}

x_indices = which(is.na(pi_s), arr.ind = TRUE)
s = pi_s
pi_s = pi_s + t(pi_s * lower.tri(pi_s, diag = FALSE))
constraints = symmetry_constraints(pi_s)


f = function(x) {
  s[x_indices] = x
  s = s + t(s) - diag(diag(s))
  -polyiden::polyiden(cum_pi_to_pi(s))[2]
}

grad = function(x) {
  numDeriv::grad(f, x)
}


hrep = cbind(rep(0, length(constraints$ci)), -constraints$ci, constraints$ui)
vrep = rcdd::scdd(hrep)$output

weight = runif(15)
candidate = colSums(rcdd::scdd(hrep)$output[, 3:8] * weight / sum(weight))
constraints$ui %*% candidate >= constraints$ci

constrOptim(
  theta = candidate,
  f = f,
  ci = constraints$ci,
  ui = constraints$ui,
  grad = grad)


#
# hrep = cbind(rep(0, length(constraints$ci)), -constraints$ci, -constraints$ui)
# rcdd::redundant(hrep, representation = "H")$output
# reduction$candidate
#
#
# weight = runif(6)
# candidate = colMeans(rcdd::scdd(hrep)$output[, 3:8] * weight / sum(weight))
#
#
# constraints$ui %*% rcdd::scdd(hrep)$output[2, 3:8] <= -constraints$ci
