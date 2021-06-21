#' Symmetrize a cumulative pi matrix
#'
#' Symmetrizes a cumulative pi matrix. Assumes the marginals have no
#'    equal elements. Returns a symmetric matrix with `NA`s for elements that
#'    are not known.
#' @param pi A cumulative pi matrix.
#' @return The symmetrized cumulative pi matrix.
symmetrize <- function(pi) {
  k1 <- ncol(pi)
  k2 <- nrow(pi)
  indices <- sort(unique(c(pi[k2, ], pi[, k1])))
  rows <- length(indices)
  mat <- matrix(NA, ncol = rows, nrow = rows)
  mat[rows, ] <- indices
  mat[, rows] <- indices

  i_index <- rep(0, k1 - 1)
  j_index <- rep(0, k2 - 1)

  for (i in seq(k1 - 1)) {
    i_index[i] <- which(indices == c(pi[k2, i]))
  }

  for (j in seq(k2 - 1)) {
    j_index[j] <- which(indices == c(pi[j, k1]))
  }

  for (i in seq(i_index)) {
    for (j in seq(j_index)) {
      mat[i_index[i], j_index[j]] <- pi[j, i]
      mat[j_index[j], i_index[i]] <- pi[j, i]
    }
  }

  mat
}


#' Generate the constraints
#' @param pi The pi matrix.
#' @return The symmetry constraints.
symmetry_constraints <- function(pi) {
  pi_s <- symmetrize(pi)
  pi_s[upper.tri(pi_s, diag = FALSE)] <- 0
  rows <- nrow(pi_s)
  aug_pi_s <- rbind(rep(0, rows + 1), cbind(rep(0, rows), pi_s))

  constraints <- matrix(
    0, nrow = choose(rows + 1, 2),
    ncol = choose(rows + 2, 2))

  counter <- 1

  rowsp <- rows + 1

  f <- function(i, j) {
    i <- i - 1
    j <- j - 1
    (rowsp * (rowsp - 1) / 2) - (rowsp - j) * ((rowsp - j) - 1) / 2 + i + 1
  }

  for (i in 2:rowsp) {
    for (j in 2:i) {
      constraints[counter, f(i, j)] <- 1
      constraints[counter, f(i - 1, j - 1)] <- 1
      if (i > j) {
        constraints[counter, f(i - 1, j)] <- -1
        constraints[counter, f(i, j - 1)] <- -1
      } else {
        constraints[counter, f(i, j - 1)] <- -2
      }
      counter <- counter + 1
    }
  }

  aug_pi_s[upper.tri(aug_pi_s, diag = FALSE)] <- 0
  x_indices <- which(is.na(aug_pi_s), arr.ind = TRUE)
  x_indices <- apply(x_indices, 1, function(u) f(u[1], u[2]))
  ui <- constraints[, x_indices]

  values <- aug_pi_s[(!is.na(aug_pi_s) & !upper.tri(aug_pi_s, diag = FALSE))]
  ci <- constraints[, setdiff(seq(ncol(constraints)), x_indices)] %*% values

  # Remove 0 rows.
  indices <- !apply(ui, 1, function(x) sum(x^2) == 0)

  list(ci = -ci[indices], ui = ui[indices, ])
}


#' Reduce a symmetric matrix problem.
#'
#' The inequality constraints in the symmetric matrix problem sometimes has
#'    an empty interior. In this case `stats::constrOptim` will not work, and
#'    we have to remove redunant variables.
#'
#' @param constraints List of constraints.
#' @return List of reduced constraints.
reduce_constraints <- function(constraints) {
  const_rhs <- constraints$ci
  const_dir <- rep(">=", length(const_rhs))
  const_mat <- constraints$ui
  rows <- ncol(const_mat)
  upper <- matrix(NA, ncol = rows, nrow = rows)
  lower <- matrix(NA, ncol = rows, nrow = rows)
  variable <- rep(TRUE, rows)

  vec <- rep(0, rows)
  for (i in seq(rows)) {
    vec[i] <- -1
    upper[i, ] <- lpSolve::lp(
      direction = "min",
      objective.in = vec,
      const.mat = const_mat,
      const.dir = const_dir,
      const.rhs = const_rhs
    )$solution
    vec[i] <- 1
    lower[i, ] <- lpSolve::lp(
      direction = "min",
      objective.in = vec,
      const.mat = const_mat,
      const.dir = const_dir,
      const.rhs = const_rhs
    )$solution
    if (lower[i, i] == upper[i, i]) variable[i] <- FALSE
    vec[i] <- 0
  }

  candidate_data <- rbind(
    upper = upper[variable, variable],
    lower = lower[variable, variable]
  )

  weight <- stats::runif(nrow(candidate_data))

  list(
    candidate = colSums(candidate_data * weight / (sum(weight))),
    candidate_data = candidate_data,
    constants = colMeans(upper[!variable, !variable]),
    variable = variable
  )
}

#' Find the upper and lower limits to the correlation assuming symmetry.
#'
#' @param pi The cumulative pi matrix.
#' @param marginals Specifies the marginal distribution. Either a named
#'    alternative or a list containing the distribution functions, quantile
#'    functions, and standard deviations for each marginal. (See the example.)
#'    Defaults to `"normal"`, which uses standard normal marginals.
#' @param method The method used for integration. One of `"direct"` or
#'    `"substitution"`. Defaults to `"substitution"`.
#' @param direction Direction of the correaltion.
#' @param max_rep Maximal repetitions when looking for an interior point.
#' @return The upper / lower limit to the correlation.
symmetry_correlation <- function(
  pi,
  marginals = c("normal", "uniform", "exponential", "laplace"),
  method = c("substitution", "direct"),
  direction = c("upper", "lower"),
  max_rep = 100000) {

  direction <- match.arg(direction)

  constraints <- symmetry_constraints(pi)
  reduction <- reduce_constraints(constraints)

  # Make new symmetrized "pi_s" matrix with fewer NAs.
  pi_s <- symmetrize(pi)
  rows <- nrow(pi_s)
  aug_pi_s <- rbind(rep(0, rows + 1), cbind(rep(0, rows), pi_s))
  aug_pi_s[upper.tri(aug_pi_s, diag = FALSE)] <- 0
  x_indices <- which(is.na(aug_pi_s), arr.ind = TRUE)[!reduction$variable, ]
  aug_pi_s[x_indices] <- reduction$constants
  pi_s <- aug_pi_s[2:(rows + 1), 2:(rows + 1)]

  # Construct the new constraints.
  indices <- rowSums(constraints$ui[, reduction$variable]^2) != 0
  ui <- constraints$ui[indices, reduction$variable]
  ci <- constraints$ci[indices] -
    constraints$ui[indices, !reduction$variable] %*% reduction$constants

  hrep <- cbind(rep(0, length(ci)), -ci, -ui)
  hrep <- rcdd::redundant(cbind(rep(0, length(ci)), -ci, -ui),
    representation = "H"
  )$output

  ci <- hrep[, 2]
  ui <- hrep[, 3:ncol(hrep)]

  # Define the functions
  x_indices <- which(is.na(pi_s), arr.ind = TRUE)

  if (direction == "upper") {
    f <- function(x) {
      s <- pi_s
      s[x_indices] <- x
      s <- s + t(s) - diag(diag(s))
      -polyiden::polyiden(cum_pi_to_pi(s), marginals, method)[2]
    }
  } else {
    f <- function(x) {
      s <- pi_s
      s[x_indices] <- x
      s <- s + t(s) - diag(diag(s))
      polyiden::polyiden(cum_pi_to_pi(s), marginals, method)[1]
    }
  }

  grad <- function(x) numDeriv::grad(f, x)

  theta <- get_interior(ui, ci, reduction$candidate, max_rep = max_rep)
  value <- stats::constrOptim(
    theta = theta,
    f = f,
    ci = -ci,
    ui = -ui,
    grad = grad
  )

  if (direction == "upper") -value$value else value$value
}


#' Find interior point.
#'
#' @param ui Constraints matrix.
#' @param ci Bounds vector.
#' @param candidate Candidate solution vector.
#' @param max_rep Maximal times to randomly search for new solution.

get_interior <- function(ui, ci, candidate, max_rep) {
  msg <- paste0(
    "No interior starting point found. try to increase max_rep above ",
    max_rep)

  index <- which(ui %*% candidate >= ci)

  if (length(index) > 0) {
    sols <- matrix(0, nrow = length(index), ncol = ncol(ui))

    for (i in seq(length(index))) {
      new_sol <- NA

      for (j in seq(max_rep)) {
        sol <- lpSolve::lp(
          direction = "min",
          objective.in = stats::rnorm(ncol(ui)),
          const.mat = ui,
          const.dir = rep("<=", length(ci)),
          const.rhs = ci
        )$solution

        if ((ui %*% sol < ci)[index[i]]) {
          new_sol <- sol
          break
        }
      }

      if (is.na(new_sol[1])) stop(msg)
      sols[i, ] <- new_sol
    }

    ok <- "no"
    for (j in seq(max_rep)) {
      weights <- stats::runif(nrow(sols) + 1)
      sol <- colSums(weights * rbind(candidate, sols)) / sum(weights)

      if (all((ui %*% sol < ci))) {
        ok <- "yes"
        break
      }
    }

    if (ok == "no") stop(msg)
  } else {
    sol <- candidate
  }

  sol
}
