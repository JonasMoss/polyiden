#' The limit copulas satisfying the constraints of pi
#'
#' @param pi Matrix of probabilities.
#' @return The lower limit copula or upper limit copula.
#' @name limit_copula
NULL

#' @rdname limit_copula
lower_limit = function(pi) {

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  copula = Vectorize(function(u, v) {
    maximand = max(outer(seq(I - 1), seq(J - 1), Vectorize(function(i, j) cum_pi[i, j] - max(cum_pi[i, J] - u, 0) - max(cum_pi[I, j] - v, 0))))
    max(0, u + v - 1, maximand)
  })

  copula

}

#' @rdname limit_copula
upper_limit = function(pi) {

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  copula = Vectorize(function(u, v) {
    minimand = min(outer(seq(I - 1), seq(J - 1), Vectorize(function(i, j) cum_pi[i, j] + max(u - cum_pi[i, J], 0) + max(v - cum_pi[I, j], 0))))
    min(u, v, minimand)
  })

  copula

}

#' @rdname limit_copula
lower_limit_2 = function(pi) {

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  copula = Vectorize(function(u, v) {

    maximand = Vectorize(function(i, j) {
      cum_pi[i, j] - max(cum_pi[i, J] - u, 0) - max(cum_pi[I, j] - v, 0)
    })

    i = which(cum_pi[seq(I), J] >= u)[1]
    j = which(cum_pi[I, seq(J)] >= v)[1]

    if (i == 1) {
      i_indices = 1
    } else if (i == I) {
      i_indices = I - 1
    } else {
      i_indices = c(i, i - 1)
    }

    if (j == 1) {
      j_indices = 1
    } else if (j == J) {
      j_indices = J - 1
    } else {
      j_indices = c(j, j - 1)
    }

    maximum = max(outer(i_indices, j_indices, maximand))
    max(0, u + v - 1, maximum)

  })

  copula

}

#' @rdname limit_copula
upper_limit_2 = function(pi) {

  cum_pi = cum_pi_matrix(pi)
  I = nrow(pi)
  J = ncol(pi)

  copula = Vectorize(function(u, v) {

    minimand = Vectorize(function(i, j) {
      cum_pi[i, j] + max(u - cum_pi[i, J], 0) + max(v - cum_pi[I, j], 0)
    })

    i = which(cum_pi[seq(I), J] >= u)[1]
    j = which(cum_pi[I, seq(J)] >= v)[1]

    if (i == 1) {
      i_indices = 1
    } else if (i == I) {
      i_indices = I - 1
    } else {
      i_indices = c(i, i - 1)
    }

    if (j == 1) {
      j_indices = 1
    } else if (j == J) {
      j_indices = J - 1
    } else {
      j_indices = c(j, j - 1)
    }

    minimum = min(outer(i_indices, j_indices, minimand))
    min(u, v, minimum)

  })

  copula

}


#' @rdname limit_copula
lower_limit_3 = function(pi) function(u, v) c(lower_limit_cpp(u, v, pi))

#' @rdname limit_copula
upper_limit_3 = function(pi) function(u, v) c(upper_limit_cpp(u, v, pi))
