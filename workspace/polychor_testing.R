# Indicator for reverse-coded items.
reverse = c(1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1)


bfi_to_pi = function(i, j) {
  data = psychTools::bfi
  data = data[!is.na(data[ , i]), ]
  data = data[!is.na(data[ , j]), ]

  x = data[ , i]
  y = data[ , j]

  if(reverse[i] == 1) {
    x = 7 - x
  }

  if(reverse[j] == 1) {
    y = 7 - y
  }

  table(x, y) / sum(table(x, y))
}


quant1 = qnorm
quant2 = qnorm
dist1 = pnorm
dist2 = pnorm
sd1 = 1
sd2 = 1

get_lim = Vectorize(function(i, j) {
  if(i == j) return(NA)
  pi = bfi_to_pi(i, j)
  lims = polychoric_pi_set(pi, dist1, dist2, quant1, quant2, sd1, sd2)
  lims[2] - lims[1]
})

get_upper= Vectorize(function(i, j) {
  if(i == j) return(NA)
  pi = bfi_to_pi(i, j)
  lims = polychoric_pi_set(pi, dist1, dist2, quant1, quant2, sd1, sd2)
  lims[2]
})

get_lower = Vectorize(function(i, j) {
  if(i == j) return(NA)
  pi = bfi_to_pi(i, j)
  lims = polychoric_pi_set(pi, dist1, dist2, quant1, quant2, sd1, sd2)
  lims[1]
})

m = 25
tictoc::tic()
lengths = outer(1:m, 1:m, get_lim)
lowers = outer(1:m, 1:m, get_lower)
uppers = outer(1:m, 1:m, get_upper)
tictoc::toc()

par(mar = c(5.1, 4.1, 4.1, 4.1)) # adapt margins
plot(lengths)
plot(lowers)
plot(uppers)

install.packages("")















cum_pi = cum_pi_matrix(pi)

c1 = lower_limit_2(pi)
c2 = lower_limit(pi)
c3 = lower_limit_3(pi)


d1 = upper_limit_2(pi)
d2 = upper_limit(pi)

u = runif(100)
v = runif(100)
c1(u, v)
c2(u, v)
#c3(u, v)
d1(u, v)
d2(u, v)
upper_limit_cpp(u, v, pi)

microbenchmark::microbenchmark(c1(u, v), c2(u, v), c3(u, v))

microbenchmark::microbenchmark(d1(u, v), d2(u, v))


maximand_j = function(j, i, v, cum_pi) {
  cum_pi[i, j] - max(cum_pi[nrow(cum_pi), j] - v, 0)
}



maximand_j = function(i, j, v, cum_pi) {
  cum_pi[i, j] - max(cum_pi[i, ncol(cum_pi)] - u, 0)
}


maximand_ij = Vectorize(function(i, j, u, v, cum_pi) {
  cum_pi[i, j] - max(cum_pi[i, ncol(cum_pi)] - u, 0) - max(cum_pi[nrow(cum_pi), j] - v, 0)
}, c("i", "j"))

u = 0.3
v = 0.4

find_max = Vectorize(function(u, v, cum_pi) {

  result = outer(seq(5), seq(5), maximand_ij, u = u, v = v, cum_pi = cum_pi)
  c(which(result == max(result), arr.ind = TRUE))[2]

}, c("u", "v"))

u = seq(0, 1, by = 0.1)
v = seq(0, 1, by = 0.1)

A = outer(u, v, find_max, cum_pi = cum_pi)



v = 0.5
i = 5
plot(cum_pi[nrow(cum_pi), ],
sapply(seq(6), maximand_j, i = i, v = v, cum_pi = cum_pi))
abline(v = v)

sapply(seq(6), maximand_j, i = i, v = v, cum_pi = cum_pi)


g = function(v) which.max(sapply(seq(6), maximand_j, i = i, v = v, cum_pi = cum_pi))

v = seq(0, 1, by = 0.01)
i = 5
#plot(v, sapply(v, g))
plot(v, cum_pi[nrow(cum_pi), sapply(v, g)])
abline(a = 0, b = 1)

v = 0.1
h = function(i) which.max(sapply(seq(6), maximand_j, i = i, v = v, cum_pi = cum_pi))
i = seq(6)
plot(i, sapply(i, h))
