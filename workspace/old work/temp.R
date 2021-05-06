rho = 0.35
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
points = c(0.5)
polyserialiden(copula, points, method = "direct")

polyserialiden(copula, points)
polyserialiden(copula, points, method = "direct")


set.seed(313)
u = runif(10)
v = runif(10)

microbenchmark::microbenchmark(m_star(u, v), m_star_(u, v), m_star_2(u, v))

set.seed(313)
u = runif(100)
v = runif(100)
microbenchmark::microbenchmark(m_star_(u, v), m_star_2(u, v))

rho = 0.35
copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
points = c(0.25, 0.5, 0.75, 1)

microbenchmark::microbenchmark({
  matrix(
    data = copula(cbind(rep_u, rep(points, each = length(u)))), nrow = length(points),
    byrow = TRUE)
}, {
  matrix(
    data = mvtnorm::pmvnorm(cbind(rep_u, rep(points, each = length(u)))), nrow = length(points),
    byrow = TRUE)
})


polyserialiden(copula, points, "")

microbenchmark::microbenchmark(
mvtnorm::pmvnorm()




# Testing cpp

m_star_2 = function(u, v) {

  rep_points = rep(points, each = length(u)) # Send in points.
  rep_u = rep(u, length(points)) # send in v

  rep_cop = matrix(
    data = copula(cbind(rep_u, rep(points, each = length(u)))), nrow = length(points),
    byrow = TRUE)

  upper_limit_polyserial_cpp(u, v, rep_cop, points)
}

rep_points = rep(points, each = length(u)) # Send in points.
rep_cop = matrix(copula(cbind(rep_u, rep_points)), nrow = length(points), byrow = TRUE) # Must send in this.
upper_limit_polyserial_cpp(u, v, rep_cop, points)
upper_limit_polyserial_cpp_point_wise(u[1], u[2], rep_cop[, 1], points)
