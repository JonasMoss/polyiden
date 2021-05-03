source("utility.R")

test_that("polyiden agrees with 2 x 2 Spearman", {

  rho = runif(1)
  copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
  pi = generate_pi(i = 2, j = 2, copula = copula)

  top = unname(polyiden(pi, marginals = "uniform", method = "direct"))
  bottom = unname(lspearman(pi))
  testthat::expect_equal(top, bottom, tolerance = 10e6)

})

test_that("polyiden agrees with 2 x 2 normal", {

  rho = runif(1)
  copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
  pi = generate_pi(i = 2, j = 2, copula = copula)

  top = unname(polyiden(pi, marginals = "normal", method = "direct"))
  bottom = unname(tetrachoric(pi))
  testthat::expect_equal(top, bottom, tolerance = 10e6)

})

test_that("latent_correlation works for normal", {

  rho = runif(1)
  copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
  pi = generate_pi(i = sample(2:20, 1), j = sample(2:20, 1), copula = copula)
  testthat::expect_equal(rho, latent_correlation(pi), tolerance = 10e6)

})




