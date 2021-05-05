test_that("polyserialiden agrees with paper 1", {
  testthat::skip("skipping heavy tests")
  # Running the code in "4.R" with a = 5 and rho yields bounds equal to
  rho = 0.35484520107474348815657094746711663901805877685546875
  bounds = c(-0.330320439386838637929599826748017221689224243164062500,
             0.782925788731825944566367070365231484174728393554687500)

  copula = function(u) copula::pCopula(u, copula::normalCopula(rho, dim = 2))
  points = c(0.5) # this is "a"

  bounds_direct = polyserialiden(copula, points, method = "direct")
  bounds_substitution = polyserialiden(copula, points, method = "substitution")
  expect_equal(bounds, bounds_direct, tolerance = 10e4)
  expect_equal(bounds, bounds_substitution, tolerance = 10e4)
})
