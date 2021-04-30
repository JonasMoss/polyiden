data = psychTools::bfi
data = data[!is.na(data$A4), ]
data = data[!is.na(data$A5), ]

pi = table(data$A4, data$A5) / sum(table(data$A4, data$A5))

pi = runif(4)
pi = pi / sum(pi)
pi = matrix(pi, nrow = 2)

quant1 = qunif
quant2 = qunif
dist1 = punif
dist2 = punif
sd1 = sqrt(1/12)
sd2 = sqrt(1/12)

quant1 = qnorm
quant2 = qnorm
dist1 = pnorm
dist2 = pnorm
sd1 = 1
sd2 = 1

tictoc::tic()
polychoric_pi_set(pi, dist1, dist2, quant1, quant2, sd1, sd2, symmetric = FALSE)
tictoc::toc()

tictoc::tic()
polychoric_pi_set_3(pi, dist1, dist2, quant1, quant2, sd1, sd2, symmetric = FALSE)
tictoc::toc()

tictoc::tic()
polychoric_pi_set_4(pi, dist1, dist2, quant1, quant2, sd1, sd2, symmetric = FALSE)
tictoc::toc()
