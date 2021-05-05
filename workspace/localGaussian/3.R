


require(pbv)

x1 = 1
x2 = 2
y1 = 0
y2 = 1.1

probMinus = Vectorize(function(rho) {
  prob <- pbv::pbvnorm(x1, y1, rho) +
    pbv::pbvnorm(x2, y2, rho) -
    pbv::pbvnorm(x2, y1, rho) -
    pbv::pbvnorm(x1, y2, rho)
  return(-prob)
})

probMatch <- 0.03

rho = seq(-1, 1, by = 0.01)
plot(rho, -probMinus(rho), type = "l")

maxOpt <- optim(par=0, fn=probMinus, lower=-1, upper=1, method="L-BFGS-B") #find the maximum
abline(v = maxOpt$par)

