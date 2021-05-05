rm(list=ls())
require(lavaan)
library(mvtnorm)
require(copula)


K <- 7
rho <- 0.7

tau1 <- c(qnorm(seq(0.1, 0.9, length.out=K-1)), Inf)
tau2 <- c(qnorm(seq(0.05, 0.97, length.out=K-1)), Inf)
S <- expand.grid((1:K), (1:K))
S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs
#dim(S)[1] == K^2-1+2-2*K #check whether the number of elements are correct.

u1 <- pnorm(tau1[S[,1]])
u2 <- pnorm(tau2[S[,2]])
#Let us rather calculate the probability using the copula
#package, so that we can more easily change the copula:
#z <- pbivnorm(tau1[S[,1]],tau2[S[,2]], rho = rho)

#z <- pCopula(cbind(u1,u2), normalCopula(param=rho, dim = 2))

z <- pCopula(cbind(u1,u2), claytonCopula(param=3, dim = 2))

P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]], u1,u2,z))
colnames(P) <- c("i","j", "tau1", "tau2", "U1", "U2", "CDF")
P


#Let's find the gaussian copula which reproduce the CDF at each point.


for(i in (1:length(P$CDF))) {
  
  diffCDF = Vectorize(function(rhoCand) {
    res <- P$CDF[i]-pCopula(c(P$U1[i],P$U2[i]), normalCopula(param=rhoCand, dim = 2))
    return(res)
  })
  
  P$rho[i] <- uniroot(diffCDF, interval=c(-1,1))$root
}

P

#Visualize the results
copObjmvdc <- mvdc(claytonCopula(3), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
copObjmvdcGaus <- mvdc(normalCopula(0.7), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
contour(copObjmvdc, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), drawlabels=F)
contour(copObjmvdcGaus, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), add=T, col="red", drawlabels=F, lty=3)
for(i in (1:(K-1))) {
  abline(v=tau1[i])
  abline(h=tau2[i])
}
for(k in (1:length(P$i) )) {
  i <- P$i[k]
  j <- P$j[k]
  text(tau1[i], tau2[j], round(P$rho[k], 2), cex=1)
}


#Let us now see how the correlations per box are.

#Code from Jonas:
