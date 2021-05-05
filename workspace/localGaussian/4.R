rm(list=ls())
require(lavaan)


# _________________
#|     |     |     |
#| A   |  B  |  C  |
#|_____|_____|_____|
#|     |     |     |
#| D   |  E  |  F  |
#|_____|_____|_____|
#|     |     |     |
#| G   |  H  |  I  |
#|_____|_____|_____|

#Practically speaking, the matches should be done
#so that we avoid near-zero probabilities


#We match the CDFs in the upper right corners of
#D,E,G,H. The lower left corners at the "edges",
#i.e., A,B,C,F,I.

#G = (1,1)
#H = (1,2)
#D = (2,1)
#E = (2,2)

K <- 3
rho <- 0.7

tau1 <- c(qnorm(seq(0.1, 0.9, length.out=K-1)), 1)
tau2 <- c(qnorm(seq(0.05, 0.97, length.out=K-1)), 1)
S <- expand.grid((1:K), (1:K))
S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs
#dim(S)[1] == K^2-1+2-2*K #check whether the number of elements are correct.

x1 <- pnorm(tau1[S[,1]])
x2 <- pnorm(tau2[S[,2]])
z <- pbivnorm(tau1[S[,1]],tau2[S[,2]], rho = rho)
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]], x1,x2,z))
colnames(P) <- c("i","j", "tau1", "tau2", "P1", "P2", "CDF")
P




#Let's find the gaussian copula which reproduce the CDF at each point.
#Since we here have a gaussian distribution, it should be the same (correct underlying)
#correlation each time.

require(copula)


diffCDF = Vectorize(function(rhoCand) {
  res <- P$CDF[i]-pCopula(c(P$P1[i],P$P2[i]), normalCopula(param=c(rhoCand), dim = 2))
  return(res)
})

for(i in (1:length(P$CDF))) {
  P$rho <- uniroot(diffCDF, interval=c(-1,1))$root
}

P

library(mvtnorm)
x.points <- seq(-3,3,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- c(0,0)
sigma <- matrix(c(1,rho,rho,1),nrow=2, byrow=T)
for (i in 1:100) {
  for (j in 1:100) {
    z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                      mean=mu,sigma=sigma)
  }
}
contour(x.points,y.points,z,  drawlabels=F)

for(i in (1:K)) {
  abline(h=tau1[i])
  abline(v=tau2[i])
}

