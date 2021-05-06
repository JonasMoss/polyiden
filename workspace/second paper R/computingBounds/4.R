rm(list=ls())

require(pbivnorm)
require(copula)
require(cubature)

rho <- 0.7 #desired correlation

#Let's start by setting up a frame-work to generate distributions of X
#relative to an underlying Gaussian Z with given thresholds.

K <- 11
tau1 <- seq(-2.7,2.7, length.out=K)
tau2 <- tau1

#Let's rather make the thresholds uniform on the probability scale
tau1 <- qnorm(seq(0.01, 0.99, length.out=K))
tau2 <- tau1


S <- expand.grid((1:K), (1:K)) #Let us generate the indices where we have restrictions.
S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs
dim(S)[1] == K^2-1+2-2*K #check whether the number of elements are correct.

#Let us now add the left hand side of eq(3) to this matrix, then two
#new columns, which are the marginal probabilities used as inputs
#on the right of eq(3).

x1 <- pnorm(tau1[S[,1]])
x2 <- pnorm(tau2[S[,2]])
z <- pbivnorm(tau1[S[,1]],tau1[S[,2]], rho = rho)

P <- as.data.frame(cbind(S, x1,x2,z))
colnames(P) <- c("j1", "j2", "PX1", "PX2", "C")
P

plot(P$PX1, P$PX2)
abline(a=0,b=1) #points on the diagonal. Will be relevant for the symmetry illustration.

M <- function(u,v) {
  res <- min(u,v,min(P$C + (u-P$PX1)*(u-P$PX1 >= 0) + (v-P$PX2)*(v-P$PX2 >= 0)))
  return(res)
}
W <- function(u,v) {
  res <- max(0,u+v-1, max(P$C - (P$PX1-u)*(P$PX1-u >= 0) - (P$PX2-v)*(P$PX2-v >= 0)))
  return(res)
}
F1 <- function(x1) {
  return(pnorm(x1))
}
F2 <- function(x2) {
  return(pnorm(x2))
}


radius <- 4 #For the normal case, the double integral in the Hoeffding representation is over an infinite domain. This domain is approximated by [-radius,radius]^2
CovIntNorm <- function(C, ...) {
  intFunct <- function(x) {
    return(C(F1(x[1]), F2(x[2])) - F1(x[1])*F2(x[2]))
  }
  return(pcubature(f=intFunct, lowerLimit=c(-radius, -radius), upperLimit=c(radius, radius), ...))
}

print(date()) #[1] "Wed Feb 24 22:24:12 2021"
upper <- CovIntNorm(M)$integral
print(upper) #[1] 0.8425461
print(date()) #[1] "Wed Feb 24 22:25:31 2021"

print(date())  #[1] "Wed Feb 24 22:25:31 2021"
lower <- CovIntNorm(W)$integral
print(lower) #[1] 0.4584238
print(date()) #[1] "Wed Feb 24 22:30:16 2021"


c(lower, upper)
