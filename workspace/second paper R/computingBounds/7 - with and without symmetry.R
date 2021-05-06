rm(list=ls())
set.seed(123)

require(pbivnorm)
require(cubature)

getBoundsNormalBase <- function(K, rho) {

  #Support functions
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
    return(hcubature(f=intFunct, lowerLimit=c(-radius, -radius), upperLimit=c(radius, radius), ...))
  }
  
  
    
  #K <- 11
  #rho <- 0.7
  
  start_time <- Sys.time()
  
  tau1 <- qnorm(seq(0.01, 1, length.out=K)) #We include 1, corresponding to the final threshold \tau_{i,K} which is always \infty for normal marginals, but we do not not include zero, which corresponds to \tau_{i,0} = -\infty, which is not needed here
  tau2 <- qnorm(seq(0.01, 1, length.out=K))
  
  S <- expand.grid((1:K), (1:K)) #Let us generate the indices where we have restrictions.
  S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs

  #1) symmetric case
  x1 <- pnorm(tau1[S[,1]] + runif(length(S[,1]), min=0, max=0.04))
  #x1 <- pnorm(tau1[S[,1]])
  x2 <- pnorm(tau2[S[,2]])

  z <- pbivnorm(tau1[S[,1]],tau2[S[,2]], rho = rho)
  P <- as.data.frame(cbind(S, x1,x2,z))
  colnames(P) <- c("j1", "j2", "PX1", "PX2", "C")

  upper <- CovIntNorm(M)$integral
  lower <- CovIntNorm(W)$integral
  
  plot(P$PX1, P$PX2, xlim=c(0,1), ylim=c(0,1))
  points(P$PX2, P$PX1, col="red")
  abline(a=0,b=1)
  
  #add the cases from symmetry to P, recompute limits
  z <- pbivnorm(tau2[S[,2]], tau1[S[,1]], rho = rho) #NB Reversed
  P <- as.data.frame(rbind(as.matrix(P), cbind(S[,2], S[,1], x2, x1,z))) #NB reversed
  
  upperSym <- CovIntNorm(M)$integral
  lowerSym <- CovIntNorm(W)$integral
  
  end_time <- Sys.time()
  res <- c(K, lower, upper, end_time - start_time)
  
  names(res) <- c("K","lower", "upper", "lowerSym", "upperSym", "time elapsed")
  return(res)
}

Kinx <- (1:10)*5
X <- NULL
for (K in Kinx) {
  res <- getBoundsNormalBase(rho=0.7, K = K)
  print(res)
  X <- rbind(X, res)
}
X <- as.data.frame(X)

plot(X$K, X$lower, type="b", ylim=c(0,1), xlab="K, number of categories", ylab="")
points(X$K, X$upper, type="b")
points(X$K, X$lowerSym, type="b", col="red")
points(X$K, X$upperSym, type="b", col="red")
abline(h=0.7, lty=3)

