#Continuing the verification from 2.R, this time with
#a direct numerical integration routine. 


#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("utility.R")

M <- 5 #For the normal case, the double integral in the Hoeffding representation is over an infinite domain. This domain is approximated by [-M,M]^2
CovIntNorm <- function(C, ...) {
  intFunct <- function(x) {
    return(C(pnorm(x[1]), pnorm(x[2])) - pnorm(x[1])*pnorm(x[2]))
  }
  return(pcubature(f=intFunct, lowerLimit=c(-M, -M), upperLimit=c(M, M), ...))
}
varUnif <- (1/12)
CovIntUnif <- function(C, ...) {
  intFunctUnif <- function(x) {
    covRes <- C(x[1],x[2]) - x[1]*x[2]
    return(covRes/(varUnif)) #We divide the covariance with the square of the standard deviations to get correlation
  }
  return(adaptIntegrate(f=intFunctUnif, lowerLimit=c(0,0), upperLimit=c(1, 1), ...))
}

#sanity check of the integration routine
# (other numerical integration routines fail due to the jump
# discontinuities) : The standard Frechet bounds with
#standard normal marginals. Bounds ought to be [-1,1]
Cu <- function(u,v) {
  return(min(u,v))
}
Cl <- function(u,v) {
  return(max(0,u+v-1))
}

CovIntNorm(Cu, maxEval=50000)$integral #approx 1
CovIntNorm(Cl, maxEval=50000)$integral #approx -1

#Continuation of the example from 2.R:

p00 <- 0.2
a <- 0.6
p01 <- a - p00
b <- 0.3
p10 <- b - p00
p11 <- 1 - (p00 + p01 + p10)
p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)

CuConstrained <- function(u,v) {
  res <- p00 + max(u- a,0) + max(v - b,0) 
  return(min(u,v,res))
}
ClConstrained <- function(u,v) {
  res <- p00 - max(a - u,0 ) - max(b - v, 0)
  return(max(0,u+v-1,res))
}

numIntBoundNorm <- c(CovIntNorm(ClConstrained)$integral, CovIntNorm(CuConstrained)$integral)
numIntBoundNorm
#-0.8752675  0.9292124
tetrachoric(p)
#-0.8758248  0.9293887 

numIntBoundUnif <- c(CovIntUnif(ClConstrained)$integral, CovIntUnif(CuConstrained)$integral)
numIntBoundUnif
#[1] -0.8200010  0.8799984
lspearman(p)
#-0.82  0.88
