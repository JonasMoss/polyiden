#Comparing the biserial and dichotomous cases.

#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("workspace/first paper/utility.R")

#Suppose the distribution of Z is in fact exactly normal,
#with standard marginals.


#We use the same distribution as earlier as a base-line:
p00 <- 0.2
a <- 0.5
p01 <- a - p00
b <- 0.3
p10 <- b - p00
p11 <- 1 - (p00 + p01 + p10)
p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)

tau1 <- qnorm(a)
round(tau1,2)
tau2 <- qnorm(b)
round(tau2, 2)
rho <- polycor::polychor(p)
round(rho,2)

cop <- normalCopula(rho)

C <- function(u,v) {
  pCopula(copula = cop, u=c(u,v))
}
F1 <- function(x1) {
  return(pnorm(x1))
}
F2 <- function(x2) {
  return(pnorm(x2))
}
h <- function(v) { #This function has to be supplied by the user.
  return(C(F1(tau1), v))
}
Mp <- function(u,v) {
    min(u,v, h(v) + max(u - F1(tau1),0))
}
Wp <- function(u,v) {
  max(0,u+v-1, h(v) - max(F1(tau1) - u,0))
}

M <- 5 #For the normal case, the double integral in the Hoeffding representation is over an infinite domain. This domain is approximated by [-M,M]^2
CovIntNorm <- function(C, ...) {
  intFunct <- function(x) {
    return(C(F1(x[1]), F2(x[2])) - F1(x[1])*F2(x[2]))
  }

  return(hcubature(f=intFunct, lowerLimit=c(-M, -M), upperLimit=c(M, M), ...))
}
numIntBoundNorm <- c(CovIntNorm(Wp)$integral, CovIntNorm(Mp)$integral)
numIntBoundNorm
#-0.4939573  0.6805938
tetrachoric(p)
#lower      upper
#-0.8758248  0.9293887

# WITH A = 0.5: [1] -0.33  0.78

round(numIntBoundNorm, 2)
#[1] -0.49  0.68
