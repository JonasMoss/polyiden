#Illustration of extremal distributions via simulation.
#The simulated extremal distributions are shown to
#be correct by reproducing an example from Nelson (2007), 
#see Theorem 3.2.3, p.71 in that book. 
#The code also includes a numerical verification of the 
#numerical simplifications underlying the lcorr-package.

#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("utility.R")


require(copBasic)

Cu <- function(u,v, p00, a, b) {
  res <- p00 + max(u-a, 0) + max(v - b, 0)
  return(min(u,v,res))
}
Cl <- function(u,v, p00, a, b) {
  res <- p00 - max(a - u,0 ) - max(b - v, 0)
  return(max(0,u+v-1,res))
}

#sanity check of Cu and Cl:
p00 <- 0.2
a <- 0.6 
b <- 0.3
#The a,b notation is from Nelson, kept for direct comparision.
# a is the first coordinate of the restriction, b the second.
# Therefore, a = p01 + p00, b = p10 + p00

Cu(a, b, p00 = p00, a = a, b = b)
#0.2
Cl(a, b, p00 = p00, a = a, b = b)
#0.2

#code based on M.R and FRECHETcop.R in the copBasic R packaged:
FrechetConstrained <- function(u,v, para, ...) {
  p00=para$p00; tau1=para$tau1; tau2=para$tau2; upperBound = para$upperBound
  
  if (is.null(p00) | is.null(tau1) | is.null(tau2)) {
    warning("Please specify p00, tau1 and tau2.")
    return(NULL)
  }
  if(length(u) > 1 & length(v) > 1 & length(u) != length(v)) {
    warning("length u = ",length(u), " and length v = ",length(v))
    warning("longer object length is not a multiple of shorter object length, no recycling in M()")
    return(NA)
  }
  if(length(u) == 1) {
    u <- rep(u, length(v))
  } else if(length(v) == 1) {
    v <- rep(v, length(u))
  }
  
  if (upperBound == TRUE) {
    return(sapply(1:length(u), function(i) { Cu(u[i], v[i], p00, pnorm(tau1), pnorm(tau2)) }))
  } else {
    return(sapply(1:length(u), function(i) { Cl(u[i], v[i], p00, pnorm(tau1), pnorm(tau2)) }))
  }
}

#Example from the Nelson book

p00 <- 0.2
tau1 <- qnorm(0.6)
tau2 <- qnorm(0.3)
a <- pnorm(tau1)
b <- pnorm(tau2)
print(tau1)
print(tau2)

#We require p00 to be compatible with a copula evaluated
#at the point a and b, which can be tested via the frechet
#bounds:
(max(a + b - 1,0) <= p00) & (p00 <= min(a,b))

set.seed(30092019)
plot(NA,NA, type="n", xlim=c(0,1), ylim=c(0,1), xlab="U", ylab="V")
for(i in 1:500) {
  v <- runif(1); t <- runif(1)
  points(derCOPinv2(cop=P,v,t),v, cex=0.5, pch=16) # black dots
}
for(i in 1:500) {
  v <- runif(1); t <- runif(1)
  u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, tau1=tau1, tau2=tau2, upperBound=TRUE))
  points(u,v, cex=2, pch=16, col=2) # red dots
}
for(i in 1:500) {
  v <- runif(1); t <- runif(1)
  u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, tau1=tau1, tau2=tau2, upperBound=FALSE))
  points(u,v, cex=2, pch=16, col=3) # red dots
}
#reproduces the plot on p 71 in Nelson

#Let us now compute the correlation bounds for this example
#in two ways, firstly, via a direct numerical integration (simulation)
#and then via the lcorr-package. 
N <- 500000
set.seed(30092019)

Xupper <- NULL
Xlower <- NULL
Uupper <- NULL
Ulower <- NULL
for(i in 1:N) {
  v <- runif(1); t <- runif(1)
  u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, tau1=tau1, tau2=tau2, upperBound=TRUE))
  
  X <- c(qnorm(u),qnorm(v))
  
  if(sum(is.finite(X)) == 2) {
    Uupper <- rbind(Uupper, c(u,v))
    Xupper <- rbind(Xupper, X)
  }
}
for(i in 1:N) {
  v <- runif(1); t <- runif(1)
  u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, tau1=tau1, tau2=tau2, upperBound=FALSE))
  
  X <- c(qnorm(u),qnorm(v))
  
  if(sum(is.finite(X)) == 2) {
    Ulower <- rbind(Ulower, c(u,v))
    Xlower <- rbind(Xlower, X)
  }
}

#To visualize the extremal distributions:
#par(mfrow=c(2,2))
#plot(Uupper, xlab="U", ylab="V", main="Upper bound, uniform marginals")
#plot(Ulower, xlab="U", ylab="V", main="Lower bound, uniform marginals")
#plot(Xupper, xlab=expression(xi[1]), ylab=expression(xi[2]), main="Upper bound, standard normal marginals")
#plot(Xlower, xlab=expression(xi[1]), ylab=expression(xi[2]), main="Lower bound, standard normal marginals")
#par(mfrow=c(1,1))

#Bounds with standard normal marginals:
boundsNormalMarginalsSim <- c(cor(Xlower)[1,2], cor(Xupper)[1,2])
#Bounds with uniform marginals:
boundsUniformMarginalsSim <- c(cor(Ulower)[1,2], cor(Uupper)[1,2])

#We now construct the implied 2x2 table 
#As mentioned above, the a,b-notation is from Nelson. 
# a is the first coordinate of the restriction, b the second.
# Therefore, a = p01 + p00, b = p10 + p00.
# Hence, p01 = a - p00, and p10 = b - p00. 
# Since the probabilities sum to one, 
# we have p11 = 1 - (p00 + p01 + p10).
p00 <- 0.2
a <- 0.6
p01 <- a - p00
p01 #0.4

b <- 0.3
p10 <- b - p00
p10 #0.1
p11 <- 1 - (p00 + p01 + p10)
p11

p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)
p

#The simulation based approximation is close to the 
#bound calculated via numerical integration (which is more precise)
tetrachoric(p)
#lower      upper 
#-0.8758248  0.9293887
boundsNormalMarginalsSim
#-0.8751501  0.9287861

lspearman(p)
#lower upper 
#-0.82  0.88 
boundsUniformMarginalsSim
#-0.8200603  0.8796400
