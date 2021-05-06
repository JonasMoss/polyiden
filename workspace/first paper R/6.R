#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("utility.R")

#Approximating the minimal length of the bounds. 
#Again, we generate 2x2 tables from the bijection
#between latent normal variables with standard marginals
#and arbitrary thresholds.

#For tables with extreme dependencies, numerical integration (for normal marginals)
#is difficult. We therefore use simulation as rough approximations.


require(copBasic)

Cu <- function(u,v, p00, a, b) {
  res <- p00 + max(u-a, 0) + max(v - b, 0)
  return(min(u,v,res))
}
Cl <- function(u,v, p00, a, b) {
  res <- p00 - max(a - u,0 ) - max(b - v, 0)
  return(max(0,u+v-1,res))
}

#code based on M.R and FRECHETcop.R in the copBasic R packaged:
FrechetConstrained <- function(u,v, para, ...) {
  p00=para$p00; a=para$a; b=para$b; upperBound = para$upperBound
  
  if (is.null(p00) | is.null(a) | is.null(b)) {
    warning("Please specify p00, a and b.")
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
    return(sapply(1:length(u), function(i) { Cu(u[i], v[i], p00, a, b) }))
  } else {
    return(sapply(1:length(u), function(i) { Cl(u[i], v[i], p00, a, b) }))
  }
}

tetrachoricSim <- function(p, N = 1000) {
  p00 <- p[1,1]
  p01 <- p[1,2]
  p10 <- p[2,1]
  a <- p01 + p00
  b <- p10 + p00
  
  Xupper <- NULL
  Xlower <- NULL
  Uupper <- NULL
  Ulower <- NULL
  for(i in 1:N) {
    v <- runif(1); t <- runif(1)
    u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, a=a, b=b, upperBound=TRUE))
    
    X <- c(qnorm(u),qnorm(v))
    
    if(sum(is.finite(X)) == 2) {
      Uupper <- rbind(Uupper, c(u,v))
      Xupper <- rbind(Xupper, X)
    }
  }
  for(i in 1:N) {
    v <- runif(1); t <- runif(1)
    u <- derCOPinv2(v,t, cop=FrechetConstrained, para=list(p00=p00, a=a, b=b, upperBound=FALSE))
    
    X <- c(qnorm(u),qnorm(v))
    
    if(sum(is.finite(X)) == 2) {
      Ulower <- rbind(Ulower, c(u,v))
      Xlower <- rbind(Xlower, X)
    }
  }
  return(c(cor(Xlower)[1,2], cor(Xupper)[1,2]))
}
res <- NULL

M <- 20
for(tau1 in seq(-2,2, length.out=M)) {
  for(tau2 in seq(-2,2, length.out=M)) {
    for(rho in seq(-1,1, length.out=M)) {
      cat("tau1=",tau1," tau2=", tau2, " rho=", rho, "\n")
      a <- pnorm(tau1)
      b <- pnorm(tau2)
      p00 <- pCopula(copula = normalCopula(rho), u=c(a,b))
      
      p01 <- a - p00
      p10 <- b - p00
      p11 <- 1 - (p00 + p01 + p10)
      p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)
      
      Nbounds <- tetrachoricSim(p)
      Ubounds <- lspearman(p)
      res <- rbind(res, c(rho, Nbounds, Nbounds[2] - Nbounds[1], Ubounds, Ubounds[2] - Ubounds[1]))
    }
  }
}

colnames(res) <- c("rho", "Nlower", "Nupper", "Nlength",  "Ulower", "Uupper", "Ulength")
res <- data.frame(res)
#save.image("6.Rdata")

idmin <- which.min(res$Nlength)
res[idmin, ]
#       rho     Nlower     Nupper   Nlength Ulower    Uupper  Ulength
#3801  -1 -0.9999994 -0.3251844 0.6748151     -1 -0.489458 0.510542

idmin <- which.min(res$Ulength)
res[idmin, ]
#       rho    Nlower    Nupper   Nlength   Ulower Uupper  Ulength
#3800   1 0.2912296 0.9999981 0.7087685 0.489458      1 0.510542

#tau1=0, tau2=0 gives a slightly smaller minimum

## Computation of bounds for the identified minimum:
tau1 <- 0
tau2 <- 0
rho <- 1
a <- pnorm(tau1)
b <- pnorm(tau2)
p00 <- pCopula(copula = normalCopula(rho), u=c(a,b))

p01 <- a - p00
p10 <- b - p00
p11 <- 1 - (p00 + p01 + p10)
p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)

Nbounds <- tetrachoricSim(p, 10000)
Nbounds
#[1] 0.3255547 0.9999874
Nbounds[2] - Nbounds[1]
#0.6744327
Ubounds <- lspearman(p)
Ubounds
#0.5   1.0 
Ubounds[2] - Ubounds[1]
#0.5