# we rather modify the (1,1)-box.
#For arbitrary positive correlation, we need one point at the (1,1) box, diverging to infinity. (done in 6.R)
#For arbitrary negative correlation, we need two points in the (1,1) box, say with \pi_{1,1} probability each,
#one at the upper left corner, and one at the bottom right corner.

#We here look at negative correlation.


rm(list=ls())
require(lavaan)
library(mvtnorm)
require(copula)
require(covsim) #Our VITA-implementation. Used here only in the bivariate case.
require(rvinecopulib)
require(psych)
require(rootSolve)

K.orig <- 7
K <- K.orig + 2
tau1 <- c(qnorm(c(
  seq(0.1, 0.9, length.out=K.orig-1), seq(0.901, 0.99, length.out=K-K.orig)
)), Inf)
tau2 <- tau1


rho <- -0.2


#Step 1; find the copula-parameter which matches up with the desired correlation.
sigma.target <- matrix(data=c(1,rho,rho,1), byrow=TRUE, ncol=2,nrow=2)
marginsnorm <- lapply(X=sqrt(diag(sigma.target)),function(X) list(distr="norm", sd=X) )
#family_set <- "clayton"
family_set <- "gauss"

set.seed(1234)
cv <- vita(marginsnorm, sigma.target =sigma.target, family_set = family_set)

#extracting the parameters is easier said than done
cvCop <- cv$copula$pair_copulas
theta <- as.numeric(cvCop[[1]][[1]]$parameters)


S <- expand.grid((1:K), (1:K)) 

#Our central object will be the matrix P, which we now start to construct.
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]]))
colnames(P) <- c("i","j", "tau1", "tau2") 
P

upperRight <- expand.grid(tau1, tau2)

#We now construct a matrix encoding the three remaining coordinates of each of the boxes
upperLeft <- upperRight*NA
lowerRight <- upperRight*NA
lowerLeft <- lowerRight*NA

#There is like a clever way to generate these automatically in the same way as
#upperRight, likely through expand.grid and reversing the elements one or both of tau1 tau2.
#I'll go for the simpler route here. I only actually use upperRight and lowerLeft
for(k in (1:dim(S)[1])) {
  indices <- S[k,]
  
  #upperLeft
  upperLeft[k,2] <- upperRight[k,2]
  if(indices[1] > 1) {
    upperLeft[k,1] <- tau1[S[k,1]-1]
  } else {
    upperLeft[k,1] <- -Inf
  }
  
  #lowerRight
  lowerRight[k,1] <- upperRight[k,1]
  if(indices[2] > 1) {
    lowerRight[k,2] <- tau2[S[k,2]-1]
  } else {
    lowerRight[k,2] <- -Inf
  }
  
  #lowerLeft
  if(indices[1] > 1) {
    lowerLeft[k,1] <- tau1[S[k,1]-1]
  } else {
    lowerLeft[k,1] <- -Inf
  }
  if(indices[2] > 1) {
    lowerLeft[k,2] <- tau2[S[k,2]-1]
  } else {
    lowerLeft[k,2] <- -Inf
  }
}


#Fill in box-probabilities in P
copObjmvdc <- mvdc(claytonCopula(theta), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
P$ProbBox <- NA
for(k in (1:dim(S)[1])) {
  
  x1 <- upperRight[k,1]
  y1 <- upperRight[k,2]
  x2 <- lowerLeft[k,1]
  y2 <- lowerLeft[k,2]
  
  P$ProbBox[k] <- pMvdc(c(x1,y1), copObjmvdc) +
    pMvdc(c(x2, y2), copObjmvdc) -
    pMvdc(c(x2, y1), copObjmvdc) -
    pMvdc(c(x1, y2), copObjmvdc)
}
P

getPoints <- function(P, M, eps = 0.001) { #M is the addition to the x and y coordinate of the upper right corner-point.
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  ypoints <- NULL
  xpoints <- NULL
  xlim <- NULL
  ylim <- NULL
  
  for(k in (1:length(P$i) )) {
    i <- P$i[k]
    j <- P$j[k]
    
    if(is.finite(tau1[i])) {
      if(i > 1) {
        xlim <- (tau1[i]+tau1[i-1])/2
      } else {
        xlim <- tau1[i] #we allow equality with tau_i in in \tau_{i,j-1} < \xi_i \leq \tau_{i,j}
      }
    } else {
      #in the tau1 and tau2 vectors, we include the +Inf cases, but not the -Inf cases.
      #Since tau1[i] is infinite, we are at the right-most edge.
      xlim <- tau1[i-1]+eps
    }
    
    #    We now do the same with the y-coordinate:
    if(is.finite(tau2[j])) {
      if(j > 1) {
        ylim <- (tau2[j]+tau2[j-1])/2
      } else {
        ylim <- tau2[j]
      }
    } else {
      ylim <- tau2[j-1]+eps
    }
    
    #special case for the lower left point:
    if ( (i == 1) & (j==1) ) {
      xlim <- xlim - M
      ylim <- ylim - M
    }
    
    xpoints <- c(xpoints, xlim)
    ypoints <- c(ypoints, ylim)
  }
  return(cbind(xpoints, ypoints))
}


getPointsNegCor <- function(P, M, eps = 0.001) {
  
  XY <- as.data.frame(getPoints(P,M=0, eps)) #will be modified later, and directly (I leave the previous code as is, as I might use it when trying to generalize this to a middle box)
  XYprob <- cbind(XY, P$ProbBox)
  colnames(XYprob) <- c("x", "y", "prob")
  
  #Let us now modify XY as planned in the introductory text
  
  #NB: hard-coded that the (1,1)-box is XY[1,]
  xorig <- XYprob$x[1]
  yorig <- XYprob$y[1]
  #assumes M > 0:
  if(M < 0) {
    stop("M < 0")
  }
  upperLeftMod <- c(-M, yorig, XYprob[1,3]/2)
  lowerRightMod <- c(xorig,-M, XYprob[1,3]/2)
  
  XYprobMod <- rbind(upperLeftMod, lowerRightMod, XYprob[-1,])
  return(XYprobMod)
}

#We modify the drawGrid-code
drawPoints <- function(P, M, eps = 0.001) {
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  
  XYprobMod <- getPointsNegCor(P, M, eps)
  XY <- data.frame(x=XYprobMod$x, y=XYprobMod$y)
  
  plot(NA, xlim=c(min(XY$x) - 5*eps, max(XY$x)+5*eps), ylim=c(min(XY$y) - 5*eps, max(XY$y)+5*eps), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
  
  
  for(i in (1:(K-1))) { #grid lines for the thresholds
    abline(v=tau1[i])
    abline(h=tau2[i])
  }
  points(XY)
}

drawPoints(P, M=0)
drawPoints(P, M=5)
drawPoints(P, M=10)
drawPoints(P, M=100)

#We then use this structure to directly compute the correlation implied by P.

#NB: this function is also ad-hoc modified. 

calcCor <- function(P,M, eps = 0.001) {
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  #XY <- as.data.frame(getPoints(P,M,eps))
  #colnames(XY) <- c("x","y")
  
  XYprobMod <- getPointsNegCor(P, M, eps)
  XY <- data.frame(x=XYprobMod$x, y=XYprobMod$y)
  
  Xexp <- 0
  Xexpsq <- 0
  Yexp <- 0
  Yexpsq <- 0
  
  #Calculate E X and E X^2, and use this to get the standard deviation
  for (xval in unique(XY$x)) {
    inx.x <- which(XY$x == xval)
    Xexp <- Xexp + xval*sum(XYprobMod$prob[inx.x])
    Xexpsq <- Xexpsq + (xval^2)*sum(XYprobMod$prob[inx.x])
  }
  stdX <- sqrt( Xexpsq - (Xexp)^2 )
  
  #Same with Y:
  for (yval in unique(XY$y)) {
    inx.y <- which(XY$y == yval)
    Yexp <- Yexp + yval*sum(XYprobMod$prob[inx.y])
    Yexpsq <- Yexpsq + (yval^2)*sum(XYprobMod$prob[inx.y])
  }
  stdY <- sqrt( Yexpsq - (Yexp)^2 )
  
  #The covariance, calculated vectorially
  covXY <- sum((XY$x - Xexp)*(XY$y - Yexp)*XYprobMod$prob)
  rhoXY <- covXY/(stdX*stdY)
  return(rhoXY)
}

calcCor(P, M=0)
calcCor(P, M=5)
calcCor(P, M=10)
calcCor(P, M=100)
calcCor(P, M=1000)
#did not work. We can add more points?