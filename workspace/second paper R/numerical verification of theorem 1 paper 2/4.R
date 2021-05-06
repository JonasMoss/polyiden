#Let us now emulate the infinite case, having many categories, but most of them have
#low probabilities near the edge

rm(list=ls())
require(lavaan)
library(mvtnorm)
require(copula)
require(covsim) #Our VITA-implementation. Used here only in the bivariate case.
require(rvinecopulib)
require(psych)
require(rootSolve)

#The thresholds are for this illustration of central importance.
K.orig <- 7
K <- K.orig + 20
tau1 <- c(qnorm(c(
  seq(0.1, 0.9, length.out=K.orig-1), seq(0.901, 0.99, length.out=K-K.orig)
  )), Inf)
tau2 <- tau1


rho <- -0.95


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



#NAMES OF THE BOXES:
P$BoxNr <- (1:length(P$i))


drawGrid <- function(P, textDraw = FALSE) {
  plot(NA, xlim=c(-3, 3), ylim=c(-3, 3), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
  for(i in (1:(K-1))) {
    abline(v=tau1[i])
    abline(h=tau2[i])
  }
  
  if(textDraw == T) {
    for(k in (1:length(P$i) )) {
      i <- P$i[k]
      j <- P$j[k]
      if(is.finite(tau1[i])) {
        xlim <- tau1[i]-0.2 #TODO: dynamic placements of text (preferably in the middle of the box)
      } else {
        xlim <- tau1[i-1]+0.2
      }
      if(is.finite(tau2[j])) {
        ylim <- tau2[j]-0.2
      } else {
        ylim <- tau2[j-1]+0.2
      }
      text(xlim, ylim,  k, cex=1)
    }
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

#We will have to arrange this so that all boxes are moved to being \eps close the straight line y = x. The values attained
#in the "outer" boxes should be chosen close to the edges of the box.


#Let us first assign values to the boxes, starting in the middle of the boxes in case where we do not
#have an "edge case", and at value which later on will be close to the y = x line in the "edge cases".

getPoints <- function(P, M, chosenPoint, eps = 0.001) { #M is the addition to the x and y coordinate of the upper right corner-point.
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

    #special case for the upper right point:
    if ( (i == K) & (j==K) ) {
      xlim <- xlim + M
      ylim <- ylim + M
    }
    
    xpoints <- c(xpoints, xlim)
    ypoints <- c(ypoints, ylim)
  }
  return(cbind(xpoints, ypoints))
}

#We modify the drawGrid-code
drawPoints <- function(P, M, eps = 0.001) {
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  plot(NA, xlim=c(tau1[1] - 5*eps, tau1[K-1]+M+5*eps), ylim=c(tau1[1]-5*eps, tau2[K-1]+M+5*eps), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
  for(i in (1:(K-1))) {
    abline(v=tau1[i])
    abline(h=tau2[i])
  }
  
  XY <- getPoints(P,M, eps)
  points(XY)
}

drawPoints(P, M=0)
drawPoints(P, M=5)
drawPoints(P, M=10)


#Should definitely be able to ignore all points outside a large square (using the same structure as here
#for those that are within the square), just going far enough.


#We then use this structure to directly compute the correlation implied by P.
calcCor <- function(P,M, eps = 0.001) {
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  XY <- as.data.frame(getPoints(P,M,eps))
  colnames(XY) <- c("x","y")
  
  Xexp <- 0
  Xexpsq <- 0
  Yexp <- 0
  Yexpsq <- 0
  
  #Calculate E X and E X^2, and use this to get the standard deviation
  for (xval in unique(XY$x)) {
    inx.x <- which(XY$x == xval)
    Xexp <- Xexp + xval*sum(P$ProbBox[inx.x])
    Xexpsq <- Xexpsq + (xval^2)*sum(P$ProbBox[inx.x])
  }
  stdX <- sqrt( Xexpsq - (Xexp)^2 )
    
  #Same with Y:
  for (yval in unique(XY$y)) {
    inx.y <- which(XY$y == yval)
    Yexp <- Yexp + yval*sum(P$ProbBox[inx.y])
    Yexpsq <- Yexpsq + (yval^2)*sum(P$ProbBox[inx.y])
  }
  stdY <- sqrt( Yexpsq - (Yexp)^2 )
  
  #The covariance, calculated vectorially
  covXY <- sum((XY$x - Xexp)*(XY$y - Yexp)*P$ProbBox)
  rhoXY <- covXY/(stdX*stdY)
  return(rhoXY)
}

calcCor(P, M=0)
calcCor(P, M=5)
calcCor(P, M=10)
calcCor(P, M=100)


#Realization: If we have an infinite number of categories in only one direction, we can do this either
#by having one point in the (1,1) box that goes to (-infty, -infty) (which gives positive correlation)
#or, more complicated, having that (1,1) box increase in volume, having correlation "locally" equal to -1.