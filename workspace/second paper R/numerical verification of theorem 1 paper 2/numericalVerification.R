rm(list=ls())
require(copula)
require(covsim) #Our VITA-implementation. Used here only in the bivariate case, and only
#in a trivial manner.
require(rvinecopulib) #To simulate from the VITA-result.

K <- 7
tau1 <- c(qnorm(seq(0.1, 0.9, length.out=K-1)), Inf)
tau2 <- tau1

S <- expand.grid((1:K), (1:K)) 

#Our central object will be the matrix P, which we now start to construct.
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]]))
colnames(P) <- c("i","j", "tau1", "tau2") 
P

upperRight <- expand.grid(tau1, tau2)

#We now construct a matrix encoding the three remaining coordinates of each of the boxes.
#(not all corners are used in the following code)
upperLeft <- upperRight*NA
lowerRight <- upperRight*NA
lowerLeft <- lowerRight*NA

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
copObjmvdc <- mvdc(normalCopula(0), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
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

getPoints <- function(P, M, eps = 0.001) { 
  #If M > 0, M is the addition to the x and y coordinate of the upper right corner point.
  #If M < 0, M(1,-1) is added to the x and y coordinate of the lower right corner point.
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
    if ( (M > 0) & (i == K) & (j==K) ) {
      xlim <- xlim + M
      ylim <- ylim + M
    }
    
    if ( (M < 0) & (i == K) & (j==1) ) {
      xlim <- xlim + abs(M)
      ylim <- ylim  - abs(M)
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
  
  XY <- as.data.frame(getPoints(P,M, eps))
  colnames(XY) <- c("x", "y")
  
  plot(NA, xlim=c(min(XY$x) - 5*eps, max(XY$x)+5*eps), ylim=c(min(XY$y) - 5*eps, max(XY$y)+5*eps), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
  
  
  for(i in (1:(K-1))) {
    abline(v=tau1[i])
    abline(h=tau2[i])
  }
  points(XY)
  
  #To draw the arrow:
  if(M > 0) {
    arrows(tau1[K-1]+eps,  tau2[K-1]+eps,  tau1[K-1]+eps+M,  tau2[K-1]+eps+M, lty=2)
  }
  if(M < 0) {
    arrows(tau1[K-1]+eps,  tau2[1]+eps,  tau1[K-1]+eps+abs(M),  tau2[1]+eps-abs(M), lty=2)
  }
  
}

drawPoints(P, M=0)
drawPoints(P, M=5)

drawPoints(P, M=0)
drawPoints(P, M=-5)

#A version which produce both "moves" in one graph:
drawPointsCombined <- function(P, M, eps = 0.001) {
  tau1 <- sort(unique(P$tau1))
  tau2 <- sort(unique(P$tau2))
  K <- length(tau1)
  
  XY <- as.data.frame(getPoints(P,0, eps)) #NB: the moves are not done, we just draw arrows
  colnames(XY) <- c("x", "y")
  
  plot(NA, xlim=c(min(XY$x) - 5*eps, max(XY$x)+5*eps + M), ylim=c(min(XY$y) - 5*eps - M, max(XY$y)+5*eps + M), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
  
  
  for(i in (1:(K-1))) {
    abline(v=tau1[i])
    abline(h=tau2[i])
  }
  points(XY)
  
  #To draw the arrow:
  if(M > 0) {
    arrows(tau1[K-1]+eps,  tau2[K-1]+eps,  tau1[K-1]+eps+M,  tau2[K-1]+eps+M, lty=2)
    arrows(tau1[K-1]+eps,  tau2[1]+eps,  tau1[K-1]+eps+M,  tau2[1]+eps-M, lty=2)
  }
}

drawPointsCombined(P,M=1)

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

calcCor(P, M=0)
calcCor(P, M=-5)
calcCor(P, M=-10)
calcCor(P, M=-100)

N <- 200
Mindx <- floor(seq(0,100, length.out=N))
rhoPos <- NULL
rhoNeg <- NULL
for (i in (1:N)) {
  rhoPos[i] <- calcCor(P, M=Mindx[i])
  rhoNeg[i] <- calcCor(P, M=-Mindx[i])
}
pdf("corrEffect.pdf")
par(mfrow=c(1,2))
drawPointsCombined(P,M=1)
plot(Mindx, rhoPos, type="l", xlab="s", ylab=expression(rho), ylim=c(-1,1))
points(Mindx, rhoNeg, type="l", lty=2, col="blue")
abline(h=1, col="red", lty=3)
abline(h=-1, col="red", lty=3)
legend("center", c("Upper", "Lower"), col=c("black","blue"), lty=c(1,2))
dev.off()
