#In each square, we use the solution which is closest to the "tetrachoric"-type solution.


#Under normality, there is one solution which is constant.
#Let us choose the solution which is closest to the "tetrachoric" solution.
#Under normality, the "tetrachoric" solution will be constant, and will result
#in consistency under normality.


rm(list=ls())
require(lavaan)
library(mvtnorm)
require(copula)
require(covsim) #Our VITA-implementation. Used here only in the bivariate case.
require(rvinecopulib)
require(psych)
require(rootSolve)

#Currenlty: only positive correlations, since otherwise we must rotate. TODO: Fix this.
#rho <- 0.6
rho <- 0.6

K <- 5


#Step 1; find the copula-parameter which matches up with the desired correlation.
sigma.target <- matrix(data=c(1,rho,rho,1), byrow=TRUE, ncol=2,nrow=2)
marginsnorm <- lapply(X=sqrt(diag(sigma.target)),function(X) list(distr="norm", sd=X) )
family_set <- "clayton"

set.seed(1234)
cv <- vita(marginsnorm, sigma.target =sigma.target, family_set = family_set)

#extracting the parameters is easier said than done
cvCop <- cv$copula$pair_copulas
theta <- as.numeric(cvCop[[1]][[1]]$parameters)

tau1 <- c(qnorm(seq(0.1, 0.9, length.out=K-1)), Inf)
tau2 <- c(qnorm(seq(0.05, 0.97, length.out=K-1)), Inf)
S <- expand.grid((1:K), (1:K)) 


u1 <- pnorm(tau1[S[,1]])
u2 <- pnorm(tau2[S[,2]])
z <- pCopula(cbind(u1,u2), claytonCopula(param=theta, dim = 2))
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]], u1,u2,z))
colnames(P) <- c("i","j", "tau1", "tau2", "U1", "U2", "CDF")
P


#Each integer pair in S corresponds to the upper right corner of a box
upperRight <- expand.grid(tau1, tau2)

#We now construct a matrix encoding the three remaining coordinates of each of the boxes
upperLeft <- upperRight*NA
lowerRight <- upperRight*NA
lowerLeft <- lowerRight*NA

#Seems like there is a clever way to generate these automatically in the same way as
#upperRight, through expand.grid and reversing the elements one or both of tau1 tau2.
#I'll go for the simpler route here:
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

#Visual check:
#Visualize the results
#copObjmvdc <- mvdc(claytonCopula(theta), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
#copObjmvdcGaus <- mvdc(normalCopula(rho), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
#contour(copObjmvdc, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), drawlabels=F)
#contour(copObjmvdcGaus, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), add=T, col="red", drawlabels=F, lty=3)
# 
# plot(NA, xlim=c(-3,3), ylim=c(-3,3), type="n")
# for(i in (1:(K-1))) {
#   abline(v=tau1[i])
#   abline(h=tau2[i])
# }
# abline(a=0,b=1)
# 
# set.seed(123)
# for(k in (1:dim(S)[1])) {
#   eps <- runif(1,min=-0.01, max=0.01)
#   points(upperRight[k,1], upperRight[k,2] + eps, col=k) 
#   points(lowerRight[k,1] + eps, lowerRight[k,2] + eps, col=k)
#   points(lowerLeft[k,1] + eps, lowerLeft[k,2] + eps, col=k)
#   points(upperLeft[k,1] + eps, upperLeft[k,2] + eps, col=k)
# }
#
# k <- 8
# x1 <- upperRight[k,1]
# y1 <- upperRight[k,2]
# 
# x2 <- lowerLeft[k,1]
# y2 <- lowerLeft[k,2]
# 
# for(k in (1:dim(S)[1])) {
#   xleft=lowerLeft[k,1]
#   ybottom=lowerLeft[k,2]
#   xright=upperRight[k,1]
#   ytop=upperRight[k,2]
#   rect(xleft, ybottom, xright, ytop, col="white", border="red")
# }


#fill in box-probabilities in P
copObjmvdc <- mvdc(claytonCopula(theta), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
P$Box <- NA
for(k in (1:dim(S)[1])) {
  
  x1 <- upperRight[k,1]
  y1 <- upperRight[k,2]
  x2 <- lowerLeft[k,1]
  y2 <- lowerLeft[k,2]
  
  P$Box[k] <- pMvdc(c(x1,y1), copObjmvdc) +
    pMvdc(c(x2, y2), copObjmvdc) -
    pMvdc(c(x2, y1), copObjmvdc) -
    pMvdc(c(x1, y2), copObjmvdc)
}
P

#For the upper block, we need the survival probabilities in order to
#calculate "reversed" tetrachorics.
P$surv <- 1 - P$U1 - P$U2 + P$CDF
round(P, 4)



PRoots <- NULL
PRoots <- as.data.frame(PRoots)
for(k in (1:dim(S)[1])) {
  print(k)
  
  x1 <- upperRight[k,1]
  y1 <- upperRight[k,2]
  
  x2 <- lowerLeft[k,1]
  y2 <- lowerLeft[k,2]
  
  probDiff = Vectorize(function(rhoCand) {
    copObjmvdcGaus <- mvdc(normalCopula(rhoCand), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))  
    probGauss <- pMvdc(c(x1,y1), copObjmvdcGaus) +
      pMvdc(c(x2, y2), copObjmvdcGaus) -
      pMvdc(c(x2, y1), copObjmvdcGaus) -
      pMvdc(c(x1, y2), copObjmvdcGaus)
    res <- P$Box[k] - probGauss
    return(res)
  })
  
  epsInterval <-  c(-1 + (2*.Machine$double.eps) ,1- (2*.Machine$double.eps))  
  rootsAll <- uniroot.all(probDiff,epsInterval)
  distNA <- 0
  
  if(length(rootsAll) == 0) { #no "exact" solution
    optRes <- optimize(function(rhoCand) { abs(probDiff(rhoCand)) }, epsInterval)
    distNA <- optRes$objective
    
    #Is the minimum unique?
    rootsAll2 <- uniroot.all(function(rhoCand) { abs(probDiff(rhoCand)) - optRes$objective}, c(-1,1))
    if(length(rootsAll2) > 1) {
      rootsAll <- rootsAll2
    } else {
      rootsAll <- optRes$minimum
    }
  }
  
  #We now choose the solution closest to the "half-local tetrachoric like" solution.
  
  best.inx <- 1 #PRESENTLY ONLY FOR MOST INDICES
  
  if((S[k,1] != K) & (S[k,2] != K)) {
    
    diffCDF = Vectorize(function(rhoCand) {
      res <- P$CDF[k]-pCopula(c(pnorm(tau1[S[k,1]]),pnorm(tau2[S[k,2]])), normalCopula(param=rhoCand, dim = 2))
      return(res)
    })
    
    tetrachoric <- uniroot(diffCDF, interval=epsInterval)$root
    best.inx <- which.min(abs(rootsAll - tetrachoric))
  } else if ((S[k,1] == K) & (S[k,2] = 1)) {
    #for the top left corner, we need something else (than the next case), since in terms of its survival
    #CDF, this corresponds to the right bottom corner of the standard case,
    #which cannot be found through the CDF
    
    
    
  } else {
    #For the top most row and right most column, we match survival copulas instead.
    #NB: The thresholds corresponds to the upper right corners of the boxes.
    #We want the table flipped, meaning we need the lower left corner of each
    #box instead of the upper right. E.g., in the Inf, Inf-box, the CDF is 1,
    #and the SurvCDF is 0. We instead want to match the survCDF at tau1[K-1], tau2[K-1].
    diffCDFSurv = Vectorize(function(rhoCand) {
      
      u1 <- pnorm(tau1[S[k,1]-1])
      u2 <- pnorm(tau2[S[k,2]-1])
      survNorm <- 1 - u1 - u2 + pCopula(cbind(u1,u2), normalCopula(param=rhoCand, dim = 2))
      
      return(P$surv[k]-survNorm)
    })
    
    plot(seq(-1, 1, by = 0.01), diffCDFSurv(seq(-1, 1, by = 0.01)), type = "l", xlim=c(-1,1), main=paste("k=", k,sep=""))
    abline(h=0)
    
    tetrachoricSurv <- uniroot(diffCDFSurv, interval=epsInterval)$root
    best.inx <- which.min(abs(rootsAll - tetrachoricSurv))
  }
  
  
  
  
  P$rhoTetrachoric[k] <- tetrachoric
  P$rhoLocal[k] <- rootsAll[best.inx]
  P$rhoNonExact[k] <- distNA
  
  #    testplot <- function(rhoCand) { abs(probDiff(rhoCand)) - optRes$objective  }
  #    plot(seq(-1, 1, by = 0.01), testplot(seq(-1, 1, by = 0.01)), type = "l", xlim=c(-1,1), main=paste("k=", k,sep=""))
  #    abline(v=optRes$minimum)
  #    abline(h=0)
  
  
  for(i in (1:length(rootsAll))) {
    PRoots <- rbind(PRoots, c(k, rootsAll[i], distNA))
  }
  
  
  #plot(seq(-1, 1, by = 0.01), probDiff(seq(-1, 1, by = 0.01)), type = "l", xlim=c(-1,1), main=paste("k=", k,sep=""))
  #abline(h=0)
}

colnames(PRoots) <- c("k", "rho", "Dist") #Dist = abs diff at minimum, > 0 if the solution is 
round(PRoots, 3)





plot(NA, xlim=c(-3, 3), ylim=c(-3, 3), type="n")
for(i in (1:(K-1))) {
  abline(v=tau1[i])
  abline(h=tau2[i])
}

for(k in (1:length(P$i) )) {
  i <- P$i[k]
  j <- P$j[k]
  
  if(is.finite(tau1[i])) {
    xlim <- tau1[i]-0.4 #TODO: dynamic placements of text (preferably in the middle of the box)
  } else {
    xlim <- tau1[i-1]+0.4
  }
  
  if(is.finite(tau1[j])) {
    ylim <- tau2[j]-0.4
  } else {
    ylim <- tau2[j-1]+0.4
  }
  
  text(xlim, ylim,  round(P$rhoLocal[k],2), cex=1)
  
  
  if(P$rhoNonExact[k] > 0) {
    #special care for the infinite rectangles:
    xleft=max(lowerLeft[k,1], -3)
    ybottom=max(lowerLeft[k,2], -3)
    xright=min(upperRight[k,1], 3)
    ytop=min(upperRight[k,2], 3)
    rect(xleft, ybottom, xright, ytop, col="red", border="blue", density=10)  
  }
}


round(P, 3)
