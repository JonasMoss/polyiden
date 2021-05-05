
#We match up a Gaussian distribution with standard marginals to each square
#In the corners, there is just one solution, as it is a "tetrachoric"-type solution.
#In the remaining cubes, there may be (and often is) multiple solutions
#for matching up a cube-probability.
#Let us choose the solution which is closest to a "tetrachoric" solution
#connected to that cube.
#
#For each box, except the corners, there are at most four ways to choose
#a "tetrachoric"-like solution. Namely looking at the CDF, the surv CDF,
#and (>, <) or (<,>). Each of these possibilities may match up with different
#gaussian distributions, though for each possibility, we have a unique solution,
#since they are "polychoric"-like.
#Let us choose to work with the direction which has the largest probability
#in total.
#
#Under normality, all "tetrachoric" solution will be constant, and the true
#underlying correlation will be one of the candidates per box. Therefore,
#choosing the solution which is closest to the "tetrachoric" will always
#identify a true solution if there is one.

#This max-direction tetrachoric correlation has the advantage of
#being (exactly) invariant to "rotations", a property shared
#by the "box probability matching" local polychorics, meaning
#that the resulting method combining the two shares this
#property.

#Therefore, the method:
#1) Does not depend on the encoding of the variables, is an ordinal method.
#2) Is consistent under normality.
#3) Is invariant under "rotation".

#Note that the "max-direction tetrachoric" correlations seems not to
#correspond to properties of latent variables which are discretized,
#and this is the motivation for the box-focus.


#Computationally: I pluss boxes, which are already computed!

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
  text(xlim, ylim,  k, cex=1)
}







#Fill in box-probabilities in P
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

tetrachoricBox <- function(k.ind) { #Finds the distributions which match up with the boxes with
  #index in k.ind and probability mass in P$box
  probDiff = Vectorize(function(rhoCand) {
    #print(rhoCand)
    copObjmvdcGaus <- mvdc(normalCopula(rhoCand), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))  
    probGauss <- 0
    for(k in k.ind) {
      x1 <- upperRight[k,1]
      y1 <- upperRight[k,2]
      x2 <- lowerLeft[k,1]
      y2 <- lowerLeft[k,2]
      probGauss <- probGauss +
        pMvdc(c(x1,y1), copObjmvdcGaus) +
        pMvdc(c(x2, y2), copObjmvdcGaus) -
        pMvdc(c(x2, y1), copObjmvdcGaus) -
        pMvdc(c(x1, y2), copObjmvdcGaus)
    }
    res <- sum(P$Box[k.ind]) - probGauss
    return(res)
  })
  epsInterval <-  c(-1 + (2*.Machine$double.eps) ,1- (2*.Machine$double.eps))  
  print(k.ind)
  plot(seq(-0.999, 0.999, by = 0.01), probDiff(seq(-0.999, 0.999, by = 0.01)), type = "l", xlim=c(-1,1), main=paste("k=", k.ind,sep="", collapse = ","))
  abline(h=0)
  return(uniroot.all(probDiff,epsInterval))
}


PRoots <- NULL
PRoots <- as.data.frame(PRoots)
for(k in (1:dim(S)[1])) {
  print(k)
  
  rootsAll <- tetrachoricBox(k)
  distSol <- 0
  
  if(length(rootsAll) == 0) { #no "exact" solution, so we use the best available.
    
    #TODO: the following should not be a separate function. Implement one single
    #optimization routine which returns all that is needed
    
    probDiff = Vectorize(function(rhoCand) { #fixed k version
      copObjmvdcGaus <- mvdc(normalCopula(rhoCand), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))  
      probGauss <- 0
      x1 <- upperRight[k,1]
      y1 <- upperRight[k,2]
      x2 <- lowerLeft[k,1]
      y2 <- lowerLeft[k,2]
      probGauss <- probGauss +
        pMvdc(c(x1,y1), copObjmvdcGaus) +
        pMvdc(c(x2, y2), copObjmvdcGaus) -
        pMvdc(c(x2, y1), copObjmvdcGaus) -
        pMvdc(c(x1, y2), copObjmvdcGaus)
      res <- sum(P$Box[k]) - probGauss
      return(res)
    })
    
    
    
    epsInterval <-  c(-1 + (2*.Machine$double.eps) ,1- (2*.Machine$double.eps))  
    optRes <- optimize(function(rhoCand) { abs(probDiff(rhoCand)) }, epsInterval)
    distSol <- optRes$objective
    
    #Is the minimum unique?
    rootsAll2 <- uniroot.all(function(rhoCand) { abs(probDiff(rhoCand)) - optRes$objective}, c(-1,1))
    if(length(rootsAll2) > 1) {
      rootsAll <- rootsAll2
    } else {
      rootsAll <- optRes$minimum
    }
  }
  
  #We now choose the solution closest to the "half-local tetrachoric like" solution.
  
  best.inx <- 1 
  
  #NB: all inequalities are inclusive, so we include the box
  #at hand.
  #S[k,]
  #1) North east tetrachorics, corresponding to survival CDFs
  k.ind.ne <- which( (S[,1] >= S[k,1]) & (S[,2] >= S[k,2]) )
  ne.prob <- sum(P$Box[k.ind.ne])
  
  
  #2) North west tetrachorics
  k.ind.nw <- which( (S[,1] <= S[k,1]) & (S[,2] >= S[k,2]) )
  nw.prob <- sum(P$Box[k.ind.nw])
  
  
  #3) South east 
  k.ind.se <- which( (S[,1] >= S[k,1]) & (S[,2] <= S[k,2]) )
  se.prob <- sum(P$Box[k.ind.se])
  
  
  
  #4) South west (standard tetrachorics based on the CDF)
  k.ind.sw <- which( (S[,1] <= S[k,1]) & (S[,2] <= S[k,2]) )
  sw.prob <- sum(P$Box[k.ind.sw])
  
  #We discard solutions which have an entire column or row in it,
  #as it will be deducable from marignal knowledge, hence
  #having no variation in \rho.
  
  if ( (length(unique(S[k.ind.ne,1])) == K) |  (length(unique(S[k.ind.ne,2])) == K) ) {
    ne.prob <- 0
  }
  if ( (length(unique(S[k.ind.nw,1])) == K) |  (length(unique(S[k.ind.nw,2])) == K) ) {
    nw.prob <- 0
  }
  if ( (length(unique(S[k.ind.se,1])) == K) |  (length(unique(S[k.ind.se,2])) == K) ) {
    se.prob <- 0
  }
  if ( (length(unique(S[k.ind.sw,1])) == K) |  (length(unique(S[k.ind.sw,2])) == K) ) {
    sw.prob <- 0
  }
  
  #we also discard trivial soltuions, e.g. 0 or 1.
  gathered <- c(ne.prob, nw.prob, se.prob, sw.prob)
  omit.inx <- which((gathered <= 10*.Machine$double.eps) | (gathered >= 1-(10*.Machine$double.eps)))
  gathered[omit.inx] <- NA
  if(sum(is.na(gathered)) == 4) {
    stop("Invalid table")
  }
  
  tetrachoric <- NA
  #NB, which.max always returns just one element, so this will work also if max is not unique (it will pick the first)
  if(which.max(gathered) == 1) {
    #NE
    tetrachoric <- tetrachoricBox(k.ind.ne)
  } else if(which.max(gathered) == 2) {
    #NW 
    tetrachoric <- tetrachoricBox(k.ind.nw)
  } else if(which.max(gathered) == 3) {
    #SE 
    tetrachoric <- tetrachoricBox(k.ind.se)
  } else {
    #SW
    tetrachoric <- tetrachoricBox(k.ind.sw)
  }
  
  cat("Found following candidates: ")  
  print(rootsAll)
  best.inx <- which.min(abs(rootsAll - tetrachoric))
  
  P$rhoTetrachoric[k] <- tetrachoric
  P$rhoLocal[k] <- rootsAll[best.inx]
  P$rhoNonExact[k] <- distSol
  
  
  #let us also save all correlation candidates.
  for(i in (1:length(rootsAll))) {
    PRoots <- rbind(PRoots, c(k, rootsAll[i], distSol))
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
