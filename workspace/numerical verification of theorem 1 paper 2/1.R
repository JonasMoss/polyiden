#Computing the correlation via simulation
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
rho <- -0.95 #large negative correlation is chosen here, since we aim at re-arranging the table to produce
#a correlation close to 1

K <- 7

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

#We now fill up a matrix P with information.
#Currently, much of this information is not presently used,
#as I use previously written code.

u1 <- pnorm(tau1[S[,1]])
u2 <- pnorm(tau2[S[,2]])
z <- pCopula(cbind(u1,u2), claytonCopula(param=theta, dim = 2))
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]], u1,u2,z))
colnames(P) <- c("i","j", "tau1", "tau2", "U1", "U2", "CDF") 
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

plot(NA, xlim=c(-3, 3), ylim=c(-3, 3), type="n", xlab=expression(Z[1]), ylab=expression(Z[2]))
for(i in (1:(K-1))) {
  abline(v=tau1[i])
  abline(h=tau2[i])
}
for(k in (1:length(P$i) )) {
  i <- P$i[k]
  j <- P$j[k]
  if(is.finite(tau1[i])) {
    xlim <- tau1[i]-0.2 #TODO: dynamic placements of text (preferably in the middle of the box)
  } else {
    xlim <- tau1[i-1]+0.2
  }
  if(is.finite(tau1[j])) {
    ylim <- tau2[j]-0.2
  } else {
    ylim <- tau2[j-1]+0.2
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

#We will have to arrange this so that all boxes are moved to being \eps close to a straight line. The values attained
#in the "outer" boxes should be chosen close to the edges of the box.


#Let us first identify the diagonal boxes, and shrink them. 

