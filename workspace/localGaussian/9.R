rm(list=ls())
require(lavaan)
library(mvtnorm)
require(copula)
require(covsim) #Our VITA-implementation. Used here only in the bivariate case.
require(rvinecopulib)
require(psych)

rho <- 0.9
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
S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs
#dim(S)[1] == K^2-1+2-2*K #check whether the number of elements are correct.

u1 <- pnorm(tau1[S[,1]])
u2 <- pnorm(tau2[S[,2]])

#Sanity check:
#z <- pCopula(cbind(u1,u2), normalCopula(param=rho, dim = 2))
#gives approximately the same "local edge" correlations as the true value.

z <- pCopula(cbind(u1,u2), claytonCopula(param=theta, dim = 2))
P <- as.data.frame(cbind(S, tau1[S[,1]], tau1[S[,2]], u1,u2,z))
colnames(P) <- c("i","j", "tau1", "tau2", "U1", "U2", "CDF")
P


#Let's find the gaussian copula which reproduce the CDF at each point.
for(i in (1:length(P$CDF))) {
  
  diffCDF = Vectorize(function(rhoCand) {
    res <- P$CDF[i]-pCopula(c(P$U1[i],P$U2[i]), normalCopula(param=rhoCand, dim = 2))
    return(res)
  })
  
  P$rho[i] <- uniroot(diffCDF, interval=c(-1,1))$root
}

P

#Visualize the results
copObjmvdc <- mvdc(claytonCopula(theta), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
copObjmvdcGaus <- mvdc(normalCopula(rho), c("norm", "norm"),list(list(mean = 0, sd =1), list(mean = 0,sd=1)))
contour(copObjmvdc, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), drawlabels=F)
contour(copObjmvdcGaus, dMvdc, xlim=c(-3, 3), ylim=c(-3, 3), add=T, col="red", drawlabels=F, lty=3)
for(i in (1:(K-1))) {
  abline(v=tau1[i])
  abline(h=tau2[i])
}
for(k in (1:length(P$i) )) {
  i <- P$i[k]
  j <- P$j[k]
  text(tau1[i], tau2[j], round(P$rho[k], 2), cex=1)
}

mean(P$rho)


#Now I need to compute the (population) polychoric correlation as well.
#Should have old code for this. We only need the table.

#Code from earlier:
th1 <- c(-Inf,tau1); th2 <- c(-Inf,tau2)
my.cop <- copObjmvdc #NB, has to be the data generating copula.

my.tab <- matrix(NA, 3,3)
for (row in 1:3){
  for(col in 1:3){
    my.tab[row, col] <- pMvdc(c(th1[row+1], th2[col+1]), my.cop)-
      pMvdc(c(th1[row], th2[col+1]), my.cop)-pMvdc(c(th1[row+1], th2[col]), my.cop)+
      pMvdc(c(th1[row], th2[col]), my.cop)
  }
}
polychoric(x=as.table(my.tab))$rho


