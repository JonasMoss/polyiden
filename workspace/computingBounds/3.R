rm(list=ls())

require(pbivnorm)
require(copula)
rho <- 0.7 #desired correlation

#Let's start by setting up a frame-work to generate distributions of X
#relative to an underlying Gaussian Z with given thresholds.

K <- 5
tau1 <- seq(-2.7,2.7, length.out=K)
tau2 <- tau1

S <- expand.grid((1:K), (1:K)) #Let us generate the indices where we have restrictions.
S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs
dim(S)[1] == K^2-1+2-2*K #check whether the number of elements are correct.

#Let us now add the left hand side of eq(3) to this matrix, then two
#new columns, which are the marginal probabilities used as inputs
#on the right of eq(3).

x1 <- pnorm(tau1[S[,1]])
x2 <- pnorm(tau2[S[,2]])
z <- pbivnorm(tau1[S[,1]],tau1[S[,2]], rho = rho)

P <- as.data.frame(cbind(S, x1,x2,z))
colnames(P) <- c("j1", "j2", "PX1", "PX2", "C")
P

plot(P$PX1, P$PX2)
abline(a=0,b=1) #four points on the diagonal. Will be relevant for the symmetry illustration.

copObj <- normalCopula(param=c(rho), dim = 2) #sanity check for the bounds

u <- 0.4
v <- 0.2

max(0,u+v-1, max(P$C - (P$PX1-u)*(P$PX1-u >= 0) - (P$PX2-v)*(P$PX2-v >= 0)))
pCopula(c(u,v), copObj)
min(u,v,min(P$C + (u-P$PX1)*(u-P$PX1 >= 0) + (v-P$PX2)*(v-P$PX2 >= 0)))
#contained
