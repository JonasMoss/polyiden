require(pbivnorm)

#Let's start by setting up a frame-work to generate distributions of X
#relative to an underlying Gaussian Z with given thresholds.

K <- 5
tau1 <- seq(-2.7,2.7, length.out=K)
tau2 <- tau1

#Let us generate the indices where we have
#restrictions.
S <- expand.grid((1:K), (1:K))
#we remove all unwanted pairs
#1) The point C(1,1)
inx <- which((S[,1] == K) & (S[,2] == K))
inx
S <- S[-inx,]
#2) the marginal distributions:
inx <- which(S[,1] == K | S[,2] == K)
inx
S <- S[-inx,]
dim(S)[1]
#Should have 
K^2-1+2-2*K
#which is a match




rho <- 0.7
pbivnorm(x,y, rho = rho)