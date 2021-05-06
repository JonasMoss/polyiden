rm(list=ls())
set.seed(123)
require(copula)


K <- 7

prob1 <- c(seq(0.05, 0.99, length.out=K-1),1)
prob2 <- c(seq(0.01, 0.97, length.out=K-1), 1)

tau1 <- qnorm(prob1)
tau2 <- qnorm(prob2)

#prob1 <- pnorm(tau1)
#prob2 <- pnorm(tau2)
  
S <- expand.grid((1:K), (1:K)) #Let us generate the indices where we have restrictions.
#S <- S[(S[,1] != K) & (S[,2] != K),] #we remove all unwanted pairs

x1 <- prob1[S[,1]]
x2 <- prob2[S[,2]]

rho <- 0.7
norm.cop <- normalCopula(rho, dim = 2)

z <- pCopula(cbind(x1,x2), norm.cop)

P <- as.data.frame(cbind(S, x1,x2,z))
colnames(P) <- c("j1", "j2", "PX1", "PX2", "C")
P
# M <- function(u,v) {
#   res <- min(u,v,min(P$C + (u-P$PX1)*(u-P$PX1 >= 0) + (v-P$PX2)*(v-P$PX2 >= 0)))
#   return(res)
# }
# W <- function(u,v) {
#   res <- max(0,u+v-1, max(P$C - (P$PX1-u)*(P$PX1-u >= 0) - (P$PX2-v)*(P$PX2-v >= 0)))
#   return(res)
# }



u <- 0.9
v <- 0.8

hDown1 <- max(prob1[(prob1 <= u)])
hUp1 <- min(prob1[prob1 > u])

hDown2 <- max(prob2[(prob2 <= v)])
hUp2 <- min(prob2[prob2 > v])

term1 <- max(
    pCopula(c(hDown1,hDown2), norm.cop) - (hDown1 - u)*(hDown1 - u > 0),
    pCopula(c(hUp1,hDown2), norm.cop) - (hUp1 - u)*(hUp1 - u > 0)
    ) - (hDown2-v)*(hDown2-v > 0)

term2 <- max(
  pCopula(c(hDown1,hUp2), norm.cop) - (hDown1 - u)*(hDown1 - u > 0),
  pCopula(c(hUp1,hUp2), norm.cop) - (hUp1 - u)*(hUp1 - u > 0)
  ) - (hUp2-v)*(hUp2-v > 0)

max(term1, term2)
max(P$C - (P$PX1-u)*(P$PX1-u >= 0) - (P$PX2-v)*(P$PX2-v >= 0))
#match.

#wow lucky!

#Also, alternative simplification:


uTest <- hUp1
vTest <- hUp1
pCopula(c(uTest,vTest), norm.cop) - (uTest - u)*(uTest - u > 0) - (vTest-v)*(vTest-v > 0)

uTest <- hDown1
vTest <- hUp1
pCopula(c(uTest,vTest), norm.cop) - (uTest - u)*(uTest - u > 0) - (vTest-v)*(vTest-v > 0)

uTest <- hUp1
vTest <- hDown1
pCopula(c(uTest,vTest), norm.cop) - (uTest - u)*(uTest - u > 0) - (vTest-v)*(vTest-v > 0)

uTest <- hDown1
vTest <- hDown1
pCopula(c(uTest,vTest), norm.cop) - (uTest - u)*(uTest - u > 0) - (vTest-v)*(vTest-v > 0)


max(P$C - (P$PX1-u)*(P$PX1-u >= 0) - (P$PX2-v)*(P$PX2-v >= 0))

