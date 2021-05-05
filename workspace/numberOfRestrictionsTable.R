#Number of restrictions for d = 2:
require(xtable)
K <- (2:10)
X <- cbind(as.integer(K),as.integer(K^2-1+2-2*K))
colnames(X) <- c("K", "Restrictions")

X

xtable(X)
