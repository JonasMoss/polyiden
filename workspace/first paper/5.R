#A systematic exploration of the widths of the intervals.
#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("utility.R")

computeBounds <- function(tau1,tau2) {
  res <- NULL
  for(rho in seq(-0.95,0.95, length.out=100)) {
    a <- pnorm(tau1)
    b <- pnorm(tau2)
    p00 <- copula::pCopula(copula = copula::normalCopula(rho), u=c(a,b))

    p01 <- a - p00
    p10 <- b - p00
    p11 <- 1 - (p00 + p01 + p10)
    p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)

    Nbounds <- tetrachoric(p)
    Ubounds <- lspearman(p)
    res <- rbind(res, c(rho, Nbounds, Nbounds[2] - Nbounds[1],
                        Ubounds, Ubounds[2] - Ubounds[1]))
  }

  colnames(res) <- c("rho", "Nlower", "Nupper", "Nlength",
                     "Ulower", "Uupper", "Ulength")
  res <- data.frame(res)

  par(las = 1)
  plot(x = res$rho,
       y = res$Nlength,
       ylim = c(0,2),
       type="l",
       ylab="Interval length",
       xlab = expression(rho),
       main = bquote(tau[1] ~ " = " ~ .(tau1) ~ "and" ~
                       tau[2] ~ " = " ~ .(tau2) ),
       lwd = 2)
  lines(x = res$rho, y = res$Ulength, lty = 2, lwd = 2)
  legend("bottom", c("Normal marginals", "Uniform marginals"), lty = c(1,2),
         bty = "n", cex = 1.4)
  abline(v = c(-0.95, 0.95), lty = 3)
  axis(side = 3, c(-0.95, 0.95))


}


computeBounds(-0.6, 0.6)
computeBounds(-0.6, 1.2)
computeBounds(0, 0)
computeBounds(0, 1.2)
computeBounds(1.2, 1.2)

