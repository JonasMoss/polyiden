a = seq(-0.88, 0, by = 0.01)
b = seq(0, 0.93, by = 0.01)

pdf("congeneric.pdf", height = 8, width = 8)
mar = c(5.1, 4.1, 4.1, 2.1) 
par(mar = mar + c(0, 0.5, 0, 0))
plot(c(-sqrt(abs(a)),sqrt(abs(b))), sqrt(abs(c(a,b))), type = "l", 
     ylim = c(0, 1), lwd = 2, lty = 1, xlab = expression(lambda[2]),
     ylab = expression(lambda[1]~and~sigma), cex.lab = 1.4, xlim = c(-1, 1))
lines(c(-sqrt(abs(a)),sqrt(abs(b))), sqrt(1 - c(a,b)^2), lwd = 2, lty = 2)
axis(side = 3, at = c(-sqrt(0.85)),
     labels = expression(sqrt(0.85)))
axis(side = 3, at = c(sqrt(0.99)),
     labels = expression(sqrt(0.99)))
legend("bottomleft", c(expression(lambda[1]), expression(sigma)),
       lwd = c(2, 2), lty = c(1, 2), bty = "n", cex = 2)
dev.off()

#Values under normality:
rho <- 0.15
lambda1 <- sqrt(abs(rho))
round(lambda1,2)
lambda2 <- sign(rho)*sqrt(abs(rho))
round(lambda2,2)
sigma <- sqrt(1-rho^2)
round(sigma, 2)
