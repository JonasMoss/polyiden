## This file plots the regions of continuity for the piece-wise contiouns 
## function C_U(u,v). See the appendix, A.4

a <- 0.1
b <- 0.3
c <- 0.4

f <- function(x) {
  u <- x[1]
  v <- x[2]
  h <- which.min(c(u, v, a + max(0, u - b) + max(0, v - c)))

  if (h == 3) {
    if (u > b) h <- 5
    if (v > c) h <- 6
    if (u > b & v > c) h <- 4
  }

  h
}

x <- seq(0, 1, by = 0.02)
xy <- expand.grid(x, x)
z <- apply(xy, 1, f)
ha <- cbind(xy, z)

par(las = 1)
par(mar = c(5,5,4,2) + 0.1)

plot(1,
     type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n",
     xlim = c(0, 1), ylim = c(0, 1)
)

for (i in 1:nrow(ha)) {
  row <- ha[i, ]
  points(row[1], row[2], col = as.character(row[3] + 1), pch = 20)
}

axis(side = 1, c(a, b, c + b - a), labels = c("a", "b", "c + b - a"))
axis(side = 2, c(a, c, c + b - a), labels = c("a", "c", "c + b - a"))
abline(h = c(0, 1, a, c, c + b - a), v = c(0, 1, a, b, c + b - a), lwd = 2)
abline(a = a - b, b = 1, lwd = 2)
abline(a = c - a, b = 1, lwd = 2)

