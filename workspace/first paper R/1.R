#Numerical illustration for Section 2.3 and Section 4.
#The numerical illustration for section 2.4 is found in 4.R

#NB: Change to your directory:
pathOfDirectory <- "C:/Users/A1110339/Dropbox/research/Njål Jonas Steffen/ordinal paper/3.01/attachments/numerical illustrations"
setwd(pathOfDirectory)
rm(list=ls())
source("utility.R")

#1.
#The numerical illustrations in Section 2.3. The numerical
#example is the same as the numerical illustration in Nelson (2017),
#Theorem 3.2.3, p.71. The extremal copulas are illustrated on the same
#page, see his Figure 3.10. These will be illustrated in 2.R.

p00 <- 0.2
a <- 0.6
p01 <- a - p00
p01
#0.4
b <- 0.3
p10 <- b - p00
p10
#0.1
p11 <- 1 - (p00 + p01 + p10)
p <- matrix(c(p00,p01, p10, p11), byrow=T, ncol=2)
p
#     [,1] [,2]
#[1,]  0.2  0.4
#[2,]  0.1  0.3


#Proposition 1 with normal marginals:
round(tetrachoric(p),2) #-0.88  0.93 

#Proposition 2, Spearman's rho with unknown Z-distribution:
round(lspearman(p),2) #-0.82  0.88 

#For comparison, Pearson's tetrachoric correlation, assuming normality:
rhohat <- polycor::polychor(p)
round(rhohat,2) #0.15

#For comparison, Spearman's rho when assuming that the copula
#is normal:
res <- spearmanGaussCop(p)
round(res, 2) #[1] 0.14

# Propostion 1 with skew t margials marginals


### ============================================================================
### Skew t examples
### ============================================================================

#' Variance function for skew t distributions.
#' @param omega,alpha,nu The parameters from sn::dskt.
#' @return The variance if the distribution.

vskt = function(omega, alpha, nu) {
  
  bv = function(nu) {
    sqrt(nu) * gamma(0.5*(nu - 1)) / (sqrt(pi) * gamma(1/2 * nu))
  }
  
  delta = function(alpha) alpha / sqrt(1 + alpha^2)
}

omega^2 * (nu / (nu - 2) - bv(nu) ^ 2 * delta(alpha) ^ 2)

}

xi = 0 # Centraliy parameter
omega = 1 # Dispersion: Equals sigma in a normal when nu = Inf and alpha = 0.
alpha = 100 # Skewness parameter
nu = 3 # Fat-tailedness, df from t distribution.

bound_length = Vectorize(function(alpha, nu) {
  sigma = sqrt(vskt(omega, alpha, nu))
  dist = function(x) sn::pst(x * sigma, xi = xi, omega = omega, alpha = alpha, nu = nu)
  quant = function(x) sn::qst(x, xi = xi, omega = omega, alpha = alpha, nu = nu) / sigma
  bound = lcorr(p, dist = dist, quant = quant, var = 1)
  unname(bound[2] - bound[1])
})

## Make the plot.

alphas = seq(-20, 20, by = 1)
nus = seq(3, 100, by = 5)
bound_lengths = outer(alphas, nus, bound_length)
rownames(bound_lengths) = alphas
colnames(bound_lengths) = nus 
bound_lengths = t(bound_lengths)
theme_set(theme_classic())
as.data.frame(bound_lengths) %>% 
  rownames_to_column(var="row") %>% 
  gather(col, value, -row) %>% 
  mutate(row=as.numeric(row), 
         col=as.numeric(col)) %>% 
  ggplot(aes(col, row, z = value)) +
  geom_contour_filled(bins = 10) + 
  xlab(expression(alpha)) + 
  ylab(expression(nu)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=16,face="bold")) + 
  guides(fill=guide_legend(title="Interval length")) ->
  plot_lengths

ggsave("../lengths.pdf", plot_lengths, width = 10, height = 8, dpi = "retina")

min(bound_lengths) # 1.756163
max(bound_lengths) # 1.946824

saveRDS(object = bound_lengths, file = "bound_lengths.Rds")

