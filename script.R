##### Set working directory #####

setwd("/Users/jarrettphillips/desktop/Bayesian DNA Barcode Gap")
# setwd("/Users/jarrettphillips/desktop/Coalescent Data and R Code/COI-5P")


##### Install required packages #####

# install.packages("rstan") # for HMC

##### Load required libraries #####

library(rstan) # for MCMC

options(mc.cores = parallel::detectCores()) # parallelize simulations
rstan_options(auto_write = TRUE) # only have to compile once unless code is changed


K <- 1
intra <- read.csv(file.choose())
inter <- read.csv(file.choose())
comb <- read.csv(file.choose())

N <- array(nrow(intra))
intra <- matrix(intra[, 2], nrow = K, ncol = N)
C <- array(nrow(comb))
comb <- matrix(comb[, 2], nrow = K, ncol = C)
M <- nrow(inter)
inter <- inter[, 2]



fit <- stan("barcode_gap.stan", 
            data = list(K = K, M = M, N = N, intra = intra, inter = inter, C = C, comb = comb),
            iter = 2000,
            control = list(adapt_delta = 0.80))
fit

traceplot(fit)


# Plot histograms of posterior samples with observed value shown

post <- as.data.frame(extract(fit))


library(ggplot2)
library(gridExtra)

par(mfrow = c(2, 1))
# First plot
plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  ggtitle(expression(p[lwr]*" vs. "*p[upr]))

# Second plot
plot2 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  ggtitle(expression(log[10](p[lwr])*" vs. "*log[10](p[upr])))

print(plot1)
print(plot2)

grid.arrange(plot1, plot2, ncol=1)




par(mfrow = c(2, 2))
hist(post$p_lwr, xlab = "p_lwr", main = expression(p[lwr]))
hist(post$p_upr, xlab = "p_upr",  main = expression(p[upr]))




