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
intra <- array(intra[, 2], dim = c(1, nrow(intra)))
C <- array(nrow(comb))
comb <- array(comb[, 2], dim = c(1, nrow(comb)))
M <- nrow(inter)
inter <- inter[, 2]


### MLEs ###

p <- mean(intra >= min(inter))
q <- mean(inter <= max(intra))
p_prime <- mean(intra >= min(comb))
q_prime <- mean(comb <= max(intra))
log10_p <- log10(p)
log10_q <- log10(q)
log10_p_prime <- log10(p_prime)
log10_q_prime <- log10(q_prime)


### Posterior Estimates ####

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

plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  geom_hline(yintercept = q, color = "blue") +
  geom_vline(xintercept = p, color = "red") +
  ggtitle(expression(p[lwr]*" vs. "*p[upr]))

plot2 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  geom_hline(yintercept = log10(q), color = "blue") +
  geom_vline(xintercept = log10(p), color = "red") +
  ggtitle(expression(log[10](p[lwr])*" vs. "*log[10](p[upr])))

print(plot1)
print(plot2)

grid.arrange(plot1, plot2, ncol=1)

# Create a ggplot object with the histograms
p1 <- ggplot(post, aes(x = p_lwr)) +
  geom_histogram() +
  labs(x =  expression(p[lwr]), title = expression(p[lwr]))

p2 <- ggplot(post, aes(x = p_upr)) +
  geom_histogram() +
  labs(x =  expression(p[upr]), title = expression(p[upr]))

p3 <- ggplot(post, aes(x = log10_p_lwr)) +
  geom_histogram() +
  labs(x = expression(log[10](p[lwr])), title = expression(log[10](p[lwr])))

p4 <- ggplot(post, aes(x = log10_p_upr)) +
  geom_histogram() +
  labs(x = expression(log[10](p[upr])), title = expression(log[10](p[upr])))

# Arrange plots in a 2x2 grid using facet_wrap
combined_plots <- list(p1, p2, p3, p4)
names(combined_plots) <- c("p_lwr", "p_upr", "log10_p_lwr", "log10_p_upr")

grid.arrange(grobs = combined_plots, ncol = 2)



