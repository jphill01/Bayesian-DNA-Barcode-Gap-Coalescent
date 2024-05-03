##### Set working directory #####

setwd("/Users/jarrettphillips/desktop/Bayesian DNA Barcode Gap")
# setwd("/Users/jarrettphillips/desktop/Coalescent Data and R Code/COI-5P")

library(ggplot2)
library(gridExtra)

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

(p <- mean(intra >= min(inter)))
(q <- mean(inter <= max(intra)))
(p_prime <- mean(intra >= min(comb)))
(q_prime <- mean(comb <= max(intra)))
(log10_p <- log10(p))
(log10_q <- log10(q))
(log10_p_prime <- log10(p_prime))
(log10_q_prime <- log10(q_prime))

N * p
M * q
N * p_prime
C * q_prime


### Posterior Estimates ####

fit <- stan("DNA_barcode_gap.stan", 
            data = list(K = K, M = M, N = N, intra = intra, inter = inter, C = C, comb = comb),
            iter = 2000,
            control = list(adapt_delta = 0.80))

print(fit, digits_summary = 6)

traceplot(fit)


# Plot histograms of posterior samples with observed value shown

post <- as.data.frame(extract(fit))

plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  geom_vline(xintercept = p, color = "red") +
  geom_hline(yintercept = q, color = "blue") +
  ggtitle(expression(p[lwr]*" vs. "*p[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime, y = p_upr_prime)) +
  geom_point() +
  xlab(expression(p[lwr]*"'")) +
  ylab(expression(p[upr]*"'")) +
  geom_vline(xintercept = p_prime, color = "red") +
  geom_hline(yintercept = q_prime, color = "blue") +
  ggtitle(expression(p[lwr]*"'"*" vs. "*p[upr]*"'"))

plot3 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  geom_vline(xintercept = log10(p), color = "red") +
  geom_hline(yintercept = log10(q), color = "blue") +
  ggtitle(expression(log[10](p[lwr])*" vs. "*log[10](p[upr])))

plot4 <- ggplot(post, aes(x = log10_p_lwr_prime, y = log10_p_upr_prime)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  geom_vline(xintercept = log10(p_prime), color = "red") +
  geom_hline(yintercept = log10(q_prime), color = "blue") +
  ggtitle(expression(log[10](p[lwr]*"'")*" vs. "*log[10](p[upr]*"'")))


print(plot1)
print(plot2)
print(plot4)
print(plot4)

grid.arrange(plot1, plot2, plot3, plot4, ncol=1)




plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(y[lwr])) +
  ylab(expression(y[upr])) +
  geom_vline(xintercept = N * p, color = "red") +
  geom_hline(yintercept = M * q, color = "blue") +
  ggtitle(expression(y[lwr]*" vs. "*y[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime, y = p_upr_prime)) +
  geom_point() +
  xlab(expression(y[lwr]*"'")) +
  ylab(expression(y[upr]*"'")) +
  geom_vline(xintercept = N * p_prime, color = "red") +
  geom_hline(yintercept = C * q_prime, color = "blue") +
  ggtitle(expression(y[lwr]*"'"*" vs. "*y[upr]*"'"))

plot3 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](y[lwr]))) +
  ylab(expression(log[10](y[upr]))) +
  geom_vline(xintercept = N * log10(p), color = "red") +
  geom_hline(yintercept = M * log10(q), color = "blue") +
  ggtitle(expression(log[10](y[lwr])*" vs. "*log[10](y[upr])))

plot4 <- ggplot(post, aes(x = log10_p_lwr_prime, y = log10_p_upr_prime)) +
  geom_point() +
  xlab(expression(log[10](y[lwr]))) +
  ylab(expression(log[10](y[upr]))) +
  geom_vline(xintercept = N * log10(p_prime), color = "red") +
  geom_hline(yintercept = C * log10(q_prime), color = "blue") +
  ggtitle(expression(log[10](y[lwr]*"'")*" vs. "*log[10](y[upr]*"'")))


print(plot1)
print(plot2)
print(plot4)
print(plot4)

grid.arrange(plot1, plot2, plot3, plot4, ncol=1)



# Create a ggplot object with the histograms
p1 <- ggplot(post, aes(x = p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = p, color = "red") +
  labs(x =  expression(p[lwr]), title = expression(p[lwr]))

p2 <- ggplot(post, aes(x = p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = q, color = "blue") +
  labs(x =  expression(p[upr]), title = expression(p[upr]))

p3 <- ggplot(post, aes(x = p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = p_prime, color = "red") +
  labs(x =  expression(p[lwr]*"'"), title = expression(p[lwr]*"'"))

p4 <- ggplot(post, aes(x = p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = q_prime, color = "blue") +
  labs(x =  expression(p[upr]*"'"), title = expression(p[upr]*"'"))

p5 <- ggplot(post, aes(x = log10_p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = log10_p, color = "red") +
  labs(x = expression(log[10](p[lwr])), title = expression(log[10](p[lwr])))

p6 <- ggplot(post, aes(x = log10_p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = log10_q, color = "blue") +
  labs(x = expression(log[10](p[upr])), title = expression(log[10](p[upr])))

p7 <- ggplot(post, aes(x = log10_p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = log10_p_prime, color = "red") +
  labs(x = expression(log[10](p[lwr]*"'")), title = expression(log[10](p[lwr]*"'")))

p8 <- ggplot(post, aes(x = log10_p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = log10_p, color = "red") +
  labs(x = expression(log[10](p[lwr]*"'")), title = expression(log[10](p[lwr]*"'")))


# Arrange plots in a 2x2 grid using facet_wrap
combined_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime", "log10_p_lwr", "log10_p_upr", "log10_p_lwr_prime", "log10_p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)



# Create a ggplot object with the histograms
p1 <- ggplot(N * post, aes(x = p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = N * p, color = "red") +
  labs(x =  expression(y[lwr]), title = expression(y[lwr]))

p2 <- ggplot(M * post, aes(x = p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = M * q, color = "blue") +
  labs(x =  expression(y[upr]), title = expression(y[upr]))

p3 <- ggplot(N * post, aes(x = p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = N * p_prime, color = "red") +
  labs(x =  expression(y[lwr]*"'"), title = expression(y[lwr]*"'"))

p4 <- ggplot(C * post, aes(x = p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = C * q_prime, color = "blue") +
  labs(x =  expression(y[upr]*"'"), title = expression(y[upr]*"'"))

p5 <- ggplot(post, aes(x = log10_p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = 10^log10_p, color = "red") +
  labs(x = expression(log[10](y[lwr])), title = expression(log[10](y[lwr])))

p6 <- ggplot(post, aes(x = log10_p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = 10^log10_q, color = "blue") +
  labs(x = expression(log[10](y[upr])), title = expression(log[10](y[upr])))

p7 <- ggplot(post, aes(x = log10_p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = 10^log10_p_prime, color = "red") +
  labs(x = expression(log[10](y[lwr]*"'")), title = expression(log[10](y[lwr]*"'")))

p8 <- ggplot(post, aes(x = log10_p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = 10^log10_p, color = "red") +
  labs(x = expression(log[10](y[lwr]*"'")), title = expression(log[10](y[lwr]*"'")))


# Arrange plots in a 2x2 grid using facet_wrap
combined_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime", "log10_p_lwr", "log10_p_upr", "log10_p_lwr_prime", "log10_p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)


