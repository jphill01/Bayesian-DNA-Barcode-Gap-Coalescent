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


K <- 2 # number of species in genus 

intra1 <- read.csv(file.choose())
comb1 <- read.csv(file.choose())

intra2 <- read.csv(file.choose())
comb2 <- read.csv(file.choose())

inter <- read.csv(file.choose())

N <- as.numeric(array(c(nrow(intra1), nrow(intra2))))
C <- as.numeric(array(c(nrow(comb1), nrow(comb2))))

intra <- c(intra1$x, intra2$x)
comb <- c(comb1$x, comb2$x)

M <- nrow(inter)
inter <- inter[, 2]




### MLEs ###

# probabilities by species

(p <- mean(intra1$x >= min(inter)))
(q <- mean(inter <= max(intra1$x)))
(p_prime <- mean(intra1$x >= min(comb1$x)))
(q_prime <- mean(comb1$x <= max(intra1$x)))

(p <- mean(intra2$x >= min(inter)))
(q <- mean(inter <= max(intra2$x)))
(p_prime <- mean(intra2$x >= min(comb2$x)))
(q_prime <- mean(comb2$x <= max(intra2$x)))



# ECDFs

df1 <- data.frame(
  x = intra1$x
)

df2 <- data.frame(
  x = intra2$x
)

df3 <- data.frame(
  x = inter
)

df4 <- data.frame(
  x = comb1$x
)

df5 <- data.frame(
  x = comb2$x
)

p1 <- ggplot(df1, aes(x)) + stat_ecdf() + xlab(expression(d[ij])) + geom_vline(xintercept = min(df1), color = "red") 
p2 <- ggplot(df2, aes(x)) + stat_ecdf() + xlab(expression(d[ij])) + geom_vline(xintercept = min(df2), color = "red") 
p3 <- ggplot(df3, aes(x)) + stat_ecdf() + xlab(expression(d[XY])) 
p4 <- ggplot(df4, aes(x)) + stat_ecdf() + xlab(expression(d*"'"[XY])) + geom_vline(xintercept = max(df1), color = "red") 
p5 <- ggplot(df5, aes(x)) + stat_ecdf() + xlab(expression(d*"'"[XY])) + geom_vline(xintercept = max(df2), color = "red") 


grid.arrange(p1, p2, p3, p4, p5)


# SEs

(SE_p <- sqrt(p * (1 - p) / length(intra)))
(SE_q <- sqrt(q * (1 - q) / length(inter)))
(SE_p_prime <- sqrt(p_prime * (1 - p_prime) / length(intra)))
(SE_q_prime <- sqrt(q_prime * (1 - q_prime) / length(comb)))


# counts

N * p
M * q
N * p_prime
C * q_prime


### Posterior Estimates ####

fit <- stan("DNA_barcode_gap.stan", 
            data = list(K = K, M = M, N = N, intra = intra, inter = inter, C = C, comb = comb), 
            chains = 4,
            iter = 2000,
            algorithm = "NUTS",
            control = list(adapt_delta = 0.80,
                           max_treedepth = 10))

print(fit, digits_summary = 6)

traceplot(fit, pars = c("p_lwr", 
                        "p_upr", 
                        "p_lwr_prime", 
                        "p_upr_prime",
                        "log10_p_lwr",
                        "log10_p_upr",
                        "log10_p_lwr_prime",
                        "log10_p_upr_prime"))


# Plots of posterior samples with observed value shown

post <- as.data.frame(extract(fit))

plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  geom_vline(xintercept = p, color = "red") + # MLE for p
  geom_hline(yintercept = q, color = "blue") + # MLE for q
  geom_vline(xintercept = mean(post$p_lwr), color = "red", lty = 2) + # posterior mean for p_lwr
  geom_hline(yintercept = mean(post$p_upr), color = "blue", lty = 2) + # posterior mean for p_upr
  ggtitle(expression(p[lwr]*" vs. "*p[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime, y = p_upr_prime)) +
  geom_point() +
  xlab(expression(p[lwr]*"'")) +
  ylab(expression(p[upr]*"'")) +
  geom_vline(xintercept = p_prime, color = "red") + 
  geom_hline(yintercept = q_prime, color = "blue") +
  geom_vline(xintercept = mean(post$p_lwr_prime), color = "red", lty = 2) + # posterior mean for p_lwr_prime
  geom_hline(yintercept = mean(post$p_upr_prime), color = "blue", lty = 2) + # posterior mean for p_upr_prime
  ggtitle(expression(p[lwr]*"'"*" vs. "*p[upr]*"'"))

plot3 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  geom_vline(xintercept = log10(p), color = "red") +
  geom_hline(yintercept = log10(q), color = "blue") +
  geom_vline(xintercept = mean(post$log10_p_lwr), color = "red", lty = 2) + # posterior mean for log10(p_lwr)
  geom_hline(yintercept = mean(post$log10_p_upr), color = "blue", lty = 2) + # posterior mean for log10(p_upr)
  ggtitle(expression(log[10](p[lwr])*" vs. "*log[10](p[upr])))

plot4 <- ggplot(post, aes(x = log10_p_lwr_prime, y = log10_p_upr_prime)) +
  geom_point() +
  xlab(expression(log[10](p[lwr]))) +
  ylab(expression(log[10](p[upr]))) +
  geom_vline(xintercept = log10(p_prime), color = "red") +
  geom_hline(yintercept = log10(q_prime), color = "blue") +
  geom_vline(xintercept = mean(post$log10_p_lwr_prime), color = "red", lty = 2) + # posterior mean for log10(p_lwr_prime)
  geom_hline(yintercept = mean(post$log10_p_upr_prime), color = "blue", lty = 2) + # posterior mean for log10(p_upr_prime)
  ggtitle(expression(log[10](p[lwr]*"'")*" vs. "*log[10](p[upr]*"'")))


print(plot1)
print(plot2)
print(plot3)
print(plot4)

grid.arrange(plot1, plot2, plot3, plot4, ncol=1)




plot1 <- ggplot(post, aes(x = p_lwr, y = p_upr)) +
  geom_point() +
  xlab(expression(y[lwr])) +
  ylab(expression(y[upr])) +
  geom_vline(xintercept = N * p, color = "red") +
  geom_hline(yintercept = M * q, color = "blue") +
  geom_vline(xintercept = mean(N * post$p_lwr), color = "red", lty = 2) + # posterior mean for y_lwr
  geom_hline(yintercept = mean(M * post$p_upr), color = "blue", lty = 2) + # posterior mean for y_upr
  ggtitle(expression(y[lwr]*" vs. "*y[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime, y = p_upr_prime)) +
  geom_point() +
  xlab(expression(y[lwr]*"'")) +
  ylab(expression(y[upr]*"'")) +
  geom_vline(xintercept = N * p_prime, color = "red") +
  geom_hline(yintercept = C * q_prime, color = "blue") +
  geom_vline(xintercept = mean(N * post$p_lwr_prime), color = "red", lty = 2) + # posterior mean for y_lwr_prime
  geom_hline(yintercept = mean(C * post$p_upr_prime), color = "blue", lty = 2) + # posterior mean for y_upr_prime
  ggtitle(expression(y[lwr]*"'"*" vs. "*y[upr]*"'"))

plot3 <- ggplot(post, aes(x = log10_p_lwr, y = log10_p_upr)) +
  geom_point() +
  xlab(expression(log[10](y[lwr]))) +
  ylab(expression(log[10](y[upr]))) +
  geom_vline(xintercept = N * log10(p), color = "red") +
  geom_hline(yintercept = M * log10(q), color = "blue") +
  geom_vline(xintercept = mean(N * post$log10_p_lwr), color = "red", lty = 2) + # posterior mean for log10(y_lwr)
  geom_hline(yintercept = mean(M * post$log10_p_upr), color = "blue", lty = 2) + # posterior mean for log10(y_lwr)
  ggtitle(expression(log[10](y[lwr])*" vs. "*log[10](y[upr])))

plot4 <- ggplot(post, aes(x = log10_p_lwr_prime, y = log10_p_upr_prime)) +
  geom_point() +
  xlab(expression(log[10](y[lwr]))) +
  ylab(expression(log[10](y[upr]))) +
  geom_vline(xintercept = N * log10(p_prime), color = "red") +
  geom_hline(yintercept = C * log10(q_prime), color = "blue") +
  geom_vline(xintercept = mean(N * post$log10_p_lwr_prime), color = "red", lty = 2) + # posterior mean for y_lwr_prime
  geom_hline(yintercept = mean(M * post$log10_p_upr_prime), color = "blue", lty = 2) + # posterior mean for y_lwr_prime
  ggtitle(expression(log[10](y[lwr]*"'")*" vs. "*log[10](y[upr]*"'")))


print(plot1)
print(plot2)
print(plot3)
print(plot4)

grid.arrange(plot1, plot2, plot3, plot4, ncol=1)


p1 <- ggplot(post, aes(x = p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = p, color = "red") +
  geom_vline(xintercept = mean(post$p_lwr), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]), title = expression(p[lwr]))

p2 <- ggplot(post, aes(x = p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = q, color = "blue") +
  geom_vline(xintercept = mean(post$p_upr), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]), title = expression(p[upr]))

p3 <- ggplot(post, aes(x = p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = p_prime, color = "red") +
  geom_vline(xintercept = mean(post$p_lwr_prime), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]*"'"), title = expression(p[lwr]*"'"))

p4 <- ggplot(post, aes(x = p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = q_prime, color = "blue") +
  geom_vline(xintercept = mean(post$p_upr_prime), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]*"'"), title = expression(p[upr]*"'"))

p5 <- ggplot(post, aes(x = log10_p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = log10_p, color = "red") +
  geom_vline(xintercept = mean(post$log10_p_lwr), color = "red", lty = 2) +
  labs(x = expression(log[10](p[lwr])), title = expression(log[10](p[lwr])))

p6 <- ggplot(post, aes(x = log10_p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = log10_q, color = "blue") +
  geom_vline(xintercept = mean(post$log10_p_upr), color = "blue", lty = 2) +
  labs(x = expression(log[10](p[upr])), title = expression(log[10](p[upr])))

p7 <- ggplot(post, aes(x = log10_p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = log10_p_prime, color = "red") +
  geom_vline(xintercept = mean(post$log10_p_lwr_prime), color = "red", lty = 2) +
  labs(x = expression(log[10](p[lwr]*"'")), title = expression(log[10](p[lwr]*"'")))

p8 <- ggplot(post, aes(x = log10_p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = log10_q_prime, color = "blue") +
  geom_vline(xintercept = mean(post$log10_p_upr_prime), color = "blue", lty = 2) +
  labs(x = expression(log[10](p[upr]*"'")), title = expression(log[10](p[upr]*"'")))


# Arrange plots in a 2x2 grid using facet_wrap
combined_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)


p1 <- ggplot(N * post, aes(x = p_lwr)) +
  geom_histogram() +
  geom_vline(xintercept = N * p, color = "red") +
  geom_vline(xintercept = mean(N * post$p_lwr), color = "red", lty = 2) +
  labs(x =  expression(y[lwr]), title = expression(y[lwr]))

p2 <- ggplot(M * post, aes(x = p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = M * q, color = "blue") +
  geom_vline(xintercept = mean(M * post$p_upr), color = "blue", lty = 2) +
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
  geom_vline(xintercept = 10^(N * log10_p), color = "red") +
  labs(x = expression(log[10](y[lwr])), title = expression(log[10](y[lwr])))

p6 <- ggplot(post, aes(x = log10_p_upr)) +
  geom_histogram() +
  geom_vline(xintercept = 10^(M * log10_q), color = "blue") +
  labs(x = expression(log[10](y[upr])), title = expression(log[10](y[upr])))

p7 <- ggplot(post, aes(x = log10_p_lwr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = 10^(N * log10_p_prime), color = "red") +
  labs(x = expression(log[10](y[lwr]*"'")), title = expression(log[10](y[lwr]*"'")))

p8 <- ggplot(post, aes(x = log10_p_upr_prime)) +
  geom_histogram() +
  geom_vline(xintercept = 10^(C * log10_p), color = "blue") +
  labs(x = expression(log[10](y[upr]*"'")), title = expression(log[10](y[upr]*"'")))


# Arrange plots in a 2x2 grid using facet_wrap
combined_plots <- list(p1, p2, p3, p4, p5, p6, p7, p8)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime", "log10_p_lwr", "log10_p_upr", "log10_p_lwr_prime", "log10_p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)


