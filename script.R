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


##### Setup #####

K <- 2

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


##### MLEs #####

(p_1 <- mean(intra1$x >= min(inter)))
(q_1 <- mean(inter <= max(intra1$x)))
(p_prime_1 <- mean(intra1$x >= min(comb1$x)))
(q_prime_1 <- mean(comb1$x <= max(intra1$x)))

(p_2 <- mean(intra2$x >= min(inter)))
(q_2 <- mean(inter <= max(intra2$x)))
(p_prime_2 <- mean(intra2$x >= min(comb2$x)))
(q_prime_2 <- mean(comb2$x <= max(intra2$x)))


##### SEs #####

(SE_p_1 <- sqrt(p_1 * (1 - p_1) / length(intra1)))
(SE_q_1 <- sqrt(q_1 * (1 - q_1) / length(inter)))
(SE_p_prime_1 <- sqrt(p_prime_1 * (1 - p_prime_1) / length(intra1)))
(SE_q_prime_1 <- sqrt(q_prime_1 * (1 - q_prime_1) / length(comb1)))


(SE_p_2 <- sqrt(p_2 * (1 - p_2) / length(intra2)))
(SE_q_2 <- sqrt(q_2 * (1 - q_2) / length(inter)))
(SE_p_prime_2 <- sqrt(p_prime_2 * (1 - p_prime_2) / length(intra2)))
(SE_q_prime_2 <- sqrt(q_prime_2 * (1 - q_prime_2) / length(comb2)))


##### 95% CIs #####

p_1 + c(-1, 1) * qnorm(0.975) * SE_p_1
q_1 + c(-1, 1) * qnorm(0.975) * SE_q_1
p_prime_1 + c(-1, 1) * qnorm(0.975) * SE_p_prime_1
q_prime_1 + c(-1, 1) * qnorm(0.975) * SE_q_prime_1


p_2 + c(-1, 1) * qnorm(0.975) * SE_p_2
q_2 + c(-1, 1) * qnorm(0.975) * SE_q_2
p_prime_2 + c(-1, 1) * qnorm(0.975) * SE_p_prime_2
q_prime_2 + c(-1, 1) * qnorm(0.975) * SE_q_prime_2




# ECDFs

ecdf_intra1 <- ecdf(intra1$x)
ecdf_intra2 <- ecdf(intra2$x)
ecdf_inter <- ecdf(inter)
ecdf_comb1 <- ecdf(comb1$x)
ecdf_comb2 <- ecdf(comb2$x)

a1 <- min(inter)
b1 <- max(intra1$x)
a_prime1 <- min(comb1$x)

a2 <- min(inter)
b2 <- max(intra2$x)
a_prime2 <- min(comb2$x)

p1_ecdf <- 1 - ecdf_intra1(a1) + mean(intra1$x == min(inter))
q1_ecdf <- ecdf_inter(b1)
p1_prime_ecdf <- 1 - ecdf_intra1(a_prime1) + mean(intra1$x == min(comb1$x))
q1_prime_ecdf <- ecdf_comb1(b1)

p2_ecdf <- 1 - ecdf_intra2(a2) + mean(intra2$x == min(inter))
q2_ecdf <- ecdf_inter(b2)
p2_prime_ecdf <- 1 - ecdf_intra2(a_prime2) + mean(intra2$x == min(comb2$x))
q2_prime_ecdf <- ecdf_comb2(b2)

plot(ecdf_intra1)
plot(ecdf_intra2)
plot(ecdf_inter)
plot(ecdf_comb1)
plot(ecdf_comb2)

##### Observed counts #####

N[1] * p_1
M * q_1
N[1] * p_prime_1
C[1] * q_prime_1

N[2] * p_2
M * q_2
N[2] * p_prime_2
C[2] * q_prime_2


##### Posterior Estimates #####

fit <- stan("DNA_barcode_gap.stan", 
            data = list(K = K, M = M, N = N, intra = intra, inter = inter, C = C, comb = comb), 
            chains = 4,
            iter = 2000,
            seed = 0673227,
            algorithm = "NUTS",
            control = list(adapt_delta = 0.80,
                           max_treedepth = 10))

print(fit, digits_summary = 3)

traceplot(fit, pars = c("p_lwr", 
                        "p_upr", 
                        "p_lwr_prime", 
                        "p_upr_prime"))

##### Plots #####

post <- as.data.frame(extract(fit))

summary(post)

plot1 <- ggplot(post, aes(x = p_lwr.1, y = p_upr.1)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  geom_vline(xintercept = p_1, color = "red") + # MLE for p
  geom_hline(yintercept = q_1, color = "blue") + # MLE for q
  geom_vline(xintercept = mean(as.numeric(post$p_lwr.1)), color = "red", lty = 2) + # posterior mean for p_lwr
  geom_hline(yintercept = mean(as.numeric(post$p_upr.1)), color = "blue", lty = 2) + # posterior mean for p_upr
  ggtitle(expression(italic("A. bipustulatus") ~ p[lwr]*" vs. "*p[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime.1, y = p_upr_prime.1)) +
  geom_point() +
  xlab(expression(p[lwr]*"'")) +
  ylab(expression(p[upr]*"'")) +
  geom_vline(xintercept = p_prime_1, color = "red") + 
  geom_hline(yintercept = q_prime_1, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_prime.1)), color = "red", lty = 2) + # posterior mean for p_lwr_prime
  geom_hline(yintercept = mean(as.numeric(post$p_upr_prime.1)), color = "blue", lty = 2) + # posterior mean for p_upr_prime
  ggtitle(expression(italic("A. bipustulatus") ~ p[lwr]*"'"*" vs. "*p[upr]*"'"))

print(plot1)
print(plot2)

grid.arrange(plot1, plot2, ncol = 1)


plot1 <- ggplot(post, aes(x = p_lwr.2, y = p_upr.2)) +
  geom_point() +
  xlab(expression(p[lwr])) +
  ylab(expression(p[upr])) +
  geom_vline(xintercept = p_2, color = "red") + # MLE for p
  geom_hline(yintercept = q_2, color = "blue") + # MLE for q
  geom_vline(xintercept = mean(as.numeric(post$p_lwr.2)), color = "red", lty = 2) + # posterior mean for p_lwr
  geom_hline(yintercept = mean(as.numeric(post$p_upr.2)), color = "blue", lty = 2) + # posterior mean for p_upr
  ggtitle(expression(italic("A. nevadensis") ~ p[lwr]*" vs. "*p[upr]))

plot2 <- ggplot(post, aes(x = p_lwr_prime.2, y = p_upr_prime.2)) +
  geom_point() +
  xlab(expression(p[lwr]*"'")) +
  ylab(expression(p[upr]*"'")) +
  geom_vline(xintercept = p_prime_2, color = "red") + 
  geom_hline(yintercept = q_prime_2, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_prime.2)), color = "red", lty = 2) + # posterior mean for p_lwr_prime
  geom_hline(yintercept = mean(as.numeric(post$p_upr_prime.2)), color = "blue", lty = 2) + # posterior mean for p_upr_prime
  ggtitle(expression(italic("A. nevadensis") ~ p[lwr]*"'"*" vs. "*p[upr]*"'"))

print(plot1)
print(plot2)

grid.arrange(plot1, plot2, ncol = 1)


p1 <- ggplot(post, aes(x = p_lwr.1)) +
  geom_density() +
  geom_vline(xintercept = p_1, color = "red") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_1)), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]), title = expression(p[lwr]))

p2 <- ggplot(post, aes(x = p_upr.1)) +
  geom_density() +
  geom_vline(xintercept = q_1, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_upr_1)), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]), title = expression(p[upr]))

p3 <- ggplot(post, aes(x = p_lwr_prime.1)) +
  geom_density() +
  geom_vline(xintercept = p_prime_1, color = "red") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_prime.1)), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]*"'"), title = expression(p[lwr]*"'"))

p4 <- ggplot(post, aes(x = p_upr_prime.1)) +
  geom_density() +
  geom_vline(xintercept = q_prime_1, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_upr_prime.1)), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]*"'"), title = expression(p[upr]*"'"))

combined_plots <- list(p1, p2, p3, p4)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)



p1 <- ggplot(post, aes(x = p_lwr.2)) +
  geom_density() +
  geom_vline(xintercept = p_2, color = "red") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_2)), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]), title = expression(p[lwr]))

p2 <- ggplot(post, aes(x = p_upr.2)) +
  geom_density() +
  geom_vline(xintercept = q_2, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_upr_2)), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]), title = expression(p[upr]))

p3 <- ggplot(post, aes(x = p_lwr_prime.2)) +
  geom_density() +
  geom_vline(xintercept = p_prime_2, color = "red") +
  geom_vline(xintercept = mean(as.numeric(post$p_lwr_prime.2)), color = "red", lty = 2) +
  labs(x =  expression(p[lwr]*"'"), title = expression(p[lwr]*"'"))

p4 <- ggplot(post, aes(x = p_upr_prime.2)) +
  geom_density() +
  geom_vline(xintercept = q_prime_2, color = "blue") +
  geom_vline(xintercept = mean(as.numeric(post$p_upr_prime.2)), color = "blue", lty = 2) +
  labs(x =  expression(p[upr]*"'"), title = expression(p[upr]*"'"))

combined_plots <- list(p1, p2, p3, p4)
names(combined_plots) <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime")

grid.arrange(grobs = combined_plots, ncol = 2)



