##### Set working directory #####

setwd("/Users/jarrettphillips/desktop/Bayesian DNA Barcode Gap")

library(ape)
library(spider)
library(boot)
library(bayesboot)

library(shinystan)

data(anoteropsis)
anoDist <- dist.dna(anoteropsis)
anoSpp <- sapply(strsplit(dimnames(anoteropsis)[[1]], split="_"),
                  function(x) paste(x[1], x[2], sep="_"))

intra <- maxInDist(anoDist, anoSpp)
inter <- nonConDist(anoDist, anoSpp)

res <- na.omit(cbind(intra, inter))
N <- nrow(res)

intra <- res[, "intra"]
inter <- res[, "inter"]

(p <- mean(intra >= min(inter)))
(q <- mean(inter <= max(intra)))

dens_intra <- density(intra)
dens_inter <- density(inter)
plot(dens_intra, col = "blue", 
     xlim = c(-0.01, 0.06), 
     ylim = c(0, 50), 
     main = "Intraspecific and interspecific genetic distance distributions")
lines(dens_inter, col = "red")
abline(v = c(min(inter), max(intra)), lty = 2, lwd = 2, col = c("red", "blue"))

legend("topright", 
       legend = c("Intraspecific", "Interspecific"), 
       col = c("blue", "red"), 
       lty = 1, 
       lwd = 2)



library(ggplot2)

# Create density data
dens_intra <- density(intra)
dens_inter <- density(inter)
data_intra <- data.frame(x = dens_intra$x, y = dens_intra$y)
data_inter <- data.frame(x = dens_inter$x, y = dens_inter$y)

library(ggplot2)

# Create density data
dens_intra <- density(intra)
dens_inter <- density(inter)
data_intra <- data.frame(x = dens_intra$x, y = dens_intra$y)
data_inter <- data.frame(x = dens_inter$x, y = dens_inter$y)

# Find x values where the blue dashed line intersects with the red curve
x_intersect_blue <- max(intra)
y_intersect_blue <- approx(data_inter$x, data_inter$y, x_intersect_blue)$y

# Find x values where the red dashed line intersects with the blue curve
x_intersect_red <- min(inter)
y_intersect_red <- approx(data_intra$x, data_intra$y, x_intersect_red)$y

# Plot using ggplot2
ggplot() +
  geom_line(data = data_intra, aes(x = x, y = y), color = "blue") +
  geom_line(data = data_inter, aes(x = x, y = y), color = "red") +
  geom_vline(xintercept = c(min(inter), max(intra)), linetype = "dashed", 
             color = c("red", "blue"), linewidth = 1) +
  geom_ribbon(data = subset(data_inter, x <= x_intersect_red), aes(x = x, ymin = 0, ymax = y),
              fill = "red", alpha = 0.3) +
  geom_ribbon(data = subset(data_intra, x >= x_intersect_blue), aes(x = x, ymin = 0, ymax = y),
              fill = "blue", alpha = 0.3) +
  xlim(-0.01, 0.06) +
  ylim(0, 50) +
  labs(title = "Intraspecific and interspecific genetic distance distributions") +
  labs(color = "Legend") +
  scale_color_manual(values = c("blue", "red"), labels = c("Intraspecific", "Interspecific")) 

### Bootstrap ###

p <- function(res, i) {
  mean(res[i, 1] >= min(res[, 2]))
}

q <- function(res, i) {
  mean(res[i, 2] <= max(res[, 1]))
}

(boot_p <- boot(res, p, 10000))
(boot_q <- boot(res, q, 10000))

plot(boot_p)
plot(boot_q)

bayesboot_p <- bayesboot(res, p, 10000)
summary(bayesboot_p)

bayesboot_q <- bayesboot(res, q, 10000)
summary(bayesboot_q)


plot(bayesboot_p)




##### Install required packages #####

# install.packages("rstan") # for HMC

##### Load required libraries #####

library(rstan) # for MCMC

options(mc.cores = parallel::detectCores()) # parallelize simulations
rstan_options(auto_write = TRUE) # only have to compile once unless code is changed

fit <- stan("barcode_gap.stan", 
            data = list(N = N, intra = intra, inter = inter),
            iter = 2000,
            control = list(adapt_delta = 0.80))
fit

traceplot(fit)


# Plot histograms of posterior samples with observed value shown

post <- extract(fit)

par(mfrow = c(1, 2))
plot(post$p_lwr, post$p_upr, xlab = "p_lwr",  ylab = "p_upr", main = paste0(expression(p[lwr], "vs.", expression(p[upr]))))
plot(post$log10_p_lwr, post$log10_p_upr, xlab = expression(log[10](p_lwr)), ylab = expression(log[10](p_upr)), main = paste0(expression(log[10](p_lwr)), " vs. " , expression(log[10](p_upr))))


par(mfrow = c(2, 2))
hist(post$p_lwr, xlab = "p_lwr", main = expression(p[lwr]))
abline(v = p, lwd = 5)
hist(post$p_upr, xlab = "p_upr",  main = expression(p[upr]))
abline(v = q, lwd = 5)
hist(post$log10_p_lwr, xlab = expression(log[10](p[lwr])), main = expression(log[10](p[lwr])))
abline(v = log10(p), lwd = 5)
hist(post$log10_p_upr, xlab = expression(log[10](p[upr])), main = expression(log[10](p[upr])))
abline(v = log10(q), lwd = 5)



