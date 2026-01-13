# ==================================================s
# Monte Carlo simulation of DNA barcode gap metrics
# ==================================================

library(ggplot2)

set.seed(1)

post_mean <- function(Y, n) (Y + 1) / (n + 2)
post_var  <- function(Y, n) ((Y + 1) * (n - Y + 1)) / ((n + 2)^2 * (n + 3))

mle_mean <- function(Y, n) Y / n
mle_var  <- function(Y, n) {
  p <- Y / n
  p * (1 - p) / n
}


n_vals <- c(1, 2, 5, 10, 20, 50, 100)
theta_true_vals <- c(0, 0.1, 0.5, 0.9, 1)
N <- 10000


sim_df <- do.call(rbind, lapply(theta_true_vals, function(theta_true) {
  do.call(rbind, lapply(n_vals, function(n) {
    Y <- rbinom(N, size = n, prob = theta_true)
    data.frame(
      n = factor(n),
      theta_true = theta_true,
      Y = Y,
      bayes_mean = post_mean(Y, n),
      bayes_var  = post_var(Y, n),
      mle_mean   = mle_mean(Y, n),
      mle_var    = mle_var(Y, n)
    )
  }))
}))


df_long <- rbind(
  transform(sim_df, mean = bayes_mean, var = bayes_var, estimator = "Bayesian estimate"),
  transform(sim_df, mean = mle_mean,   var = mle_var,   estimator = "Frequentist MLE")
)
df_long$estimator <- factor(df_long$estimator,
                            levels = c("Bayesian estimate", "Frequentist MLE"))


exact_df <- do.call(rbind, lapply(n_vals, function(n) {
  Y <- 0:n
  data.frame(
    n = factor(n),
    Y = Y,
    mean_bayes = post_mean(Y, n),
    var_bayes  = post_var(Y, n),
    mean_mle   = mle_mean(Y, n),
    var_mle    = mle_var(Y, n)
  )
}))

exact_long <- rbind(
  transform(exact_df, mean = mean_bayes, var = var_bayes, estimator = "Bayesian estimate"),
  transform(exact_df, mean = mean_mle,   var = var_mle,   estimator = "Frequentist MLE")
)
exact_long$estimator <- factor(exact_long$estimator,
                               levels = c("Bayesian estimate", "Frequentist MLE"))


exact_long <- do.call(rbind, lapply(theta_true_vals, function(tt) {
  tmp <- exact_long
  tmp$theta_true <- tt
  tmp
}))


p_overlay_all <- ggplot() +
  # MC points
  geom_point(
    data = df_long,
    aes(x = mean, y = var, colour = n),
    alpha = 0.10, size = 0.9
  ) +
  # analytic curves
  geom_path(
    data = exact_long,
    aes(x = mean, y = var, colour = n, linetype = estimator,
        group = interaction(n, estimator)),
    alpha = 0.95
  ) +
  facet_wrap(~ theta_true, nrow = 1, labeller = label_bquote(theta[true] == .(theta_true))) +
  labs(
    x = "Expected Value",
    y = "Variance",
    colour = expression(paste("Sample size (", italic(n), ")")),
    linetype = "Estimator",
    title = "Monte Carlo and exact mean-variance relationship for the DNA barcode gap metrics"
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_linetype_manual(values = c("Bayesian estimate" = "dashed",
                                   "Frequentist MLE" = "solid"))

print(p_overlay_all)

# p_overlay_all + scale_y_continuous(trans = "log10")
