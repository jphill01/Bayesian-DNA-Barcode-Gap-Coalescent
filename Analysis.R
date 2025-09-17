#############################################
# R script to run analyses in Phillips et al.
# Created by: Jarrett D. Phillips
# Last updated: September 4, 2025
############################################


# install.packages("ggplot2")
# install.packages("rstan")
# install.packages("dplyr")
# install.packages("parallel")
# install.packages("doParallel")
# install.packages("reshape2")
# install.packages("rstudioapi")

library(ggplot2)
library(rstan)
library(dplyr)
library(parallel)
library(doParallel)
library(reshape2)
library(rstudioapi)

options(mc.cores = detectCores())
rstan_options(auto_write = TRUE)

setwd("/Users/jarrettphillips/desktop/Bayesian DNA Barcode Gap Analysis")


run_DNA_barcode_gap_analysis <- function(data_dir = NULL,
                                         stan_file = "DNA_barcode_gap.stan") {


  if (is.null(data_dir)) {
    message("Select base data directory (with marker folders)...")
    data_dir <- selectDirectory(caption = "Select data directory with marker folders")
  }

  if (!file.exists(stan_file)) stop("Stan model file not found: ", stan_file)

  markers <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)

  for (marker_dir in markers) {
    marker <- basename(marker_dir)
    message("Processing marker: ", marker)

    inter_file <- file.path(marker_dir, "inter.csv")
    if (!file.exists(inter_file)) {
      warning("'inter.csv' not found in: ", marker_dir)
      next
    }

    inter <- read.csv(inter_file)[, 2]
    M <- length(inter)

    intra_files <- list.files(marker_dir, pattern = "_intra\\.csv$", full.names = TRUE)
    comb_files  <- list.files(marker_dir, pattern = "_inter\\.csv$", full.names = TRUE)

    species_intra <- gsub("_intra\\.csv$", "", basename(intra_files))
    species_comb  <- gsub("_inter\\.csv$", "", basename(comb_files))

    species <- intersect(species_intra, species_comb)
    K <- length(species)

    if (K == 0) {
      warning("No matching intraspecific and interspecific species files found in ", marker_dir)
      next
    }

    intra_list <- comb_list <- list()
    N <- C <- numeric(K)

    for (i in seq_along(species)) {
      intra_path <- file.path(marker_dir, paste0(species[i], "_intra.csv"))
      comb_path  <- file.path(marker_dir, paste0(species[i], "_inter.csv"))
      intra_list[[i]] <- read.csv(intra_path)
      comb_list[[i]]  <- read.csv(comb_path)
      N[i] <- nrow(intra_list[[i]])
      C[i] <- nrow(comb_list[[i]])
    }

    intra_combined <- unlist(lapply(intra_list, function(df) df$x))
    comb_combined  <- unlist(lapply(comb_list, function(df) df$x))

    stan_data <- list(K = K,
                      M = M,
                      N = N,
                      intra = intra_combined,
                      inter = inter,
                      C = C,
                      comb = comb_combined)

    marker_out <- file.path(marker_dir, "Results")
    dir.create(marker_out, recursive = TRUE, showWarnings = FALSE)

    message("Running Stan model for marker: ", marker)
    fit <- stan(file = stan_file, data = stan_data,
                chains = 4, iter = 2000, seed = 673227,
                control = list(adapt_delta = 0.85, max_treedepth = 10))

    saveRDS(fit, file.path(marker_out, "stan_fit.rds"))
    post <- as.data.frame(extract(fit))

    message("Saving results per species...")

    param_labels <- list(
      p_lwr        = expression(p[lwr]),
      p_upr        = expression(p[upr]),
      p_lwr_prime  = expression(p[lwr]^"'"),
      p_upr_prime  = expression(p[upr]^"'")
    )

    for (i in seq_along(species)) {
      sp <- species[i]
      sp_dir <- file.path(marker_out, sp)
      dir.create(sp_dir, showWarnings = FALSE)

      intra_x <- intra_list[[i]]$x
      comb_x  <- comb_list[[i]]$x

      p        <- mean(intra_x >= min(inter))
      q        <- mean(inter <= max(intra_x))
      p_prime  <- mean(intra_x >= min(comb_x))
      q_prime  <- mean(comb_x <= max(intra_x))

      se_p        <- sqrt(p * (1 - p) / length(intra_x))
      se_q        <- sqrt(q * (1 - q) / length(inter))
      se_p_prime  <- sqrt(p_prime * (1 - p_prime) / length(intra_x))
      se_q_prime  <- sqrt(q_prime * (1 - q_prime) / length(comb_x))

      ci_p        <- p + c(-1, 1) * qnorm(0.975) * se_p
      ci_q        <- q + c(-1, 1) * qnorm(0.975) * se_q
      ci_p_prime  <- p_prime + c(-1, 1) * qnorm(0.975) * se_p_prime
      ci_q_prime  <- q_prime + c(-1, 1) * qnorm(0.975) * se_q_prime

      p_lwr_col <- paste0("p_lwr.", i)
      p_upr_col <- paste0("p_upr.", i)
      p_lwr_prime_col <- paste0("p_lwr_prime.", i)
      p_upr_prime_col <- paste0("p_upr_prime.", i)

      p_lwr_mean <- mean(post[[p_lwr_col]])
      p_upr_mean <- mean(post[[p_upr_col]])
      p_lwr_prime_mean <- mean(post[[p_lwr_prime_col]])
      p_upr_prime_mean <- mean(post[[p_upr_prime_col]])

      est_df <- data.frame(
        species = sp,
        p = p, se_p = se_p, ci_p_lower = ci_p[1], ci_p_upper = ci_p[2],
        q = q, se_q = se_q, ci_q_lower = ci_q[1], ci_q_upper = ci_q[2],
        p_prime = p_prime, se_p_prime = se_p_prime, ci_p_prime_lower = ci_p_prime[1], ci_p_prime_upper = ci_p_prime[2],
        q_prime = q_prime, se_q_prime = se_q_prime, ci_q_prime_lower = ci_q_prime[1], ci_q_prime_upper = ci_q_prime[2]
      )

      write.csv(est_df, file.path(sp_dir, "estimates.csv"), row.names = FALSE)

      species_post <- select(post, matches(paste0("\\.", i, "$")))
      saveRDS(species_post, file.path(sp_dir, "posterior_samples.rds"))

      params <- c("p_lwr", "p_upr", "p_lwr_prime", "p_upr_prime")
      values <- c(p, q, p_prime, q_prime)
      colors <- c("red", "blue", "red", "blue")

      plot1 <- ggplot(post, aes_string(x = p_lwr_col, y = p_upr_col)) +
        geom_point() +
        xlab(param_labels[["p_lwr"]]) +
        ylab(param_labels[["p_upr"]]) +
        geom_vline(xintercept = p, color = "red") +
        geom_hline(yintercept = q, color = "blue") +
        geom_vline(xintercept = p_lwr_mean, color = "red", linetype = 2) +
        geom_hline(yintercept = p_upr_mean, color = "blue", linetype = 2) +
        ggtitle(bquote(italic(.(paste("A.", sp)))))

      plot2 <- ggplot(post, aes_string(x = p_lwr_prime_col, y = p_upr_prime_col)) +
        geom_point() +
        xlab(param_labels[["p_lwr_prime"]]) +
        ylab(param_labels[["p_upr_prime"]]) +
        geom_vline(xintercept = p_prime, color = "red") +
        geom_hline(yintercept = q_prime, color = "blue") +
        geom_vline(xintercept = p_lwr_prime_mean, color = "red", linetype = 2) +
        geom_hline(yintercept = p_upr_prime_mean, color = "blue", linetype = 2) +
        ggtitle(bquote(italic(.(paste("A.", sp)))))

      ggsave(file.path(sp_dir, paste0("posterior_lwr_upr_", sp, ".png")),
             plot1, width = 6, height = 4)

      ggsave(file.path(sp_dir, paste0("posterior_lwr_prime_upr_prime_", sp, ".png")),
             plot2, width = 6, height = 4)

      for (j in seq_along(params)) {
        var <- paste0(params[j], ".", i)

        trace <- traceplot(fit, pars = params[j])

        ggsave(file.path(marker_out, paste0("traceplot_", params[j], ".png")),
               plot = trace, width = 14, height = 10)

        density_plot <- ggplot(post, aes_string(x = var)) +
          geom_density() +
          geom_vline(xintercept = values[j], color = colors[j]) +
          geom_vline(xintercept = mean(post[[var]]), color = colors[j], linetype = 2) +
          labs(
            title = bquote(italic("A. ") ~ italic(.(sp)) ~ .(param_labels[[params[j]]])),
            x = param_labels[[params[j]]]
          )

        ggsave(filename = file.path(sp_dir, paste0("density_", params[j], ".png")),
               plot = density_plot, width = 6, height = 4)
      }

      ecdf_intra <- data.frame(x = sort(intra_x), y = 1 - ecdf(intra_x)(sort(intra_x)) + mean(intra_x == min(inter)))
      ecdf_inter <- data.frame(x = sort(inter), y = ecdf(inter)(sort(inter)))
      ecdf_comb  <- data.frame(x = sort(comb_x), y = ecdf(comb_x)(sort(comb_x)))

      p1 <- ggplot(ecdf_intra, aes(x = x, y = y)) +
        geom_step() + geom_vline(xintercept = min(inter), linetype = "dashed") +
        labs(title = bquote(italic("A. ") ~ italic(.(sp))),
             x = expression(d[ij]),
             y = expression(1 - hat(F)(d[ij]) + P(d[ij] == a)))

      p2 <- ggplot(ecdf_inter, aes(x = x, y = y)) +
        geom_step() + geom_vline(xintercept = max(intra_x), linetype = "dashed") +
        labs(title = bquote(italic("A. ") ~ italic(.(sp))),
             x = expression(d[XY]),
             y = expression(hat(F)(d[XY])))

      p3 <- ggplot(ecdf_intra, aes(x = x, y = y)) +
        geom_step() + geom_vline(xintercept = min(comb_x), linetype = "dashed") +
        labs(title = bquote(italic("A. ") ~ italic(.(sp))),
             x = expression(d[ij]),
             y = expression(1 - hat(F)(d[ij]) + P(d[ij] == a^"'")))

      p4 <- ggplot(ecdf_comb, aes(x = x, y = y)) +
        geom_step() + geom_vline(xintercept = max(intra_x), linetype = "dashed") +
        labs(title = bquote(italic("A. ") ~ italic(.(sp))),
             x = expression(d[XY]^"'"),
             y = expression(hat(F)(d[XY]^"'")))

      ggsave(file.path(sp_dir, "ecdf_intra_inter.png"), p1, width = 6, height = 4)
      ggsave(file.path(sp_dir, "ecdf_inter_intra.png"), p2, width = 6, height = 4)
      ggsave(file.path(sp_dir, "ecdf_intra_comb.png"), p3, width = 6, height = 4)
      ggsave(file.path(sp_dir, "ecdf_comb_intra.png"), p4, width = 6, height = 4)
    }

    message("Completed marker: ", marker)
  }

  message("All analyses complete. Results saved inside each marker folder.")
}



# Run

run_DNA_barcode_gap_analysis()

