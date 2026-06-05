lines <- readLines("/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")
idx <- grep("class\\(result\\) <- \"because\"", lines)
idx <- idx[idx > 3100 & idx < 3250][1]

patch <- c(
  "    # Compute summary statistics",
  "    sum_stats <- if (!is.null(mcmc_samples)) summary(mcmc_samples) else NULL",
  "    if (!is.null(mcmc_samples) && n.chains > 1) {",
  "      tryCatch({",
  "        n_ch <- length(mcmc_samples)",
  "        first_chain <- as.matrix(mcmc_samples[[1]])",
  "        pnames <- colnames(first_chain)",
  "        n_params <- length(pnames)",
  "        rhat_vals <- numeric(n_params)",
  "        n_iter <- nrow(first_chain)",
  "        for (p in 1:n_params) {",
  "          chain_means <- numeric(n_ch)",
  "          chain_vars <- numeric(n_ch)",
  "          for (c in 1:n_ch) {",
  "            vals <- as.matrix(mcmc_samples[[c]])[, p]",
  "            chain_means[c] <- mean(vals)",
  "            chain_vars[c] <- var(vals)",
  "          }",
  "          grand_mean <- mean(chain_means)",
  "          B <- n_iter * var(chain_means)",
  "          W <- mean(chain_vars)",
  "          if (W > 0) {",
  "            var_plus <- ((n_iter - 1) / n_iter) * W + (1 / n_iter) * B",
  "            rhat_vals[p] <- sqrt(var_plus / W)",
  "          } else {",
  "            rhat_vals[p] <- 1.0",
  "          }",
  "        }",
  "        sum_stats$statistics <- cbind(sum_stats$statistics, Rhat = rhat_vals)",
  "      }, error = function(e) {})",
  "    }",
  "    result$summary <- sum_stats",
  ""
)

lines <- c(lines[1:(idx-1)], patch, lines[idx:length(lines)])
writeLines(lines, "/Users/achazhardenberg/Library/CloudStorage/Dropbox/Repos/because/R/because.R")
