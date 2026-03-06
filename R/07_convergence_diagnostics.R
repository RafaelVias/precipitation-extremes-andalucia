# 07_convergence_diagnostics.R — MCMC convergence diagnostics for Stage 2
#
# Produces trace plots, Rhat/ESS summaries, and pairs plots for the
# 12 GP hyperparameters from the Matérn(5/2) + PC prior Stan model.
#
# Outputs saved to diagnostics/ (not shown in README).
#
# Run from project root: Rscript R/07_convergence_diagnostics.R

library(cmdstanr)
library(bayesplot)
library(ggplot2)
library(patchwork)

color_scheme_set("mix-blue-red")

cat("========================================\n")
cat("Convergence diagnostics\n")
cat("========================================\n")

# ---- 1. Load fit ----
cat("Loading Stage 2 results...\n")

s2  <- readRDS("data/stage2_matern_pc_results.rds")
fit <- s2$stan_fit

# Hyperparameter names
hyper_vars <- c(
  "mu[1]", "mu[2]", "mu[3]",
  "sigma_gp[1]", "sigma_gp[2]", "sigma_gp[3]",
  "phi_gp[1]", "phi_gp[2]", "phi_gp[3]",
  "nugget[1]", "nugget[2]", "nugget[3]"
)

# Nice labels: parameter[index] -> GEV param name
hyper_labels <- c(
  "mu[1]" = expression(mu[psi]), "mu[2]" = expression(mu[tau]), "mu[3]" = expression(mu[phi]),
  "sigma_gp[1]" = expression(sigma[psi]), "sigma_gp[2]" = expression(sigma[tau]),
  "sigma_gp[3]" = expression(sigma[phi]),
  "phi_gp[1]" = expression(rho[psi]), "phi_gp[2]" = expression(rho[tau]),
  "phi_gp[3]" = expression(rho[phi]),
  "nugget[1]" = expression(nu[psi]), "nugget[2]" = expression(nu[tau]),
  "nugget[3]" = expression(nu[phi])
)

draws_array <- fit$draws(variables = hyper_vars, format = "draws_array")

cat("  Chains:", dim(draws_array)[2], "\n")
cat("  Draws per chain:", dim(draws_array)[1], "\n")

# Print summary table
cat("\n  Hyperparameter summary:\n")
summ <- fit$summary(variables = hyper_vars)
print(summ[, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk", "ess_tail")])

diag_summ <- fit$diagnostic_summary()
cat("\n  Divergences per chain:", diag_summ$num_divergent, "\n")
cat("  Max treedepth per chain:", diag_summ$num_max_treedepth, "\n")

dir.create("diagnostics", showWarnings = FALSE)

# ---- 2. Trace plots ----
cat("\nGenerating trace plots...\n")

p_trace <- mcmc_trace(draws_array, pars = hyper_vars, n_warmup = 0,
                      facet_args = list(ncol = 3, labeller = label_parsed)) +
  labs(title = "Trace plots: GP hyperparameters",
       subtitle = "4 chains \u00d7 1000 post-warmup draws") +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, colour = "grey40"),
        strip.text = element_text(size = 10))

ggsave("diagnostics/trace_plots.png", p_trace,
       width = 14, height = 10, dpi = 200, bg = "white")
cat("  Saved diagnostics/trace_plots.png\n")

# ---- 3. Rhat and ESS ----
cat("Generating Rhat/ESS plots...\n")

# Get summary for ALL parameters (hyperparams + latent field)
all_summ <- fit$summary()

rhat_vals <- all_summ$rhat
names(rhat_vals) <- all_summ$variable

ess_ratio <- all_summ$ess_bulk / (dim(draws_array)[1] * dim(draws_array)[2])
names(ess_ratio) <- all_summ$variable

# Filter out lp__ and eta_raw (internal)
keep <- !grepl("^(lp__|eta_raw)", names(rhat_vals))
rhat_vals <- rhat_vals[keep]
ess_ratio <- ess_ratio[keep]

p_rhat <- mcmc_rhat(rhat_vals) +
  labs(title = expression(hat(R) ~ "for all parameters"),
       subtitle = "Target: < 1.01") +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

p_ess <- mcmc_neff(ess_ratio) +
  labs(title = expression(N[eff] / N ~ "ratio for all parameters"),
       subtitle = "Target: > 0.1") +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

p_rhat_ess <- p_rhat + p_ess +
  plot_annotation(
    title = "Convergence diagnostics: all model parameters",
    subtitle = sprintf("Hyperparameters + %d latent field parameters",
                       sum(grepl("^eta\\[", all_summ$variable))),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("diagnostics/rhat_ess.png", p_rhat_ess,
       width = 14, height = 7, dpi = 200, bg = "white")
cat("  Saved diagnostics/rhat_ess.png\n")

# Print worst Rhat and ESS
cat(sprintf("  Max Rhat: %.4f\n", max(rhat_vals, na.rm = TRUE)))
cat(sprintf("  Min ESS ratio: %.3f\n", min(ess_ratio, na.rm = TRUE)))

# ---- 4. Pairs plots ----
cat("Generating pairs plots...\n")

# Split into 3 groups (one per GEV parameter) for readability
gev_names <- c("psi", "tau", "phi")
groups <- list(
  c("mu[1]", "sigma_gp[1]", "phi_gp[1]", "nugget[1]"),
  c("mu[2]", "sigma_gp[2]", "phi_gp[2]", "nugget[2]"),
  c("mu[3]", "sigma_gp[3]", "phi_gp[3]", "nugget[3]")
)

# Get divergence info
np <- nuts_params(fit)

# mcmc_pairs returns a bayesplot object — save directly in the loop
for (g in seq_along(groups)) {
  fname <- sprintf("diagnostics/pairs_%s.png", gev_names[g])
  png(fname, width = 10, height = 10, units = "in", res = 200, bg = "white")
  print(mcmc_pairs(draws_array, pars = groups[[g]], np = np,
                   off_diag_args = list(size = 0.5, alpha = 0.3)))
  dev.off()
  cat("  Saved", fname, "\n")
}

# ---- 5. Autocorrelation plots ----
cat("Generating autocorrelation plots...\n")

p_acf <- mcmc_acf(draws_array, pars = hyper_vars, lags = 40) +
  labs(title = "Autocorrelation: GP hyperparameters",
       subtitle = "4 chains, 40 lags") +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, colour = "grey40"))

ggsave("diagnostics/autocorrelation.png", p_acf,
       width = 14, height = 10, dpi = 200, bg = "white")
cat("  Saved diagnostics/autocorrelation.png\n")

cat("\n========================================\n")
cat("Done. All diagnostics saved to diagnostics/\n")
cat("========================================\n")
