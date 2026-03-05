# =============================================================================
# 03b_mcmc_rerun.R
# Re-run MCMC with warm start and longer chain for better hyperparameter mixing
# =============================================================================

library(dplyr)
library(INLA)
library(Matrix)
library(spam)
library(fields)
library(sf)
library(ggplot2)
library(coda)

source("vendor/max_and_smooth/stage1_functions.R")
source("vendor/max_and_smooth/stage2_functions.R")

# ---- Load Stage 1 results and previous MCMC ---
stage1 <- readRDS("data/stage1_results.rds")
prev   <- readRDS("data/stage2_mcmc_results.rds")

mles.covmats <- stage1$mles.covmats
loc          <- stage1$loc
ns           <- nrow(loc)

# ---- Rebuild INLA mesh and matrices (same settings as before) ---
mesh.domain <- inla.mesh.2d(
  loc = loc,
  max.edge = c(0.3, 1.0),
  offset = c(0.2, 0.5),
  cutoff = 0.15
)

A <- inla.spde.make.A(mesh = mesh.domain, loc = loc)
fem.mesh <- inla.mesh.fem(mesh.domain, order = 2)
inla.mats <- list(
  c.mat = fem.mesh$c0,
  g1.mat = fem.mesh$g1,
  g2.mat = fem.mesh$g2,
  A = A
)

nmesh <- mesh.domain$n
cat("Mesh:", nmesh, "nodes\n")

# ---- Extract warm-start values from previous run ---
# Use posterior means / last values from previous MCMC
n_prev <- length(prev$beta_psi)

beta_psi_init  <- mean(prev$beta_psi)
beta_tau_init  <- mean(prev$beta_tau)
beta_phi_init  <- mean(prev$beta_phi)
sigmaSq_psi_init <- mean(prev$sigma_psi)^2
sigmaSq_tau_init <- mean(prev$sigma_tau)^2
sigmaSq_phi_init <- mean(prev$sigma_phi)^2
s2_psi_init <- mean(prev$s_psi)^2
s2_tau_init <- mean(prev$s_tau)^2
rho_psi_init <- mean(prev$rho_psi)
rho_tau_init <- mean(prev$rho_tau)
w_psi_init <- prev$w_psi.posmean
w_tau_init <- prev$w_tau.posmean

cat("Warm-start values:\n")
cat(sprintf("  beta_psi=%.3f, beta_tau=%.3f, beta_phi=%.4f\n",
            beta_psi_init, beta_tau_init, beta_phi_init))
cat(sprintf("  sigma_psi=%.3f, sigma_tau=%.3f, sigma_phi=%.4f\n",
            sqrt(sigmaSq_psi_init), sqrt(sigmaSq_tau_init), sqrt(sigmaSq_phi_init)))
cat(sprintf("  s_psi=%.3f, s_tau=%.3f\n", sqrt(s2_psi_init), sqrt(s2_tau_init)))
cat(sprintf("  rho_psi=%.3f, rho_tau=%.3f\n", rho_psi_init, rho_tau_init))

# ---- Run longer MCMC ---
cat("\nStarting MCMC: 50,000 iterations, burn=5000, thin=5...\n")
cat("Expected time: ~50-60 minutes\n\n")

fit.mcmc <- mcmc.maxNsmooth(
  mles.covmats = mles.covmats,
  loc = loc,
  inla.mats = inla.mats,
  alpha = 2,
  # warm start
  beta_psi.init = beta_psi_init,
  beta_tau.init = beta_tau_init,
  beta_phi.init = beta_phi_init,
  w_psi.init = w_psi_init,
  w_tau.init = w_tau_init,
  sigmaSq_psi.init = sigmaSq_psi_init,
  sigmaSq_tau.init = sigmaSq_tau_init,
  sigmaSq_phi.init = sigmaSq_phi_init,
  s2_psi.init = s2_psi_init,
  s2_tau.init = s2_tau_init,
  rho_psi.init = rho_psi_init,
  rho_tau.init = rho_tau_init,
  # priors (same as before)
  sd_beta_psi = 1e2,
  sd_beta_tau = 1e2,
  sd_beta_phi = 1e2,
  lambda_sigma_psi = 0.1,
  lambda_sigma_tau = 0.1,
  lambda_sigma_phi = 0.1,
  lambda_s_psi = 0.1,
  lambda_s_tau = 0.1,
  lambda_rho_psi = 0.1,
  lambda_rho_tau = 0.1,
  # longer chain
  iters = 50000,
  burn = 5000,
  thin = 5
)

cat("\nMCMC completed in", round(fit.mcmc$minutes, 1), "minutes\n")

# ---- Diagnostics ---
cat("\n=== Effective sample sizes ===\n")
hyper_names <- c("beta_psi", "beta_tau", "beta_phi",
                 "sigma_psi", "sigma_tau", "sigma_phi",
                 "s_psi", "s_tau", "rho_psi", "rho_tau")
for (h in hyper_names) {
  chain <- fit.mcmc[[h]]
  ess <- effectiveSize(chain)
  cat(sprintf("  %-15s ESS = %6.0f / %d   mean=%.4f  95%%CI=(%.4f, %.4f)\n",
              h, ess, length(chain), mean(chain),
              quantile(chain, 0.025), quantile(chain, 0.975)))
}

# ---- Save results ---
saveRDS(fit.mcmc, "data/stage2_mcmc_results_long.rds")
cat("\nSaved to data/stage2_mcmc_results_long.rds\n")

# ---- Compute return levels ---
cat("\nComputing return levels...\n")
mu.chain <- exp(fit.mcmc$psi.selected)
sigma.chain <- exp(fit.mcmc$psi.selected + fit.mcmc$tau.selected)
xi.chain <- g(fit.mcmc$phi.selected)

station_meta <- stage1$station_meta

for (M in c(20, 50, 100)) {
  rl.chain <- mu.chain + sigma.chain / xi.chain *
    ((-log(1 - 1/M))^(-xi.chain) - 1)
  rl.posmean <- apply(rl.chain, 2, mean)
  rl.possd <- apply(rl.chain, 2, sd)
  assign(paste0("rl", M, ".mean"), rl.posmean)
  assign(paste0("rl", M, ".sd"), rl.possd)
  cat(sprintf("  %d-year: range %.1f - %.1f mm (mean), SD %.1f - %.1f\n",
              M, min(rl.posmean), max(rl.posmean), min(rl.possd), max(rl.possd)))
}

# Smoothed parameters
mu.posmean <- exp(fit.mcmc$psi.posmean)
sigma.posmean <- exp(fit.mcmc$psi.posmean + fit.mcmc$tau.posmean)
xi.posmean <- g(fit.mcmc$phi.posmean)

cat(sprintf("\nSmoothed parameters:\n  mu: %.1f - %.1f\n  sigma: %.1f - %.1f\n  xi: %.3f - %.3f\n",
            min(mu.posmean), max(mu.posmean),
            min(sigma.posmean), max(sigma.posmean),
            min(xi.posmean), max(xi.posmean)))

# ---- Trace plots ---
library(patchwork)

plots <- lapply(hyper_names, function(nm) {
  chain <- fit.mcmc[[nm]]
  df <- data.frame(iter = seq_along(chain), value = chain)
  ggplot(df, aes(iter, value)) +
    geom_line(alpha = 0.4, linewidth = 0.3) +
    labs(title = nm, x = "", y = "") +
    theme_minimal(base_size = 9)
})

p_trace <- wrap_plots(plots, ncol = 2) +
  plot_annotation(title = "MCMC trace plots (50k iterations, warm start)")
ggsave("figures/mcmc_trace_plots_long.pdf", p_trace, width = 12, height = 14)

# ---- MLE vs smoothed ----
mles.orig <- stage1$mles.original

comp_df <- data.frame(
  mu_mle = mles.orig[,1], mu_smooth = mu.posmean,
  sigma_mle = mles.orig[,2], sigma_smooth = sigma.posmean,
  xi_mle = mles.orig[,3], xi_smooth = xi.posmean
)

p_comp1 <- ggplot(comp_df, aes(mu_mle, mu_smooth)) +
  geom_point(alpha = 0.6) + geom_abline(linetype = "dashed", colour = "red") +
  labs(title = expression(mu), x = "MLE", y = "Posterior mean") + theme_minimal()

p_comp2 <- ggplot(comp_df, aes(sigma_mle, sigma_smooth)) +
  geom_point(alpha = 0.6) + geom_abline(linetype = "dashed", colour = "red") +
  labs(title = expression(sigma), x = "MLE", y = "Posterior mean") + theme_minimal()

p_comp3 <- ggplot(comp_df, aes(xi_mle, xi_smooth)) +
  geom_point(alpha = 0.6) + geom_abline(linetype = "dashed", colour = "red") +
  labs(title = expression(xi), x = "MLE", y = "Posterior mean") + theme_minimal()

p_comp <- p_comp1 + p_comp2 + p_comp3 +
  plot_annotation(title = "Max-and-Smooth: Site-wise MLE vs smoothed posterior mean (50k run)")
ggsave("figures/mle_vs_smooth_long.pdf", p_comp, width = 14, height = 5)

# ---- Return level maps ----
rl_df <- data.frame(
  station_meta,
  mu_smooth = mu.posmean,
  sigma_smooth = sigma.posmean,
  xi_smooth = xi.posmean,
  rl20_mean = rl20.mean, rl20_sd = rl20.sd,
  rl50_mean = rl50.mean, rl50_sd = rl50.sd,
  rl100_mean = rl100.mean, rl100_sd = rl100.sd
)

rl_sf <- st_as_sf(rl_df, coords = c("longitud", "latitud"), crs = 4326)

p20 <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl20_mean), size = 2.5) +
  scale_colour_viridis_c(option = "B") +
  labs(title = "20-year return level", colour = "mm") + theme_minimal()

p50 <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl50_mean), size = 2.5) +
  scale_colour_viridis_c(option = "B") +
  labs(title = "50-year return level", colour = "mm") + theme_minimal()

p100 <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl100_mean), size = 2.5) +
  scale_colour_viridis_c(option = "B") +
  labs(title = "100-year return level", colour = "mm") + theme_minimal()

p20sd <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl20_sd), size = 2.5) +
  scale_colour_viridis_c(option = "E") +
  labs(title = "20-year SD", colour = "mm") + theme_minimal()

p50sd <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl50_sd), size = 2.5) +
  scale_colour_viridis_c(option = "E") +
  labs(title = "50-year SD", colour = "mm") + theme_minimal()

p100sd <- ggplot(rl_sf) +
  geom_sf(aes(colour = rl100_sd), size = 2.5) +
  scale_colour_viridis_c(option = "E") +
  labs(title = "100-year SD", colour = "mm") + theme_minimal()

p_rl <- (p20 + p50 + p100) / (p20sd + p50sd + p100sd) +
  plot_annotation(
    title = "Max-and-Smooth: Return level estimates (top) and posterior SD (bottom)",
    subtitle = sprintf("%d stations, 50k MCMC iterations with warm start", ns)
  )
ggsave("figures/return_level_maps_long.pdf", p_rl, width = 16, height = 10)

saveRDS(rl_df, "data/return_levels_andalucia_long.rds")
write.csv(rl_df, "data/return_levels_andalucia_long.csv", row.names = FALSE)

cat("\n=== Long MCMC run complete ===\n")
