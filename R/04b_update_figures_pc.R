# 04b_update_figures_pc.R — Diagnostic figures for Matérn GP + PC priors model
# Run from project root: Rscript R/04b_update_figures_pc.R

library(ggplot2)
library(sf)
library(coda)
library(patchwork)

source("vendor/max_and_smooth/stage1_functions.R")

# Load saved results
s1  <- readRDS("data/stage1_results.rds")
fit <- readRDS("data/stage2_matern_pc_results.rds")

loc        <- s1$loc
meta       <- s1$station_meta
mles.orig  <- s1$mles.original
mles       <- s1$mles
ns         <- nrow(loc)

# Posterior parameter chains
mu.chain    <- exp(fit$psi.selected)
sigma.chain <- exp(fit$psi.selected + fit$tau.selected)
xi.chain    <- g(fit$phi.selected)

# Posterior means
mu.pm    <- exp(fit$psi.posmean)
sigma.pm <- exp(fit$psi.posmean + fit$tau.posmean)
xi.pm    <- g(fit$phi.posmean)

cat("Stations:", ns, "\n")
cat("MCMC samples:", nrow(mu.chain), "\n")

# =============================================================================
# 1. MCMC TRACE PLOTS (GP hyperparameters)
# =============================================================================
cat("Generating trace plots...\n")

hyper_names <- c("beta_psi", "beta_tau", "beta_phi",
                 "sigma_gp_psi", "sigma_gp_tau", "sigma_gp_phi",
                 "phi_gp_psi", "phi_gp_tau", "phi_gp_phi",
                 "nugget_psi", "nugget_tau", "nugget_phi")

pdf("figures/mcmc_trace_plots_pc.pdf", width = 14, height = 14)
par(mfrow = c(6, 2), mar = c(3, 4, 2, 1))
for (h in hyper_names) {
  chain <- fit[[h]]
  ess <- round(effectiveSize(chain))
  plot(chain, type = "l", col = "steelblue",
       main = sprintf("%s  (ESS = %d / %d)", h, ess, length(chain)),
       xlab = "Iteration", ylab = h, cex.main = 1.1)
  abline(h = mean(chain), col = "red", lwd = 2)
  abline(h = quantile(chain, c(0.025, 0.975)), col = "red", lty = 2)
}
dev.off()

# =============================================================================
# 2. MLE vs SMOOTHED SCATTER PLOTS
# =============================================================================
cat("Generating MLE vs smoothed scatter...\n")

df_compare <- data.frame(
  mu_mle    = mles.orig[, "mu"],
  sigma_mle = mles.orig[, "sigma"],
  xi_mle    = mles.orig[, "xi"],
  mu_sm     = mu.pm,
  sigma_sm  = sigma.pm,
  xi_sm     = xi.pm
)

p1 <- ggplot(df_compare, aes(mu_mle, mu_sm)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(x = expression(hat(mu)[MLE]), y = expression(hat(mu)[smooth]),
       title = expression(mu ~ "(location)")) +
  theme_minimal(base_size = 13)

p2 <- ggplot(df_compare, aes(sigma_mle, sigma_sm)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(x = expression(hat(sigma)[MLE]), y = expression(hat(sigma)[smooth]),
       title = expression(sigma ~ "(scale)")) +
  theme_minimal(base_size = 13)

p3 <- ggplot(df_compare, aes(xi_mle, xi_sm)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(x = expression(hat(xi)[MLE]), y = expression(hat(xi)[smooth]),
       title = expression(xi ~ "(shape)")) +
  theme_minimal(base_size = 13)

p_scatter <- p1 + p2 + p3 +
  plot_annotation(title = "MLE vs Bayesian-smoothed parameter estimates (PC priors)",
                  subtitle = sprintf("%d stations, Max-and-Smooth with Matern(5/2) GP + PC priors", ns))
ggsave("figures/mle_vs_smooth_pc.pdf", p_scatter, width = 14, height = 5)

# =============================================================================
# 3. RETURN LEVEL MAPS
# =============================================================================
cat("Generating return level maps...\n")

rl_df <- data.frame(meta)
for (M in c(20, 50, 100)) {
  rl.chain <- mu.chain + sigma.chain / xi.chain *
    ((-log(1 - 1/M))^(-xi.chain) - 1)
  rl_df[[paste0("rl", M, "_mean")]] <- apply(rl.chain, 2, mean)
  rl_df[[paste0("rl", M, "_sd")]]   <- apply(rl.chain, 2, sd)
}

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
    title = "Return level estimates (top) and posterior SD (bottom) — PC priors",
    subtitle = sprintf("%d stations, Matern(5/2) GP + PC priors", ns)
  )
ggsave("figures/return_level_maps_pc.pdf", p_rl, width = 16, height = 10)

# =============================================================================
# 4. XI MLE vs SMOOTHED MAP (side-by-side)
# =============================================================================
cat("Generating xi comparison map...\n")

xi_df <- data.frame(meta,
                    xi_mle = mles.orig[, "xi"],
                    xi_smooth = xi.pm)
xi_sf <- st_as_sf(xi_df, coords = c("longitud", "latitud"), crs = 4326)

xi_range <- range(c(xi_df$xi_mle, xi_df$xi_smooth), na.rm = TRUE)

p_xi_mle <- ggplot(xi_sf) +
  geom_sf(aes(colour = xi_mle), size = 2.5) +
  scale_colour_viridis_c(option = "D", limits = xi_range) +
  labs(title = expression("MLE " * hat(xi)), colour = expression(xi)) +
  theme_minimal(base_size = 13)

p_xi_sm <- ggplot(xi_sf) +
  geom_sf(aes(colour = xi_smooth), size = 2.5) +
  scale_colour_viridis_c(option = "D", limits = xi_range) +
  labs(title = expression("Smoothed " * hat(xi)), colour = expression(xi)) +
  theme_minimal(base_size = 13)

p_xi_compare <- p_xi_mle + p_xi_sm +
  plot_annotation(
    title = expression("Shape parameter " * xi * ": MLE vs Bayesian-smoothed (PC priors)"),
    subtitle = sprintf("%d stations", ns)
  )
ggsave("figures/xi_mle_vs_smooth_map_pc.pdf", p_xi_compare, width = 12, height = 5)

# =============================================================================
# 5. SMOOTHED PARAMETER MAPS (posterior means)
# =============================================================================
cat("Generating smoothed parameter maps...\n")

smooth_df <- data.frame(meta,
                        mu = mu.pm,
                        sigma = sigma.pm,
                        xi = xi.pm)
smooth_sf <- st_as_sf(smooth_df, coords = c("longitud", "latitud"), crs = 4326)

p_mu_s <- ggplot(smooth_sf) +
  geom_sf(aes(colour = mu), size = 2.5) +
  scale_colour_viridis_c(option = "C") +
  labs(title = expression("Smoothed " * hat(mu)), colour = "mm") +
  theme_minimal()

p_sigma_s <- ggplot(smooth_sf) +
  geom_sf(aes(colour = sigma), size = 2.5) +
  scale_colour_viridis_c(option = "C") +
  labs(title = expression("Smoothed " * hat(sigma)), colour = "mm") +
  theme_minimal()

p_xi_s <- ggplot(smooth_sf) +
  geom_sf(aes(colour = xi), size = 2.5) +
  scale_colour_viridis_c(option = "D") +
  labs(title = expression("Smoothed " * hat(xi)), colour = "") +
  theme_minimal()

p_smooth <- p_mu_s + p_sigma_s + p_xi_s +
  plot_annotation(title = "Bayesian-smoothed GEV parameters (PC priors)",
                  subtitle = sprintf("%d stations, Matern(5/2) GP + PC priors", ns))
ggsave("figures/smoothed_param_maps_pc.pdf", p_smooth, width = 14, height = 5)

cat("All PC prior figures generated.\n")
