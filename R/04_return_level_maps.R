# 04_return_level_maps.R — 2x2 panel: return level posterior means
#
# Simple, clean figure: 10, 20, 50, 100 year return levels side by side.
# Run from project root: Rscript R/04_return_level_maps.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("1x4 Return level panel\n")
cat("========================================\n")

# ---- 1. Load data ----
cat("Loading data...\n")

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

loc  <- s1$loc
ns   <- nrow(loc)

eta_psi_draws <- s2$psi.selected
eta_tau_draws <- s2$tau.selected
eta_phi_draws <- s2$phi.selected
n_draws <- nrow(eta_psi_draws)

# Covariate model: beta_psi_draws is n_draws x n_cov_psi matrix
beta_psi_draws <- s2$beta_psi_draws   # n_draws x 4
mu_tau_draws   <- s2$beta_tau         # length n_draws
mu_phi_draws   <- s2$beta_phi         # length n_draws
cov_std        <- s2$cov_standardisation

# phi has no GP — only intercept + iid noise
phi_has_gp     <- isTRUE(s2$phi_has_gp)
nugget_phi_draws <- s2$nugget_phi     # length n_draws (iid noise SD)

cat(sprintf("  %d draws, %d stations\n", n_draws, ns))
cat(sprintf("  Covariates: %d coefficients for psi\n", ncol(beta_psi_draws)))
cat(sprintf("  phi has GP: %s\n", phi_has_gp))

# ---- 2. Load grid covariates & build design matrices ----
cat("Loading grid covariates...\n")

grid_cov <- readRDS("data/grid_covariates.rds")
pred_pts <- data.frame(lon = grid_cov$lon, lat = grid_cov$lat)
n_pred   <- nrow(pred_pts)

# Station-level design matrix (same standardisation)
stn_exp <- readRDS("data/station_exposure.rds")
s1_ids  <- s1$station_meta$indicativo
exp_idx <- match(s1_ids, stn_exp$indicativo)
alt_std_stn <- (stn_exp$alt_dem[exp_idx] - cov_std$alt_mean) / cov_std$alt_sd
exp_std_stn <- (stn_exp$exposure_mean[exp_idx] - cov_std$exp_mean) / cov_std$exp_sd
X_stn <- cbind(1, alt_std_stn, exp_std_stn, alt_std_stn * exp_std_stn)

# Standardise grid covariates with Clausius-Clapeyron attenuation above
# highest station.  Daily precipitation extremes increase with altitude
# (Formetta et al. 2022: 7.5-10 % per 1000 m), but the enhancement must
# taper as precipitable water decreases ~exp(-h/H_w).  We use the
# atmospheric moisture scale height H_w = 2000 m (standard mid-latitude
# value) to attenuate the altitude effect above the highest station.
h_peak <- max(stn_exp$alt_dem[exp_idx])   # highest station altitude (m)
H_w    <- 2000                            # moisture scale height (m)

alt_raw_pred <- grid_cov$alt_dem
above <- pmax(0, alt_raw_pred - h_peak)
decay <- exp(-above / H_w)
# Effective altitude: actual altitude up to h_peak, then attenuated above
alt_eff <- ifelse(alt_raw_pred <= h_peak, alt_raw_pred,
                  h_peak + above * decay)
alt_std_pred <- (alt_eff - cov_std$alt_mean) / cov_std$alt_sd

exp_std_pred <- (grid_cov$exposure_mean - cov_std$exp_mean) / cov_std$exp_sd
# Clamp exposure only (altitude handled by attenuation)
exp_std_pred <- pmin(pmax(exp_std_pred, min(exp_std_stn)), max(exp_std_stn))
X_pred <- cbind(1, alt_std_pred, exp_std_pred, alt_std_pred * exp_std_pred)

cat(sprintf("  Altitude attenuation: h_peak=%dm, H_w=%dm\n", h_peak, H_w))
cat(sprintf("  Alt_eff at summit (%.0fm): %.0fm  (decay=%.2f)\n",
            max(alt_raw_pred), max(alt_eff), min(decay[above > 0])))

cat(sprintf("  Grid: %d points, X_pred: %d x %d\n", n_pred, nrow(X_pred), ncol(X_pred)))

# ---- 3. Kriging weights (psi and tau only — phi is iid) ----
cat("Computing kriging weights (psi and tau)...\n")

pred_res <- grid_cov$lon[2] - grid_cov$lon[1]
if (is.na(pred_res) || pred_res <= 0) pred_res <- 0.025

matern52 <- function(d, sigma2, phi) {
  s5 <- sqrt(5) * d / phi
  sigma2 * (1 + s5 + s5^2 / 3) * exp(-s5)
}

# GP hyperparameters (posterior means) — only for psi and tau
sigma_gp <- c(mean(s2$sigma_gp_psi), mean(s2$sigma_gp_tau))
phi_gp   <- c(mean(s2$phi_gp_psi), mean(s2$phi_gp_tau))
nugget   <- c(mean(s2$nugget_psi), mean(s2$nugget_tau))

dist_ss <- as.matrix(dist(loc))
# Compute prediction-to-station distances directly (avoids huge intermediate matrix)
dist_pg <- matrix(0, n_pred, ns)
for (j in seq_len(ns)) {
  dist_pg[, j] <- sqrt((pred_pts$lon - loc[j, "lon"])^2 +
                        (pred_pts$lat - loc[j, "lat"])^2)
}

compute_kriging <- function(k) {
  sig2 <- sigma_gp[k]^2
  nug2 <- nugget[k]^2
  phi  <- phi_gp[k]
  Sigma <- matern52(dist_ss, sig2, phi)
  diag(Sigma) <- diag(Sigma) + nug2
  Sigma_inv <- solve(Sigma)
  gamma_mat <- matern52(dist_pg, sig2, phi)
  W <- gamma_mat %*% Sigma_inv
  cond_var <- sig2 - rowSums(W * gamma_mat)
  cond_var <- pmax(cond_var, 1e-10)
  list(W = W, cond_sd = sqrt(cond_var))
}

krig <- lapply(1:2, compute_kriging)  # only psi and tau

# ---- 4. Return levels at T = 10, 20, 50, 100 ----
cat("Computing return levels...\n")

return_periods <- c(10, 20, 50, 100)
n_rp <- length(return_periods)

rl_mean <- matrix(0, n_pred, n_rp)
rl_sd   <- matrix(0, n_pred, n_rp)

# GP residuals: eta - covariate mean
# psi: matrix mean (covariates)
mean_psi_at_stations <- beta_psi_draws %*% t(X_stn)   # n_draws x ns
resid_psi <- eta_psi_draws - mean_psi_at_stations
# tau: scalar intercept
resid_tau <- sweep(eta_tau_draws, 1, mu_tau_draws, "-")
# phi: no GP residuals (iid noise, not spatially correlated)

chunk_size <- 2000
n_chunks <- ceiling(n_pred / chunk_size)
set.seed(42)

for (ch in seq_len(n_chunks)) {
  idx_start <- (ch - 1) * chunk_size + 1
  idx_end   <- min(ch * chunk_size, n_pred)
  idx <- idx_start:idx_end
  nc <- length(idx)

  if (ch %% 5 == 1 || ch == n_chunks)
    cat(sprintf("  Chunk %d/%d...\n", ch, n_chunks))

  # Predict psi: covariate mean + kriged GP residual + conditional noise
  predict_psi_cov <- function() {
    W_ch <- krig[[1]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[1]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    kriged <- t(W_ch %*% t(resid_psi))                  # n_draws x nc
    mean_pred <- beta_psi_draws %*% t(X_pred[idx, , drop = FALSE])  # n_draws x nc
    pred <- mean_pred + kriged
    pred <- pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
    pred
  }

  # Predict tau: scalar intercept + kriged GP residual + conditional noise
  predict_tau <- function() {
    W_ch <- krig[[2]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[2]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    pred <- sweep(t(W_ch %*% t(resid_tau)), 1, mu_tau_draws, "+")
    pred <- pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
    pred
  }

  # Predict phi: intercept + iid noise (no spatial kriging)
  predict_phi_iid <- function() {
    noise_sd <- matrix(nugget_phi_draws, nrow = n_draws, ncol = nc)
    sweep(matrix(rnorm(n_draws * nc), n_draws, nc) * noise_sd,
          1, mu_phi_draws, "+")
  }

  psi_ch <- predict_psi_cov()
  tau_ch <- predict_tau()
  phi_ch <- predict_phi_iid()

  mu_gev    <- exp(psi_ch)
  sigma_gev <- exp(psi_ch + tau_ch)
  xi_gev    <- g(phi_ch)

  for (r in seq_len(n_rp)) {
    Tr <- return_periods[r]
    rl_ch <- mu_gev + sigma_gev / xi_gev * ((-log(1 - 1 / Tr))^(-xi_gev) - 1)
    rl_mean[idx, r] <- colMeans(rl_ch)
    rl_sd[idx, r]   <- apply(rl_ch, 2, sd)
  }
}

# ---- 5. Build data frames ----
grid_list <- list()
for (r in seq_len(n_rp)) {
  grid_list[[r]] <- data.frame(
    lon = pred_pts$lon, lat = pred_pts$lat,
    rl_mean = rl_mean[, r],
    rl_sd   = rl_sd[, r],
    rp = return_periods[r]
  )
}
grid_df <- do.call(rbind, grid_list)

# Station-level return levels
stn_list <- list()
for (r in seq_len(n_rp)) {
  Tr <- return_periods[r]
  mu_stn    <- exp(eta_psi_draws)
  sigma_stn <- exp(eta_psi_draws + eta_tau_draws)
  xi_stn    <- g(eta_phi_draws)
  rl_stn_draws <- mu_stn + sigma_stn / xi_stn * ((-log(1 - 1 / Tr))^(-xi_stn) - 1)
  stn_list[[r]] <- data.frame(
    lon = loc[, "lon"], lat = loc[, "lat"],
    rl_mean = colMeans(rl_stn_draws),
    rp = return_periods[r]
  )
}
stn_df <- do.call(rbind, stn_list)

for (r in seq_len(n_rp)) {
  rr <- grid_df[grid_df$rp == return_periods[r], ]
  cat(sprintf("  RL%d: %.1f - %.1f mm\n",
    return_periods[r], min(rr$rl_mean), max(rr$rl_mean)))
}

# ---- 6. Plot 1x4 ----
cat("Generating 1x4 panel...\n")

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

rl_min <- min(grid_df$rl_mean)
rl_max <- max(grid_df$rl_mean)

# AEMET alarm thresholds (precip in 12h, typical for Andalucia)
alarm_thresholds <- c(80, 120)        # orange, red (mm)
alarm_colours    <- c("#FF8C00", "red") # orange, red

legend_breaks <- seq(0, ceiling(rl_max / 50) * 50, by = 50)

make_panel <- function(rp_val) {
  g_sub <- grid_df[grid_df$rp == rp_val, ]
  s_sub <- stn_df[stn_df$rp == rp_val, ]
  data_range <- range(g_sub$rl_mean, na.rm = TRUE)

  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = g_sub, aes(x = lon, y = lat, fill = rl_mean),
              width = pred_res, height = pred_res) +
    scale_fill_viridis_c(option = "B", name = "mm",
                         limits = c(rl_min, rl_max),
                         breaks = legend_breaks)

  p +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = paste0(rp_val, "-year return level")) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.key.height = unit(1.0, "cm"),
          plot.margin = margin(2, 2, 2, 2))
}

p_all <- make_panel(10) + make_panel(20) + make_panel(50) + make_panel(100) +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = "Maximum daily rainfall return levels: posterior mean",
    subtitle = sprintf(
      "%d stations | %d grid points | %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"),
                  plot.margin = margin(2, 2, 2, 2))
  )

ggsave("figures/return_level_maps.png",
       p_all, width = 14, height = 8.5, dpi = 300, bg = "white")

cat("Saved figures/return_level_maps.png\n")

# ---- 7. SD panel ----
cat("Generating SD panel...\n")

sd_min <- min(grid_df$rl_sd)
sd_max <- max(grid_df$rl_sd)
sd_breaks <- seq(0, ceiling(sd_max / 10) * 10, by = 10)

make_sd_panel <- function(rp_val) {
  g_sub <- grid_df[grid_df$rp == rp_val, ]

  ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = g_sub, aes(x = lon, y = lat, fill = rl_sd),
              width = pred_res, height = pred_res) +
    scale_fill_viridis_c(option = "A", name = "SD (mm)",
                         limits = c(sd_min, sd_max),
                         breaks = sd_breaks) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = paste0(rp_val, "-year return level SD")) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.key.height = unit(1.0, "cm"),
          plot.margin = margin(2, 2, 2, 2))
}

p_sd <- make_sd_panel(10) + make_sd_panel(20) + make_sd_panel(50) + make_sd_panel(100) +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = "Maximum daily rainfall return levels: posterior standard deviation",
    subtitle = sprintf(
      "%d stations | %d grid points | %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"),
                  plot.margin = margin(2, 2, 2, 2))
  )

ggsave("figures/return_level_maps_sd.png",
       p_sd, width = 14, height = 8.5, dpi = 300, bg = "white")

cat("Saved figures/return_level_maps_sd.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
