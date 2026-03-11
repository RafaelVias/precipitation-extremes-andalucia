# 05_exceedance_maps.R — 3x3 panel: exceedance probability (posterior mean)
#
# Columns: T = 20, 50, 100 year planning horizons
# Rows: rainfall thresholds = 100, 150, 200 mm/day
# Each cell: P(annual max > threshold at least once in T years)
#            = 1 - (1 - p_annual)^T  where p_annual = P(annual max > threshold)
#
# Run from project root: Rscript R/05_exceedance_maps.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Plots 2-3: Exceedance probability panels\n")
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

# phi GP flag
phi_has_gp <- isTRUE(s2$phi_has_gp)

cat("  Stations:", ns, "| Draws:", n_draws, "\n")
cat(sprintf("  Covariates: %d coefficients for psi\n", ncol(beta_psi_draws)))
cat(sprintf("  phi has GP: %s\n", phi_has_gp))

# ---- 2. Load grid covariates & build design matrices ----
cat("Loading grid covariates...\n")

grid_cov <- readRDS("data/grid_covariates.rds")
pred_pts <- data.frame(lon = grid_cov$lon, lat = grid_cov$lat)
n_pred   <- nrow(pred_pts)

pred_res <- abs(sort(unique(grid_cov$lon))[2] - sort(unique(grid_cov$lon))[1])

# Station-level design matrix (same standardisation)
stn_exp <- readRDS("data/station_exposure.rds")
s1_ids  <- s1$station_meta$indicativo
exp_idx <- match(s1_ids, stn_exp$indicativo)
alt_std_stn <- (stn_exp$alt_dem[exp_idx] - cov_std$alt_mean) / cov_std$alt_sd
exp_std_stn <- (stn_exp$exposure_mean[exp_idx] - cov_std$exp_mean) / cov_std$exp_sd
X_stn <- cbind(1, alt_std_stn, exp_std_stn, alt_std_stn * exp_std_stn)

# Clausius-Clapeyron attenuation of altitude effect above highest station
# (see Formetta et al. 2022; H_w = moisture scale height)
h_peak <- max(stn_exp$alt_dem[exp_idx])
H_w    <- 2000
alt_raw_pred <- grid_cov$alt_dem
above <- pmax(0, alt_raw_pred - h_peak)
decay <- exp(-above / H_w)
alt_eff <- ifelse(alt_raw_pred <= h_peak, alt_raw_pred, h_peak + above * decay)
alt_std_pred <- (alt_eff - cov_std$alt_mean) / cov_std$alt_sd

exp_std_pred <- (grid_cov$exposure_mean - cov_std$exp_mean) / cov_std$exp_sd
exp_std_pred <- pmin(pmax(exp_std_pred, min(exp_std_stn)), max(exp_std_stn))
X_pred <- cbind(1, alt_std_pred, exp_std_pred, alt_std_pred * exp_std_pred)

cat("  Grid:", n_pred, "points at", pred_res, "deg\n")

# ---- 3. Kriging weights ----
cat("Computing kriging weights...\n")

matern52 <- function(d, sigma2, phi) {
  s5 <- sqrt(5) * d / phi
  sigma2 * (1 + s5 + s5^2 / 3) * exp(-s5)
}

# GP hyperparameters (posterior means) — all 3 GPs
sigma_gp <- c(mean(s2$sigma_gp_psi), mean(s2$sigma_gp_tau), mean(s2$sigma_gp_phi))
phi_gp   <- c(mean(s2$phi_gp_psi), mean(s2$phi_gp_tau), mean(s2$phi_gp_phi))
nugget   <- c(mean(s2$nugget_psi), mean(s2$nugget_tau), mean(s2$nugget_phi))

dist_ss <- as.matrix(dist(loc))
# Compute prediction-to-station distances directly (avoids huge intermediate matrix)
dist_pg <- matrix(0, n_pred, ns)
for (j in seq_len(ns)) {
  dist_pg[, j] <- sqrt((pred_pts$lon - loc[j, "lon"])^2 +
                        (pred_pts$lat - loc[j, "lat"])^2)
}

compute_cond_gp <- function(k) {
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

krig <- lapply(1:3, compute_cond_gp)  # psi, tau, phi
cat("  Conditional GP weights computed (psi, tau, phi)\n")

# ---- 4. Compute exceedance probabilities ----
cat("Computing exceedance probabilities at grid points...\n")

thresholds <- c(100, 150, 200)       # mm/day
horizons   <- c(20, 50, 100)         # years
n_thr <- length(thresholds)
n_hor <- length(horizons)

# For each grid point and threshold, we need P(annual max > threshold)
# from each posterior draw, then convert to T-year probability.
# Store: p_annual_mean[n_pred, n_thr], p_annual_sd[n_pred, n_thr]
# Then derive horizon probabilities analytically.

# But to get the SD of the T-year probability correctly, we need to
# propagate through 1-(1-p)^T for each draw. So store per-draw
# annual exceedance probability summaries.

# Strategy: for each chunk of grid points, compute p_annual for each draw,
# then for each horizon compute 1-(1-p)^T, and accumulate mean + sd.

# Allocators: n_pred x n_thr x n_hor
exc_mean <- array(0, dim = c(n_pred, n_thr, n_hor))
exc_sd   <- array(0, dim = c(n_pred, n_thr, n_hor))

# GP residuals: eta - covariate/intercept mean
mean_psi_at_stations <- beta_psi_draws %*% t(X_stn)   # n_draws × ns
resid_psi <- eta_psi_draws - mean_psi_at_stations
resid_tau <- sweep(eta_tau_draws, 1, mu_tau_draws, "-")
resid_phi <- sweep(eta_phi_draws, 1, mu_phi_draws, "-")

chunk_size <- 2000
n_chunks <- ceiling(n_pred / chunk_size)
set.seed(42)

for (ch in seq_len(n_chunks)) {
  idx_start <- (ch - 1) * chunk_size + 1
  idx_end   <- min(ch * chunk_size, n_pred)
  idx <- idx_start:idx_end
  nc <- length(idx)

  if (ch %% 5 == 1 || ch == n_chunks)
    cat(sprintf("  Chunk %d/%d (points %d-%d)...\n", ch, n_chunks, idx_start, idx_end))

  # Predict psi: covariate mean + kriged GP residual + conditional noise
  predict_psi_cov <- function() {
    W_ch <- krig[[1]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[1]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    kriged <- t(W_ch %*% t(resid_psi))
    mean_pred <- beta_psi_draws %*% t(X_pred[idx, , drop = FALSE])
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

  # Predict phi: scalar intercept + conditional GP residual + conditional noise
  predict_phi <- function() {
    W_ch <- krig[[3]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[3]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    pred <- sweep(t(W_ch %*% t(resid_phi)), 1, mu_phi_draws, "+")
    pred <- pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
    pred
  }

  psi_ch <- predict_psi_cov()
  tau_ch <- predict_tau()
  phi_ch <- predict_phi()

  mu_gev    <- exp(psi_ch)                  # n_draws x nc
  sigma_gev <- exp(psi_ch + tau_ch)
  xi_gev    <- g(phi_ch)

  # GEV exceedance probability: P(Z > x) = 1 - GEV_cdf(x)
  # GEV_cdf(x) = exp(-(1 + xi*(x - mu)/sigma)^(-1/xi))  when 1 + xi*(x-mu)/sigma > 0
  for (t in seq_len(n_thr)) {
    x <- thresholds[t]
    z <- 1 + xi_gev * (x - mu_gev) / sigma_gev

    # Where z <= 0, x is beyond the support: p_annual = 0 (for xi < 0)
    # or in the extreme tail. For xi > 0, z should always be > 0.
    z <- pmax(z, 1e-10)

    p_annual <- 1 - exp(-z^(-1 / xi_gev))   # n_draws x nc

    for (h in seq_len(n_hor)) {
      T_h <- horizons[h]
      p_horizon <- 1 - (1 - p_annual)^T_h
      exc_mean[idx, t, h] <- colMeans(p_horizon)
      exc_sd[idx, t, h]   <- apply(p_horizon, 2, sd)
    }
  }
}

cat("  Exceedance probabilities computed\n")

# Print diagnostic ranges
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    cat(sprintf("  P(max > %dmm in %dy): mean %.3f-%.3f | SD %.3f-%.3f\n",
      thresholds[t], horizons[h],
      min(exc_mean[, t, h]), max(exc_mean[, t, h]),
      min(exc_sd[, t, h]), max(exc_sd[, t, h])))
  }
}

# ---- 5. Assemble data frames ----
cat("Assembling plot data...\n")

grid_rows <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      lon = pred_pts$lon, lat = pred_pts$lat,
      prob_mean = exc_mean[, t, h],
      prob_sd = exc_sd[, t, h],
      threshold = thresholds[t],
      horizon = horizons[h]
    )
  }
}
grid_df <- do.call(rbind, grid_rows)

grid_df$thr_label <- paste0("> ", grid_df$threshold, " mm/day")
grid_df$thr_label <- factor(grid_df$thr_label,
  levels = paste0("> ", thresholds, " mm/day"))

grid_df$hor_label <- paste0(grid_df$horizon, "-year horizon")
grid_df$hor_label <- factor(grid_df$hor_label,
  levels = paste0(horizons, "-year horizon"))

# ---- 6. Plot 2: Posterior mean exceedance probability ----
cat("Generating Plot 2 (exceedance probability mean)...\n")

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)
neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

tile_res <- pred_res

base_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        legend.key.height = unit(0.7, "cm"),
        plot.margin = margin(2, 2, 2, 2))

make_panel <- function(data, fill_var, fill_label, title, option = "D",
                       limits = c(0, 1), oob = scales::squish,
                       contour_breaks = NULL) {
  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    scale_fill_viridis_c(option = option, name = fill_label,
                         limits = limits, oob = oob)

  p + geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

# Build 3x3 grid: rows = thresholds, columns = horizons
panels_mean <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    g_sub <- grid_df[grid_df$threshold == thresholds[t] &
                     grid_df$horizon == horizons[h], ]
    title <- paste0("P(max > ", thresholds[t], " mm in ", horizons[h], " yr)")
    panels_mean[[length(panels_mean) + 1]] <- make_panel(
      g_sub, "prob_mean", "prob", title, option = "B",
      contour_breaks = c(0.25, 0.5, 0.75))
  }
}

p_mean <- (panels_mean[[1]] + panels_mean[[2]] + panels_mean[[3]]) /
          (panels_mean[[4]] + panels_mean[[5]] + panels_mean[[6]]) /
          (panels_mean[[7]] + panels_mean[[8]] + panels_mean[[9]]) +
  plot_annotation(
    title = "Exceedance probability: posterior mean",
    subtitle = sprintf(
      "%d stations | %d grid points | %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"),
                  plot.margin = margin(2, 2, 2, 2))
  )

ggsave("figures/exceedance_prob.png", p_mean, width = 16, height = 10, dpi = 300, bg = "white")
cat("  Saved figures/exceedance_prob.png\n")

# ---- 7. Plot 3: Posterior SD exceedance probability ----
cat("Generating Plot 3 (exceedance probability SD)...\n")

sd_max <- max(grid_df$prob_sd)

panels_sd <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    g_sub <- grid_df[grid_df$threshold == thresholds[t] &
                     grid_df$horizon == horizons[h], ]
    title <- paste0("SD: P(max > ", thresholds[t], " mm in ", horizons[h], " yr)")
    panels_sd[[length(panels_sd) + 1]] <- make_panel(
      g_sub, "prob_sd", "SD", title, option = "A",
      limits = c(0, sd_max), contour_breaks = NULL)
  }
}

p_sd <- (panels_sd[[1]] + panels_sd[[2]] + panels_sd[[3]]) /
        (panels_sd[[4]] + panels_sd[[5]] + panels_sd[[6]]) /
        (panels_sd[[7]] + panels_sd[[8]] + panels_sd[[9]]) +
  plot_annotation(
    title = "Exceedance probability: posterior standard deviation",
    subtitle = sprintf(
      "%d stations | %d grid points | %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"),
                  plot.margin = margin(2, 2, 2, 2))
  )

ggsave("figures/exceedance_prob_sd.png", p_sd, width = 16, height = 10, dpi = 300, bg = "white")
cat("  Saved figures/exceedance_prob_sd.png\n")

cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
