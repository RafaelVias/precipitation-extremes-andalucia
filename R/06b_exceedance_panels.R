# 06b_exceedance_panels.R — 3x3 panels: exceedance probability + SD
#
# Columns: T = 20, 50, 100 year planning horizons
# Rows: rainfall thresholds = 100, 150, 200 mm/day
# Each cell: P(annual max > threshold at least once in T years)
#            = 1 - (1 - p_annual)^T  where p_annual = P(annual max > threshold)
#
# Produces two figures:
#   Plot 2: posterior mean exceedance probability
#   Plot 3: posterior SD of exceedance probability
#
# Run from project root: Rscript R/06b_exceedance_panels.R

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

mu_psi_draws <- s2$beta_psi
mu_tau_draws <- s2$beta_tau
mu_phi_draws <- s2$beta_phi

cat("  Stations:", ns, "| Draws:", n_draws, "\n")

# ---- 2. Prediction grid ----
cat("Creating prediction grid...\n")

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

pred_res <- 0.025
bbox <- st_bbox(andalucia)
pred_pts <- expand.grid(
  lon = seq(bbox["xmin"] + pred_res / 2, bbox["xmax"], by = pred_res),
  lat = seq(bbox["ymin"] + pred_res / 2, bbox["ymax"], by = pred_res)
)
pred_sf <- st_as_sf(pred_pts, coords = c("lon", "lat"), crs = 4326)
inside <- st_intersects(pred_sf, andalucia, sparse = FALSE)[, 1]
pred_pts <- pred_pts[inside, ]
n_pred <- nrow(pred_pts)

cat("  Grid:", n_pred, "points at", pred_res, "deg\n")

# ---- 3. Kriging weights ----
cat("Computing kriging weights...\n")

matern52 <- function(d, sigma2, phi) {
  s5 <- sqrt(5) * d / phi
  sigma2 * (1 + s5 + s5^2 / 3) * exp(-s5)
}

sigma_gp <- c(mean(s2$sigma_gp_psi), mean(s2$sigma_gp_tau), mean(s2$sigma_gp_phi))
phi_gp   <- c(mean(s2$phi_gp_psi), mean(s2$phi_gp_tau), mean(s2$phi_gp_phi))
nugget   <- c(mean(s2$nugget_psi), mean(s2$nugget_tau), mean(s2$nugget_phi))

dist_ss <- as.matrix(dist(loc))
dist_pg <- as.matrix(
  dist(rbind(as.matrix(pred_pts), as.matrix(data.frame(lon = loc[, "lon"], lat = loc[, "lat"]))))
)[1:n_pred, (n_pred + 1):(n_pred + ns)]

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

krig <- lapply(1:3, compute_kriging)
cat("  Kriging weights computed\n")

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

resid_psi <- sweep(eta_psi_draws, 1, mu_psi_draws, "-")
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

  predict_param <- function(k, resid, mu_draws) {
    W_ch <- krig[[k]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[k]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    pred <- sweep(t(W_ch %*% t(resid)), 1, mu_draws, "+")
    pred <- pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
    pred
  }

  psi_ch <- predict_param(1, resid_psi, mu_psi_draws)
  tau_ch <- predict_param(2, resid_tau, mu_tau_draws)
  phi_ch <- predict_param(3, resid_phi, mu_phi_draws)

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
        legend.key.height = unit(0.7, "cm"))

make_panel <- function(data, fill_var, fill_label, title, option = "D",
                       limits = c(0, 1), oob = scales::squish) {
  ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    scale_fill_viridis_c(option = option, name = fill_label,
                         limits = limits, oob = oob) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
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
      g_sub, "prob_mean", "prob", title, option = "D")
  }
}

p_mean <- (panels_mean[[1]] + panels_mean[[2]] + panels_mean[[3]]) /
          (panels_mean[[4]] + panels_mean[[5]] + panels_mean[[6]]) /
          (panels_mean[[7]] + panels_mean[[8]] + panels_mean[[9]]) +
  plot_annotation(
    title = "Exceedance probability: posterior mean",
    subtitle = sprintf(
      "Mat\u00e9rn(5/2) GP + PC priors | %d stations, %d grid points, %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/exceedance_prob_mean.pdf", p_mean, width = 16, height = 13, bg = "white")
ggsave("figures/exceedance_prob_mean.png", p_mean, width = 16, height = 13, dpi = 200, bg = "white")
cat("  Saved figures/exceedance_prob_mean.pdf and .png\n")

# ---- 7. Plot 3: Posterior SD of exceedance probability ----
cat("Generating Plot 3 (exceedance probability SD)...\n")

panels_sd <- list()
sd_max <- max(grid_df$prob_sd) * 1.02

for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    g_sub <- grid_df[grid_df$threshold == thresholds[t] &
                     grid_df$horizon == horizons[h], ]
    title <- paste0("SD: P(max > ", thresholds[t], " mm in ", horizons[h], " yr)")
    panels_sd[[length(panels_sd) + 1]] <- make_panel(
      g_sub, "prob_sd", "SD", title, option = "A",
      limits = c(0, sd_max))
  }
}

p_sd <- (panels_sd[[1]] + panels_sd[[2]] + panels_sd[[3]]) /
        (panels_sd[[4]] + panels_sd[[5]] + panels_sd[[6]]) /
        (panels_sd[[7]] + panels_sd[[8]] + panels_sd[[9]]) +
  plot_annotation(
    title = "Exceedance probability: posterior uncertainty (SD)",
    subtitle = sprintf(
      "Mat\u00e9rn(5/2) GP + PC priors | %d stations, %d grid points, %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/exceedance_prob_sd.pdf", p_sd, width = 16, height = 13, bg = "white")
ggsave("figures/exceedance_prob_sd.png", p_sd, width = 16, height = 13, dpi = 200, bg = "white")
cat("  Saved figures/exceedance_prob_sd.pdf and .png\n")

cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
