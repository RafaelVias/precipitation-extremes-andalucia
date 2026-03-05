# 05_interpolation_maps.R — Posterior predictive maps of Andalucía
#
# Uses the Matérn(5/2) GP posterior (from 03c_smooth_matern.R) to predict
# at unobserved grid locations via the GP conditional distribution.
#
# For each posterior draw l, at each grid point s0:
#   eta(s0)^(l) = mu^(l) + W * (eta^(l) - mu^(l)) + N(0, cond_sd)
#
# where W = gamma^T Sigma^{-1} are kriging weights computed from
# the posterior mean GP hyperparameters (sigma_gp, phi_gp, nugget).
#
# Run from project root: Rscript R/05_interpolation_maps.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Posterior predictive maps (Matern GP)\n")
cat("========================================\n")

# ---- 1. Load Matérn GP posterior ----
cat("Loading Matern GP posterior...\n")

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_results.rds")

loc  <- s1$loc
meta <- s1$station_meta
ns   <- nrow(loc)

# eta draws [n_draws x ns]
eta_psi_draws <- s2$psi.selected
eta_tau_draws <- s2$tau.selected
eta_phi_draws <- s2$phi.selected
n_draws <- nrow(eta_psi_draws)

# GP hyperparameter draws
mu_psi_draws <- s2$beta_psi
mu_tau_draws <- s2$beta_tau
mu_phi_draws <- s2$beta_phi

cat("  Stations:", ns, "| Draws:", n_draws, "\n")

# ---- 2. Create prediction grid inside Andalucía ----
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

# ---- 3. Compute kriging weights from posterior mean GP hyperparameters ----
cat("Computing Matern(5/2) kriging weights...\n")

matern52 <- function(d, sigma2, phi) {
  s5 <- sqrt(5) * d / phi
  sigma2 * (1 + s5 + s5^2 / 3) * exp(-s5)
}

# Posterior mean GP hyperparameters
sigma_gp <- c(mean(s2$sigma_gp_psi), mean(s2$sigma_gp_tau), mean(s2$sigma_gp_phi))
phi_gp   <- c(mean(s2$phi_gp_psi), mean(s2$phi_gp_tau), mean(s2$phi_gp_phi))
nugget   <- c(mean(s2$nugget_psi), mean(s2$nugget_tau), mean(s2$nugget_phi))

cat(sprintf("  psi: sigma=%.3f, phi=%.3f, nugget=%.4f\n", sigma_gp[1], phi_gp[1], nugget[1]))
cat(sprintf("  tau: sigma=%.3f, phi=%.3f, nugget=%.4f\n", sigma_gp[2], phi_gp[2], nugget[2]))
cat(sprintf("  phi: sigma=%.3f, phi=%.3f, nugget=%.4f\n", sigma_gp[3], phi_gp[3], nugget[3]))

# Station-station distance matrix [ns x ns]
dist_ss <- as.matrix(dist(loc))

# Station-grid distance matrix [n_pred x ns]
dist_pg <- as.matrix(
  dist(rbind(as.matrix(pred_pts), as.matrix(data.frame(lon = loc[, "lon"], lat = loc[, "lat"]))))
)[1:n_pred, (n_pred + 1):(n_pred + ns)]

# Build covariance, invert, compute kriging weights for each parameter
compute_kriging <- function(k) {
  sig2 <- sigma_gp[k]^2
  nug2 <- nugget[k]^2
  phi  <- phi_gp[k]

  # Station-station covariance [ns x ns]
  Sigma <- matern52(dist_ss, sig2, phi)
  diag(Sigma) <- diag(Sigma) + nug2
  Sigma_inv <- solve(Sigma)

  # Cross-covariance [n_pred x ns] (no nugget at prediction points)
  gamma_mat <- matern52(dist_pg, sig2, phi)

  # Kriging weights [n_pred x ns]
  W <- gamma_mat %*% Sigma_inv

  # Conditional variance [n_pred]
  cond_var <- sig2 - rowSums(W * gamma_mat)
  cond_var <- pmax(cond_var, 1e-10)

  list(W = W, cond_sd = sqrt(cond_var))
}

krig <- lapply(1:3, compute_kriging)

cat("  Kriging weights computed\n")
cat(sprintf("  Conditional SD range — psi: [%.4f, %.4f], tau: [%.4f, %.4f], phi: [%.4f, %.4f]\n",
            min(krig[[1]]$cond_sd), max(krig[[1]]$cond_sd),
            min(krig[[2]]$cond_sd), max(krig[[2]]$cond_sd),
            min(krig[[3]]$cond_sd), max(krig[[3]]$cond_sd)))

# ---- 4. Sample from GP conditional at prediction locations ----
cat("Sampling from posterior predictive at grid points...\n")

chunk_size <- 2000
n_chunks <- ceiling(n_pred / chunk_size)

mu_gev_mean    <- numeric(n_pred)
sigma_gev_mean <- numeric(n_pred)
xi_gev_mean    <- numeric(n_pred)
rl100_mean     <- numeric(n_pred)
rl100_sd       <- numeric(n_pred)

M <- 100
set.seed(42)

# Precompute residuals (eta - mu) for all draws
resid_psi <- sweep(eta_psi_draws, 1, mu_psi_draws, "-")  # [n_draws x ns]
resid_tau <- sweep(eta_tau_draws, 1, mu_tau_draws, "-")
resid_phi <- sweep(eta_phi_draws, 1, mu_phi_draws, "-")

for (ch in seq_len(n_chunks)) {
  idx_start <- (ch - 1) * chunk_size + 1
  idx_end   <- min(ch * chunk_size, n_pred)
  idx <- idx_start:idx_end
  nc <- length(idx)

  if (ch %% 5 == 1 || ch == n_chunks)
    cat(sprintf("  Chunk %d/%d (points %d-%d)...\n", ch, n_chunks, idx_start, idx_end))

  # For each parameter: eta_pred = mu + W * (eta - mu) + N(0, cond_sd)
  predict_param <- function(k, resid, mu_draws) {
    W_ch <- krig[[k]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[k]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)

    # GP conditional mean + noise
    pred <- sweep(t(W_ch %*% t(resid)), 1, mu_draws, "+")  # [n_draws x nc]
    pred <- pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
    pred
  }

  psi_ch <- predict_param(1, resid_psi, mu_psi_draws)
  tau_ch <- predict_param(2, resid_tau, mu_tau_draws)
  phi_ch <- predict_param(3, resid_phi, mu_phi_draws)

  # Back-transform to GEV parameters
  mu_gev_ch    <- exp(psi_ch)
  sigma_gev_ch <- exp(psi_ch + tau_ch)
  xi_gev_ch    <- g(phi_ch)

  # 100-year return level
  rl100_ch <- mu_gev_ch + sigma_gev_ch / xi_gev_ch *
    ((-log(1 - 1 / M))^(-xi_gev_ch) - 1)

  mu_gev_mean[idx]    <- colMeans(mu_gev_ch)
  sigma_gev_mean[idx] <- colMeans(sigma_gev_ch)
  xi_gev_mean[idx]    <- colMeans(xi_gev_ch)
  rl100_mean[idx]     <- colMeans(rl100_ch)
  rl100_sd[idx]       <- apply(rl100_ch, 2, sd)
}

pred_pts$mu_mean    <- mu_gev_mean
pred_pts$sigma_mean <- sigma_gev_mean
pred_pts$xi_mean    <- xi_gev_mean
pred_pts$rl100_mean <- rl100_mean
pred_pts$rl100_sd   <- rl100_sd

cat("  Prediction done\n")

# ---- 5. Station posterior summaries ----
cat("Computing station summaries...\n")

psi_stn_draws <- s2$psi.selected
tau_stn_draws <- s2$tau.selected
phi_stn_draws <- s2$phi.selected

mu_stn    <- colMeans(exp(psi_stn_draws))
sigma_stn <- colMeans(exp(psi_stn_draws + tau_stn_draws))
xi_stn    <- colMeans(g(phi_stn_draws))
rl100_stn_draws <- exp(psi_stn_draws) +
  exp(psi_stn_draws + tau_stn_draws) / g(phi_stn_draws) *
  ((-log(1 - 1 / M))^(-g(phi_stn_draws)) - 1)
rl100_stn <- colMeans(rl100_stn_draws)

cat("  Station RL100 range:", round(min(rl100_stn), 1), "-",
    round(max(rl100_stn), 1), "mm\n")

# ---- 6. Prepare plot data ----
cat("Preparing plot data...\n")

grid_df <- data.frame(
  lon = pred_pts$lon, lat = pred_pts$lat,
  mu_mean = pred_pts$mu_mean, sigma_mean = pred_pts$sigma_mean,
  xi_mean = pred_pts$xi_mean, rl100_mean = pred_pts$rl100_mean
)

stn_df <- data.frame(
  lon = loc[, "lon"], lat = loc[, "lat"],
  mu_mean = mu_stn, sigma_mean = sigma_stn,
  xi_mean = xi_stn, rl100_mean = rl100_stn
)

cat("  Grid:", n_pred, "tiles | Stations:", ns, "circles\n")

# ---- 7. Plot filled maps ----
cat("Generating maps...\n")

stn_sf <- st_as_sf(meta, coords = c("longitud", "latitud"), crs = 4326)

neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real",
                                          "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

tile_res <- pred_res

base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

make_panel <- function(grid_data, stn_data, fill_var, fill_label, title, option = "C") {
  ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = grid_data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    geom_point(data = stn_data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
               shape = 21, size = 1.5, colour = "grey70", stroke = 0.15) +
    scale_fill_viridis_c(option = option, name = fill_label) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

p1 <- make_panel(grid_df, stn_df, "mu_mean", "mm",
                 expression(hat(mu) ~ "(location)"), "C")
p2 <- make_panel(grid_df, stn_df, "sigma_mean", "mm",
                 expression(hat(sigma) ~ "(scale)"), "C")
p3 <- make_panel(grid_df, stn_df, "xi_mean", expression(xi),
                 expression(hat(xi) ~ "(shape)"), "D")
p4 <- make_panel(grid_df, stn_df, "rl100_mean", "mm",
                 "100-year return level", "B")

p_all <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Posterior predictive: extreme precipitation across Andalucia",
    subtitle = sprintf(
      "Matern(5/2) GP | %d stations, %d grid points, %d draws | %.3f deg tiles",
      ns, n_pred, n_draws, pred_res),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/interpolated_maps.pdf", p_all, width = 14, height = 10)
ggsave("figures/interpolated_maps.png", p_all, width = 14, height = 10, dpi = 200)

cat("Saved figures/interpolated_maps.pdf and .png\n")

# ---- 8. Version 2: grid-only maps with hollow station markers ----
cat("Generating v2 maps (hollow station markers)...\n")

make_panel_v2 <- function(grid_data, stn_data, fill_var, fill_label, title, option = "C") {
  ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = grid_data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    geom_point(data = stn_data, aes(x = lon, y = lat),
               shape = 21, size = 1.5, fill = NA, colour = "grey40", stroke = 0.3) +
    scale_fill_viridis_c(option = option, name = fill_label) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

q1 <- make_panel_v2(grid_df, stn_df, "mu_mean", "mm",
                    expression(hat(mu) ~ "(location)"), "C")
q2 <- make_panel_v2(grid_df, stn_df, "sigma_mean", "mm",
                    expression(hat(sigma) ~ "(scale)"), "C")
q3 <- make_panel_v2(grid_df, stn_df, "xi_mean", expression(xi),
                    expression(hat(xi) ~ "(shape)"), "D")
q4 <- make_panel_v2(grid_df, stn_df, "rl100_mean", "mm",
                    "100-year return level", "B")

p_v2 <- (q1 + q2) / (q3 + q4) +
  plot_annotation(
    title = "Posterior predictive: extreme precipitation across Andalucia",
    subtitle = sprintf(
      "Matern(5/2) GP | %d stations, %d grid points, %d draws | %.3f deg tiles",
      ns, n_pred, n_draws, pred_res),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/interpolated_maps_v2.pdf", p_v2, width = 14, height = 10)
ggsave("figures/interpolated_maps_v2.png", p_v2, width = 14, height = 10, dpi = 200)

cat("Saved figures/interpolated_maps_v2.pdf and .png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
