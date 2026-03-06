# 06c_return_level_1x3.R — 1×4 panel: return level posterior means
#
# Simple, clean figure: 10, 20, 50, 100 year return levels side by side.
# Run from project root: Rscript R/06c_return_level_1x3.R

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

mu_psi_draws <- s2$beta_psi
mu_tau_draws <- s2$beta_tau
mu_phi_draws <- s2$beta_phi

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

cat("  Grid:", n_pred, "points\n")

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

# ---- 4. Return levels at T = 10, 20, 50, 100 ----
cat("Computing return levels...\n")

return_periods <- c(10, 20, 50, 100)
n_rp <- length(return_periods)

rl_mean <- matrix(0, n_pred, n_rp)

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
    cat(sprintf("  Chunk %d/%d...\n", ch, n_chunks))

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

  mu_gev    <- exp(psi_ch)
  sigma_gev <- exp(psi_ch + tau_ch)
  xi_gev    <- g(phi_ch)

  for (r in seq_len(n_rp)) {
    Tr <- return_periods[r]
    rl_ch <- mu_gev + sigma_gev / xi_gev * ((-log(1 - 1 / Tr))^(-xi_gev) - 1)
    rl_mean[idx, r] <- colMeans(rl_ch)
  }
}

# ---- 5. Build data frames ----
grid_list <- list()
for (r in seq_len(n_rp)) {
  grid_list[[r]] <- data.frame(
    lon = pred_pts$lon, lat = pred_pts$lat,
    rl_mean = rl_mean[, r],
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

# ---- 6. Plot 1x3 ----
cat("Generating 1x4 panel...\n")

neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

rl_min <- min(grid_df$rl_mean)
rl_max <- max(grid_df$rl_mean)

# AEMET alarm thresholds (precip in 12h, typical for Andalucía)
alarm_thresholds <- c(80, 120)        # orange, red (mm)
alarm_colours    <- c("#FF8C00", "red") # orange, red

legend_breaks <- seq(40, 200, by = 20)

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

  # Add AEMET alarm contour lines
  mid_lon <- mean(range(g_sub$lon))
  mid_lat <- mean(range(g_sub$lat))
  for (i in seq_along(alarm_thresholds)) {
    thr <- alarm_thresholds[i]
    if (thr > data_range[1] && thr < data_range[2]) {
      p <- p +
        geom_contour(data = g_sub, aes(x = lon, y = lat, z = rl_mean),
                     breaks = thr,
                     colour = "white", linewidth = 0.3, alpha = 0.85)

      # One label near map center
      diffs <- abs(g_sub$rl_mean - thr)
      candidates <- which(diffs < quantile(diffs, 0.01))
      if (length(candidates) == 0) candidates <- which.min(diffs)
      dist_center <- (g_sub$lon[candidates] - mid_lon)^2 +
                     (g_sub$lat[candidates] - mid_lat)^2
      best <- candidates[which.min(dist_center)]
      lbl <- data.frame(lon = g_sub$lon[best], lat = g_sub$lat[best],
                        label = thr)
      p <- p +
        geom_label(data = lbl, aes(x = lon, y = lat, label = label),
                   size = 1.8, fill = alarm_colours[i], colour = "white",
                   fontface = "bold", label.padding = unit(0.15, "lines"),
                   label.r = unit(0.1, "lines"))
    }
  }

  p +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    geom_point(data = s_sub, aes(x = lon, y = lat),
               shape = 21, size = 1.0, fill = NA, colour = "grey50", stroke = 0.3) +
    labs(title = paste0(rp_val, "-year return level")) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.key.height = unit(1.0, "cm"))
}

p_all <- make_panel(10) + make_panel(20) + make_panel(50) + make_panel(100) +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = "Maximum daily rainfall return levels across Andaluc\u00eda",
    subtitle = sprintf(
      "Posterior mean | Mat\u00e9rn(5/2) GP + PC priors | %d stations",
      ns),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/presentation/return_levels_2x2.png",
       p_all, width = 14, height = 12, dpi = 200, bg = "white")

cat("Saved figures/presentation/return_levels_2x2.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
