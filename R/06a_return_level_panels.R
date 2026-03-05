# 06a_return_level_panels.R — 2x3 panel: return level surfaces + posterior SD
#
# Columns: T = 20, 50, 100 year return periods
# Top row: posterior mean return level (mm)
# Bottom row: posterior SD (mm)
#
# Uses GP kriging from PC prior model results.
# Run from project root: Rscript R/06a_return_level_panels.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Plot 1: Return level panel (2x3)\n")
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

# ---- 4. Posterior predictive: return levels at T = 20, 50, 100 ----
cat("Computing return levels at grid points...\n")

return_periods <- c(20, 50, 100)
n_rp <- length(return_periods)

# Allocators: n_pred x n_rp matrices
rl_mean <- matrix(0, n_pred, n_rp)
rl_sd   <- matrix(0, n_pred, n_rp)

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

  mu_gev    <- exp(psi_ch)
  sigma_gev <- exp(psi_ch + tau_ch)
  xi_gev    <- g(phi_ch)

  for (r in seq_len(n_rp)) {
    M <- return_periods[r]
    rl_ch <- mu_gev + sigma_gev / xi_gev * ((-log(1 - 1 / M))^(-xi_gev) - 1)
    rl_mean[idx, r] <- colMeans(rl_ch)
    rl_sd[idx, r]   <- apply(rl_ch, 2, sd)
  }
}

cat("  Return levels computed\n")

# ---- 5. Assemble data frames ----
grid_list <- list()
for (r in seq_len(n_rp)) {
  grid_list[[r]] <- data.frame(
    lon = pred_pts$lon, lat = pred_pts$lat,
    rl_mean = rl_mean[, r], rl_sd = rl_sd[, r],
    rp = return_periods[r]
  )
}
grid_df <- do.call(rbind, grid_list)
grid_df$rp_label <- paste0(grid_df$rp, "-year return level")
grid_df$rp_label <- factor(grid_df$rp_label,
  levels = paste0(return_periods, "-year return level"))

# Station summaries for overlay dots
stn_list <- list()
for (r in seq_len(n_rp)) {
  M <- return_periods[r]
  mu_stn    <- exp(eta_psi_draws)
  sigma_stn <- exp(eta_psi_draws + eta_tau_draws)
  xi_stn    <- g(eta_phi_draws)
  rl_stn_draws <- mu_stn + sigma_stn / xi_stn * ((-log(1 - 1 / M))^(-xi_stn) - 1)
  stn_list[[r]] <- data.frame(
    lon = loc[, "lon"], lat = loc[, "lat"],
    rl_mean = colMeans(rl_stn_draws),
    rl_sd = apply(rl_stn_draws, 2, sd),
    rp = return_periods[r]
  )
}
stn_df <- do.call(rbind, stn_list)
stn_df$rp_label <- paste0(stn_df$rp, "-year return level")
stn_df$rp_label <- factor(stn_df$rp_label,
  levels = paste0(return_periods, "-year return level"))

# Print ranges for diagnostics
for (r in seq_len(n_rp)) {
  rr <- grid_df[grid_df$rp == return_periods[r], ]
  cat(sprintf("  RL%d mean: %.1f - %.1f mm | SD: %.1f - %.1f mm\n",
    return_periods[r],
    min(rr$rl_mean), max(rr$rl_mean),
    min(rr$rl_sd), max(rr$rl_sd)))
}

# ---- 6. Plot ----
cat("Generating 2x3 panel plot...\n")

neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

tile_res <- pred_res

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11),
        legend.key.height = unit(0.8, "cm"))

make_panel <- function(data, fill_var, fill_label, title,
                       stn_data = NULL, option = "C",
                       limits = NULL) {
  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    scale_fill_viridis_c(option = option, name = fill_label, limits = limits) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4)

  if (!is.null(stn_data)) {
    p <- p + geom_point(data = stn_data, aes(x = lon, y = lat),
                         shape = 21, size = 0.8, fill = NA, colour = "grey50", stroke = 0.2)
  }

  p + labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

# Shared colour limits across columns for comparability
rl_max <- max(grid_df$rl_mean) * 1.02
sd_max <- max(grid_df$rl_sd) * 1.02

panels <- list()
for (r in seq_len(n_rp)) {
  M <- return_periods[r]
  g_sub <- grid_df[grid_df$rp == M, ]
  s_sub <- stn_df[stn_df$rp == M, ]

  panels[[r]] <- make_panel(g_sub, "rl_mean", "mm",
    paste0(M, "-year return level"),
    stn_data = s_sub, option = "C",
    limits = c(0, rl_max))

  panels[[n_rp + r]] <- make_panel(g_sub, "rl_sd", "mm",
    paste0(M, "-year SD"),
    stn_data = s_sub, option = "A",
    limits = c(0, sd_max))
}

p_all <- (panels[[1]] + panels[[2]] + panels[[3]]) /
         (panels[[4]] + panels[[5]] + panels[[6]]) +
  plot_annotation(
    title = "Return level estimates and posterior uncertainty across Andaluc\u00eda",
    subtitle = sprintf(
      "Mat\u00e9rn(5/2) GP + PC priors | %d stations, %d grid points, %d posterior draws",
      ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/return_level_panels.pdf", p_all, width = 16, height = 9, bg = "white")
ggsave("figures/return_level_panels.png", p_all, width = 16, height = 9, dpi = 200, bg = "white")

cat("Saved figures/return_level_panels.pdf and .png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
