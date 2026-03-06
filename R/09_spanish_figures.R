# 09_spanish_figures.R — Spanish versions of key figures for Jesús
#
# Produces _es.png copies of:
#   - niveles_retorno.png
#   - probabilidad_superacion.png
#   - maximos_anuales_estaciones.png
#
# Run from project root: Rscript R/09_spanish_figures.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)

source("vendor/max_and_smooth/stage1_functions.R")

Sys.setlocale("LC_ALL", "en_US.UTF-8")

cat("========================================\n")
cat("Figuras en espa\u00f1ol\n")
cat("========================================\n")

# ---- Common data ----
cat("Cargando datos...\n")

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

loc  <- s1$loc
ns   <- nrow(loc)
meta <- s1$station_meta

eta_psi_draws <- s2$psi.selected
eta_tau_draws <- s2$tau.selected
eta_phi_draws <- s2$phi.selected
n_draws <- nrow(eta_psi_draws)

mu_psi_draws <- s2$beta_psi
mu_tau_draws <- s2$beta_tau
mu_phi_draws <- s2$beta_phi

# ---- Common geography ----
states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)
neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

# ---- Prediction grid ----
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

# ---- Kriging weights ----
cat("Calculando pesos de kriging...\n")

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

resid_psi <- sweep(eta_psi_draws, 1, mu_psi_draws, "-")
resid_tau <- sweep(eta_tau_draws, 1, mu_tau_draws, "-")
resid_phi <- sweep(eta_phi_draws, 1, mu_phi_draws, "-")

predict_chunk <- function(idx) {
  nc <- length(idx)
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
  list(mu = exp(psi_ch), sigma = exp(psi_ch + tau_ch), xi = g(phi_ch))
}

# ============================================================
# FIGURE 1: Return levels 2x2  (ES)
# ============================================================
cat("\n--- Figura 1: Niveles de retorno 2x2 ---\n")

return_periods <- c(10, 20, 50, 100)
n_rp <- length(return_periods)
rl_mean <- matrix(0, n_pred, n_rp)

chunk_size <- 2000
n_chunks <- ceiling(n_pred / chunk_size)
set.seed(42)

for (ch in seq_len(n_chunks)) {
  idx <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_pred)
  gev <- predict_chunk(idx)
  for (r in seq_len(n_rp)) {
    Tr <- return_periods[r]
    rl_ch <- gev$mu + gev$sigma / gev$xi * ((-log(1 - 1 / Tr))^(-gev$xi) - 1)
    rl_mean[idx, r] <- colMeans(rl_ch)
  }
}

grid_list <- list()
for (r in seq_len(n_rp)) {
  grid_list[[r]] <- data.frame(
    lon = pred_pts$lon, lat = pred_pts$lat,
    rl_mean = rl_mean[, r], rp = return_periods[r])
}
grid_df <- do.call(rbind, grid_list)

stn_list <- list()
for (r in seq_len(n_rp)) {
  Tr <- return_periods[r]
  mu_stn    <- exp(eta_psi_draws)
  sigma_stn <- exp(eta_psi_draws + eta_tau_draws)
  xi_stn    <- g(eta_phi_draws)
  rl_stn_draws <- mu_stn + sigma_stn / xi_stn * ((-log(1 - 1 / Tr))^(-xi_stn) - 1)
  stn_list[[r]] <- data.frame(
    lon = loc[, "lon"], lat = loc[, "lat"],
    rl_mean = colMeans(rl_stn_draws), rp = return_periods[r])
}
stn_df <- do.call(rbind, stn_list)

rl_min <- min(grid_df$rl_mean)
rl_max <- max(grid_df$rl_mean)

alarm_thresholds <- c(80, 120)
alarm_colours    <- c("#FF8C00", "red")
legend_breaks <- seq(40, 200, by = 20)

make_rl_panel <- function(rp_val) {
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

  mid_lon <- mean(range(g_sub$lon))
  mid_lat <- mean(range(g_sub$lat))
  for (i in seq_along(alarm_thresholds)) {
    thr <- alarm_thresholds[i]
    if (thr > data_range[1] && thr < data_range[2]) {
      p <- p +
        geom_contour(data = g_sub, aes(x = lon, y = lat, z = rl_mean),
                     breaks = thr, colour = "white", linewidth = 0.3, alpha = 0.85)
      diffs <- abs(g_sub$rl_mean - thr)
      candidates <- which(diffs < quantile(diffs, 0.01))
      if (length(candidates) == 0) candidates <- which.min(diffs)
      dist_center <- (g_sub$lon[candidates] - mid_lon)^2 +
                     (g_sub$lat[candidates] - mid_lat)^2
      best <- candidates[which.min(dist_center)]
      lbl <- data.frame(lon = g_sub$lon[best], lat = g_sub$lat[best], label = thr)
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
    labs(title = paste0("Nivel de retorno: ", rp_val, " a\u00f1os")) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.key.height = unit(1.0, "cm"))
}

p_rl <- make_rl_panel(10) + make_rl_panel(20) + make_rl_panel(50) + make_rl_panel(100) +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = "Niveles de retorno de precipitaci\u00f3n m\u00e1xima diaria en Andaluc\u00eda",
    subtitle = sprintf("Media posterior | %d estaciones", ns),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/jesus-figures/niveles_retorno.png",
       p_rl, width = 14, height = 12, dpi = 200, bg = "white")
cat("  Guardado niveles_retorno.png\n")

# ============================================================
# FIGURE 2: Exceedance probability  (ES)
# ============================================================
cat("\n--- Figura 2: Probabilidad de superaci\u00f3n ---\n")

thresholds <- c(100, 150, 200)
horizons   <- c(20, 50, 100)
n_thr <- length(thresholds)
n_hor <- length(horizons)

exc_mean <- array(0, dim = c(n_pred, n_thr, n_hor))
set.seed(42)

for (ch in seq_len(n_chunks)) {
  idx <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_pred)
  if (ch %% 5 == 1 || ch == n_chunks) cat(sprintf("  Bloque %d/%d...\n", ch, n_chunks))
  gev <- predict_chunk(idx)
  for (t in seq_len(n_thr)) {
    x <- thresholds[t]
    z <- 1 + gev$xi * (x - gev$mu) / gev$sigma
    z <- pmax(z, 1e-10)
    p_annual <- 1 - exp(-z^(-1 / gev$xi))
    for (h in seq_len(n_hor)) {
      p_horizon <- 1 - (1 - p_annual)^horizons[h]
      exc_mean[idx, t, h] <- colMeans(p_horizon)
    }
  }
}

grid_rows <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      lon = pred_pts$lon, lat = pred_pts$lat,
      prob_mean = exc_mean[, t, h],
      threshold = thresholds[t], horizon = horizons[h])
  }
}
exc_df <- do.call(rbind, grid_rows)

base_theme <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        legend.key.height = unit(0.7, "cm"))

make_exc_panel <- function(data, title) {
  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = prob_mean),
              width = pred_res, height = pred_res) +
    scale_fill_viridis_c(option = "B", name = "prob",
                         limits = c(0, 1), oob = scales::squish)

  # Contour lines at 0.25, 0.5, 0.75
  data_range <- range(data$prob_mean, na.rm = TRUE)
  contour_breaks <- c(0.25, 0.5, 0.75)
  valid_breaks <- contour_breaks[contour_breaks > data_range[1] &
                                 contour_breaks < data_range[2]]
  if (length(valid_breaks) > 0) {
    p <- p +
      geom_contour(data = data, aes(x = lon, y = lat, z = prob_mean),
                   breaks = valid_breaks, colour = "white", linewidth = 0.3, alpha = 0.85)
    mid_lon <- mean(range(data$lon))
    mid_lat <- mean(range(data$lat))
    label_df <- do.call(rbind, lapply(valid_breaks, function(lev) {
      diffs <- abs(data$prob_mean - lev)
      candidates <- which(diffs < quantile(diffs, 0.01))
      if (length(candidates) == 0) candidates <- which.min(diffs)
      dist_center <- (data$lon[candidates] - mid_lon)^2 +
                     (data$lat[candidates] - mid_lat)^2
      best <- candidates[which.min(dist_center)]
      data.frame(lon = data$lon[best], lat = data$lat[best], label = lev)
    }))
    p <- p +
      geom_text(data = label_df, aes(x = lon, y = lat, label = label),
                size = 1.8, colour = "white", fontface = "bold")
  }

  p + geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

panels_exc <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    g_sub <- exc_df[exc_df$threshold == thresholds[t] &
                    exc_df$horizon == horizons[h], ]
    title <- paste0("P(m\u00e1x > ", thresholds[t], " mm en ", horizons[h], " a\u00f1os)")
    panels_exc[[length(panels_exc) + 1]] <- make_exc_panel(g_sub, title)
  }
}

p_exc <- (panels_exc[[1]] + panels_exc[[2]] + panels_exc[[3]]) /
         (panels_exc[[4]] + panels_exc[[5]] + panels_exc[[6]]) /
         (panels_exc[[7]] + panels_exc[[8]] + panels_exc[[9]]) +
  plot_annotation(
    title = "Probabilidad de superaci\u00f3n: media posterior",
    subtitle = sprintf("%d estaciones | %d puntos de malla", ns, n_pred),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/jesus-figures/probabilidad_superacion.png",
       p_exc, width = 16, height = 13, dpi = 200, bg = "white")
cat("  Guardado probabilidad_superacion.png\n")

# ============================================================
# FIGURE 3: Station annual maxima  (ES)
# ============================================================
cat("\n--- Figura 3: M\u00e1ximos anuales por estaci\u00f3n ---\n")

am <- readRDS("data/annual_maxima_andalucia.rds")

rl_from_params <- function(psi, tau, phi, T) {
  mu    <- exp(psi)
  sigma <- exp(psi + tau)
  xi    <- g(phi)
  mu + sigma / xi * ((-log(1 - 1 / T))^(-xi) - 1)
}

selected <- c("6155A", "5783", "6058I", "6325O", "5402", "5514")
station_labels <- c(
  "6155A" = "M\u00e1laga Aeropuerto",
  "5783"  = "Sevilla Aeropuerto",
  "6058I" = "Estepona",
  "6325O" = "Almer\u00eda Aeropuerto",
  "5402"  = "C\u00f3rdoba Aeropuerto",
  "5514"  = "Granada Base A\u00e9rea"
)

rp_vals <- c(20, 50, 100)
rl_colours <- c("20-year" = "#2196F3", "50-year" = "#FF9800", "100-year" = "#E91E63")

plots_es <- list()
for (i in seq_along(selected)) {
  id <- selected[i]
  idx <- which(meta$indicativo == id)
  label <- station_labels[id]
  prov <- meta$provincia[idx]

  obs <- am[am$indicativo == id, ]
  obs <- obs[order(obs$year), ]

  psi_draws <- s2$psi.selected[, idx]
  tau_draws <- s2$tau.selected[, idx]
  phi_draws <- s2$phi.selected[, idx]

  rl_vals <- sapply(rp_vals, function(T) mean(rl_from_params(psi_draws, tau_draws, phi_draws, T)))
  names(rl_vals) <- paste0(rp_vals, "-year")

  obs$colour_cat <- "Normal"
  obs$colour_cat[obs$max_prec >= rl_vals["20-year"]] <- ">NR20"
  obs$colour_cat[obs$max_prec >= rl_vals["50-year"]] <- ">NR50"
  obs$colour_cat[obs$max_prec >= rl_vals["100-year"]] <- ">NR100"
  obs$colour_cat <- factor(obs$colour_cat,
    levels = c("Normal", ">NR20", ">NR50", ">NR100"))

  bar_colours <- c("Normal" = "grey60", ">NR20" = "#2196F3",
                    ">NR50" = "#FF9800", ">NR100" = "#E91E63")

  yr_range <- paste0(min(obs$year), "\u2013", max(obs$year))

  plots_es[[i]] <- ggplot(obs, aes(x = year, y = max_prec)) +
    geom_col(aes(fill = colour_cat), width = 0.8) +
    scale_fill_manual(values = bar_colours, name = NULL, drop = FALSE) +
    geom_hline(yintercept = rl_vals["20-year"],
               colour = "#2196F3", linewidth = 0.6, linetype = "dashed") +
    geom_hline(yintercept = rl_vals["50-year"],
               colour = "#FF9800", linewidth = 0.6, linetype = "dashed") +
    geom_hline(yintercept = rl_vals["100-year"],
               colour = "#E91E63", linewidth = 0.6, linetype = "dashed") +
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["20-year"],
             label = "NR20", colour = "#2196F3", size = 2.8, hjust = 0, fontface = "bold") +
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["50-year"],
             label = "NR50", colour = "#FF9800", size = 2.8, hjust = 0, fontface = "bold") +
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["100-year"],
             label = "NR100", colour = "#E91E63", size = 2.8, hjust = 0, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = label,
      subtitle = paste0(prov, " \u00b7 ", nrow(obs), " a\u00f1os (", yr_range, ")"),
      x = NULL,
      y = "Precip. m\u00e1xima diaria anual (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

p_stn <- (plots_es[[1]] + plots_es[[2]] + plots_es[[3]]) /
         (plots_es[[4]] + plots_es[[5]] + plots_es[[6]]) +
  plot_annotation(
    title = "Precipitaci\u00f3n m\u00e1xima diaria anual con niveles de retorno",
    subtitle = paste0(
      "Barras coloreadas por superaci\u00f3n: gris = normal, ",
      "azul = >NR20, naranja = >NR50, rojo = >NR100"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/jesus-figures/maximos_anuales_estaciones.png",
       p_stn, width = 18, height = 10, dpi = 200, bg = "white")
cat("  Guardado maximos_anuales_estaciones.png\n")

cat("\n========================================\n")
cat("Todas las figuras en espa\u00f1ol generadas.\n")
cat("========================================\n")
