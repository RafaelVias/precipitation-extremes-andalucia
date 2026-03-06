# 08_spanish_figures.R — Spanish versions of all figures
#
# Regenerates the 4 main figures (station map, return level maps,
# exceedance probability, station diagnostics) with Spanish labels.
# Output: figures/es/
#
# Run from project root: Rscript R/08_spanish_figures.R

library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(dplyr)

source("vendor/max_and_smooth/stage1_functions.R")

Sys.setlocale("LC_ALL", "es_ES.UTF-8")

cat("========================================\n")
cat("Figuras en espa\u00f1ol\n")
cat("========================================\n")

dir.create("figures/es", showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Shared geographic data
# =============================================================================

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)
neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

prov_labels <- data.frame(
  name = c("Huelva", "Sevilla", "C\u00f3rdoba", "Ja\u00e9n",
           "C\u00e1diz", "M\u00e1laga", "Granada", "Almer\u00eda"),
  lon  = c(-6.95, -5.85, -4.78, -3.75, -5.85, -4.55, -3.45, -2.40),
  lat  = c(37.65, 37.50, 37.95, 37.85, 36.55, 36.80, 37.25, 37.05)
)

# =============================================================================
# Figure 0: Station map
# =============================================================================
cat("\n--- Mapa de estaciones ---\n")

am <- readRDS("data/annual_maxima_andalucia.rds")
am <- am %>% filter(indicativo != "6381")

stn_summary <- am %>%
  group_by(indicativo) %>%
  summarise(
    n_years  = n(),
    max_prec = max(max_prec, na.rm = TRUE),
    nombre   = first(nombre),
    provincia = first(provincia),
    lon      = first(longitud),
    lat      = first(latitud),
    .groups  = "drop"
  )

p0 <- ggplot() +
  geom_sf(data = port_crop, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = mor_crop, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = neighbours, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = andalucia_provs, fill = "grey98", colour = "grey60", linewidth = 0.3) +
  geom_text(data = prov_labels, aes(x = lon, y = lat, label = name),
            size = 3.0, colour = "grey50", fontface = "italic") +
  geom_point(data = stn_summary,
             aes(x = lon, y = lat, fill = max_prec, size = n_years),
             shape = 21, colour = "grey30", stroke = 0.3) +
  scale_fill_viridis_c(option = "B", name = "M\u00e1ximo\nobservado (mm)",
                       breaks = c(50, 100, 150, 200, 250)) +
  scale_size_continuous(name = "Longitud del\nregistro (a\u00f1os)",
                        range = c(1.5, 5),
                        breaks = c(10, 20, 30, 40, 50, 60, 70)) +
  geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
  labs(
    title = "Red de estaciones AEMET en Andaluc\u00eda",
    subtitle = sprintf("%d estaciones | %d\u2013%d a\u00f1os de registro | CC: \u226590%% completitud diaria",
                       nrow(stn_summary), min(am$year), max(am$year))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave("figures/es/station_map.png", p0, width = 12, height = 7, dpi = 200, bg = "white")
cat("  Guardado figures/es/station_map.png\n")

# =============================================================================
# Shared model data for figures 1-3
# =============================================================================

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

loc  <- s1$loc
ns   <- nrow(loc)
meta <- s1$station_meta
mles <- s1$mles

eta_psi_draws <- s2$psi.selected
eta_tau_draws <- s2$tau.selected
eta_phi_draws <- s2$phi.selected
n_draws <- nrow(eta_psi_draws)

mu_psi_draws <- s2$beta_psi
mu_tau_draws <- s2$beta_tau
mu_phi_draws <- s2$beta_phi

# =============================================================================
# Shared prediction grid and kriging
# =============================================================================
cat("\nCalculando pesos de kriging...\n")

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

matern52 <- function(d, sigma2, phi) {
  s5 <- sqrt(5) * d / phi
  sigma2 * (1 + s5 + s5^2 / 3) * exp(-s5)
}

sigma_gp <- c(mean(s2$sigma_gp_psi), mean(s2$sigma_gp_tau), mean(s2$sigma_gp_phi))
phi_gp   <- c(mean(s2$phi_gp_psi), mean(s2$phi_gp_tau), mean(s2$phi_gp_phi))
nugget_v <- c(mean(s2$nugget_psi), mean(s2$nugget_tau), mean(s2$nugget_phi))

dist_ss <- as.matrix(dist(loc))
dist_pg <- as.matrix(
  dist(rbind(as.matrix(pred_pts), as.matrix(data.frame(lon = loc[, "lon"], lat = loc[, "lat"]))))
)[1:n_pred, (n_pred + 1):(n_pred + ns)]

compute_kriging <- function(k) {
  sig2 <- sigma_gp[k]^2
  nug2 <- nugget_v[k]^2
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

predict_at <- function(idx) {
  nc <- length(idx)
  predict_param <- function(k, resid, mu_draws) {
    W_ch <- krig[[k]]$W[idx, , drop = FALSE]
    sd_ch <- matrix(krig[[k]]$cond_sd[idx], nrow = n_draws, ncol = nc, byrow = TRUE)
    pred <- sweep(t(W_ch %*% t(resid)), 1, mu_draws, "+")
    pred + matrix(rnorm(n_draws * nc), n_draws, nc) * sd_ch
  }
  list(
    psi = predict_param(1, resid_psi, mu_psi_draws),
    tau = predict_param(2, resid_tau, mu_tau_draws),
    phi = predict_param(3, resid_phi, mu_phi_draws)
  )
}

# =============================================================================
# Figure 1: Return level maps
# =============================================================================
cat("\n--- Mapas de niveles de retorno ---\n")

return_periods_map <- c(10, 20, 50, 100)
n_rp <- length(return_periods_map)
rl_mean <- matrix(0, n_pred, n_rp)

chunk_size <- 2000
n_chunks <- ceiling(n_pred / chunk_size)
set.seed(42)

for (ch in seq_len(n_chunks)) {
  idx <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_pred)
  nc <- length(idx)
  pp <- predict_at(idx)
  mu_gev <- exp(pp$psi); sigma_gev <- exp(pp$psi + pp$tau); xi_gev <- g(pp$phi)
  for (r in seq_len(n_rp)) {
    Tr <- return_periods_map[r]
    rl_ch <- mu_gev + sigma_gev / xi_gev * ((-log(1 - 1 / Tr))^(-xi_gev) - 1)
    rl_mean[idx, r] <- colMeans(rl_ch)
  }
}

grid_list <- list()
for (r in seq_len(n_rp)) {
  grid_list[[r]] <- data.frame(lon = pred_pts$lon, lat = pred_pts$lat,
                                rl_mean = rl_mean[, r], rp = return_periods_map[r])
}
grid_df <- do.call(rbind, grid_list)

rl_min <- min(grid_df$rl_mean); rl_max <- max(grid_df$rl_mean)
alarm_thresholds <- c(80, 120)
alarm_colours <- c("#FF8C00", "red")
legend_breaks <- seq(40, 200, by = 20)

make_rl_panel <- function(rp_val) {
  g_sub <- grid_df[grid_df$rp == rp_val, ]
  data_range <- range(g_sub$rl_mean, na.rm = TRUE)

  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = g_sub, aes(x = lon, y = lat, fill = rl_mean),
              width = pred_res, height = pred_res) +
    scale_fill_viridis_c(option = "B", name = "mm",
                         limits = c(rl_min, rl_max), breaks = legend_breaks)

  mid_lon <- mean(range(g_sub$lon)); mid_lat <- mean(range(g_sub$lat))
  for (i in seq_along(alarm_thresholds)) {
    thr <- alarm_thresholds[i]
    if (thr > data_range[1] && thr < data_range[2]) {
      p <- p + geom_contour(data = g_sub, aes(x = lon, y = lat, z = rl_mean),
                             breaks = thr, colour = "white", linewidth = 0.3, alpha = 0.85)
      diffs <- abs(g_sub$rl_mean - thr)
      candidates <- which(diffs < quantile(diffs, 0.01))
      if (length(candidates) == 0) candidates <- which.min(diffs)
      dist_center <- (g_sub$lon[candidates] - mid_lon)^2 + (g_sub$lat[candidates] - mid_lat)^2
      best <- candidates[which.min(dist_center)]
      lbl <- data.frame(lon = g_sub$lon[best], lat = g_sub$lat[best], label = thr)
      p <- p + geom_label(data = lbl, aes(x = lon, y = lat, label = label),
                           size = 1.8, fill = alarm_colours[i], colour = "white",
                           fontface = "bold", label.padding = unit(0.15, "lines"),
                           label.r = unit(0.1, "lines"))
    }
  }

  p + geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = paste0("Nivel de retorno a ", rp_val, " a\u00f1os")) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          legend.key.height = unit(1.0, "cm"),
          plot.margin = margin(2, 2, 2, 2))
}

p_rl <- make_rl_panel(10) + make_rl_panel(20) + make_rl_panel(50) + make_rl_panel(100) +
  plot_layout(nrow = 2, ncol = 2) +
  plot_annotation(
    title = "Niveles de retorno de precipitaci\u00f3n m\u00e1xima diaria en Andaluc\u00eda",
    subtitle = sprintf("Media a posteriori | Mat\u00e9rn(5/2) GP + priors PC | %d estaciones", ns),
    theme = theme(plot.title = element_text(face = "bold", size = 15),
                  plot.subtitle = element_text(size = 11, colour = "grey40"),
                  plot.margin = margin(2, 2, 2, 2))
  )

ggsave("figures/es/return_level_maps.png", p_rl, width = 14, height = 10, dpi = 200, bg = "white")
cat("  Guardado figures/es/return_level_maps.png\n")

# =============================================================================
# Figure 2: Exceedance probability
# =============================================================================
cat("\n--- Probabilidad de excedencia ---\n")

thresholds_exc <- c(100, 150, 200)
horizons <- c(20, 50, 100)
n_thr <- length(thresholds_exc); n_hor <- length(horizons)
exc_mean <- array(0, dim = c(n_pred, n_thr, n_hor))

set.seed(42)
for (ch in seq_len(n_chunks)) {
  idx <- ((ch - 1) * chunk_size + 1):min(ch * chunk_size, n_pred)
  pp <- predict_at(idx)
  mu_gev <- exp(pp$psi); sigma_gev <- exp(pp$psi + pp$tau); xi_gev <- g(pp$phi)
  for (t in seq_len(n_thr)) {
    x <- thresholds_exc[t]
    z <- pmax(1 + xi_gev * (x - mu_gev) / sigma_gev, 1e-10)
    p_annual <- 1 - exp(-z^(-1 / xi_gev))
    for (h in seq_len(n_hor)) {
      exc_mean[idx, t, h] <- colMeans(1 - (1 - p_annual)^horizons[h])
    }
  }
}

grid_rows <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      lon = pred_pts$lon, lat = pred_pts$lat,
      prob_mean = exc_mean[, t, h],
      threshold = thresholds_exc[t], horizon = horizons[h]
    )
  }
}
exc_df <- do.call(rbind, grid_rows)

base_theme_es <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 10),
        legend.key.height = unit(0.7, "cm"))

make_exc_panel <- function(thr_val, hor_val) {
  g_sub <- exc_df[exc_df$threshold == thr_val & exc_df$horizon == hor_val, ]
  title <- paste0("P(m\u00e1x > ", thr_val, " mm en ", hor_val, " a\u00f1os)")

  p <- ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = g_sub, aes(x = lon, y = lat, fill = prob_mean),
              width = pred_res, height = pred_res) +
    scale_fill_viridis_c(option = "B", name = "prob",
                         limits = c(0, 1), oob = scales::squish)

  data_range <- range(g_sub$prob_mean, na.rm = TRUE)
  contour_breaks <- c(0.25, 0.5, 0.75)
  valid_breaks <- contour_breaks[contour_breaks > data_range[1] & contour_breaks < data_range[2]]
  if (length(valid_breaks) > 0) {
    p <- p + geom_contour(data = g_sub, aes(x = lon, y = lat, z = prob_mean),
                           breaks = valid_breaks, colour = "white", linewidth = 0.3, alpha = 0.85)
    mid_lon <- mean(range(g_sub$lon)); mid_lat <- mean(range(g_sub$lat))
    label_df <- do.call(rbind, lapply(valid_breaks, function(lev) {
      diffs <- abs(g_sub$prob_mean - lev)
      candidates <- which(diffs < quantile(diffs, 0.01))
      if (length(candidates) == 0) candidates <- which.min(diffs)
      dist_center <- (g_sub$lon[candidates] - mid_lon)^2 + (g_sub$lat[candidates] - mid_lat)^2
      best <- candidates[which.min(dist_center)]
      data.frame(lon = g_sub$lon[best], lat = g_sub$lat[best], label = lev)
    }))
    p <- p + geom_text(data = label_df, aes(x = lon, y = lat, label = label),
                        size = 1.8, colour = "white", fontface = "bold")
  }

  p + geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme_es
}

panels_es <- list()
for (t in seq_len(n_thr)) {
  for (h in seq_len(n_hor)) {
    panels_es[[length(panels_es) + 1]] <- make_exc_panel(thresholds_exc[t], horizons[h])
  }
}

p_exc <- (panels_es[[1]] + panels_es[[2]] + panels_es[[3]]) /
         (panels_es[[4]] + panels_es[[5]] + panels_es[[6]]) /
         (panels_es[[7]] + panels_es[[8]] + panels_es[[9]]) +
  plot_annotation(
    title = "Probabilidad de excedencia: media a posteriori",
    subtitle = sprintf("Mat\u00e9rn(5/2) GP + priors PC | %d estaciones, %d puntos de malla, %d muestras a posteriori",
                       ns, n_pred, n_draws),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/es/exceedance_prob.png", p_exc, width = 16, height = 13, dpi = 200, bg = "white")
cat("  Guardado figures/es/exceedance_prob.png\n")

# =============================================================================
# Figure 3: Station diagnostics
# =============================================================================
cat("\n--- Curvas de nivel de retorno por estaci\u00f3n ---\n")

selected <- c("6155A", "5783", "6058I", "6325O", "5402", "5514")
station_labels <- c(
  "6155A" = "M\u00e1laga Aeropuerto", "5783"  = "Sevilla Aeropuerto",
  "6058I" = "Estepona",              "6325O" = "Almer\u00eda Aeropuerto",
  "5402"  = "C\u00f3rdoba Aeropuerto", "5514"  = "Granada Base A\u00e9rea"
)

rl_from_params <- function(psi, tau, phi, T) {
  mu <- exp(psi); sigma <- exp(psi + tau); xi <- g(phi)
  mu + sigma / xi * ((-log(1 - 1 / T))^(-xi) - 1)
}

return_periods_diag <- exp(seq(log(1.001), log(200), length.out = 150))
plots_es <- list()

for (i in seq_along(selected)) {
  id <- selected[i]
  idx <- which(meta$indicativo == id)
  label <- station_labels[id]

  obs <- am[am$indicativo == id, ]
  obs <- obs[order(obs$max_prec), ]
  n <- nrow(obs)
  obs$rp <- rev((n + 0.12) / (seq_len(n) - 0.44 + 0.12))

  psi_mle <- mles[idx, "psi"]; tau_mle <- mles[idx, "tau"]; phi_mle <- mles[idx, "phi"]
  rl_mle <- sapply(return_periods_diag, function(T) rl_from_params(psi_mle, tau_mle, phi_mle, T))

  psi_d <- s2$psi.selected[, idx]; tau_d <- s2$tau.selected[, idx]; phi_d <- s2$phi.selected[, idx]
  rl_draws <- matrix(NA, n_draws, length(return_periods_diag))
  for (j in seq_along(return_periods_diag))
    rl_draws[, j] <- rl_from_params(psi_d, tau_d, phi_d, return_periods_diag[j])

  curve_df <- data.frame(
    rp = return_periods_diag, mle = rl_mle,
    post_mean = colMeans(rl_draws),
    post_lo = apply(rl_draws, 2, quantile, 0.05),
    post_hi = apply(rl_draws, 2, quantile, 0.95)
  )

  prov <- meta$provincia[idx]
  yr_range <- paste0(min(obs$year), "\u2013", max(obs$year))

  plots_es[[i]] <- ggplot() +
    geom_ribbon(data = curve_df, aes(x = rp, ymin = post_lo, ymax = post_hi),
                fill = "#3366CC", alpha = 0.2) +
    geom_line(data = curve_df, aes(x = rp, y = mle, colour = "EMV"),
              linewidth = 0.7, linetype = "dashed") +
    geom_line(data = curve_df, aes(x = rp, y = post_mean, colour = "Suavizado"),
              linewidth = 0.9) +
    geom_point(data = obs, aes(x = rp, y = max_prec),
               shape = 21, fill = "grey20", colour = "white", size = 2, stroke = 0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = c("2", "5", "10", "20", "50", "100", "200")) +
    scale_colour_manual(values = c("EMV" = "#CC3333", "Suavizado" = "#3366CC"), name = NULL) +
    labs(
      title = label,
      subtitle = paste0(prov, " \u00b7 ", n, " a\u00f1os (", yr_range, ")"),
      x = "Per\u00edodo de retorno (a\u00f1os)",
      y = "Precipitaci\u00f3n diaria (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank()
    )
}

p_diag <- (plots_es[[1]] + plots_es[[2]] + plots_es[[3]]) /
          (plots_es[[4]] + plots_es[[5]] + plots_es[[6]]) +
  plot_annotation(
    title = "Curvas de nivel de retorno: datos observados vs ajustes del modelo",
    subtitle = "Puntos = m\u00e1ximos anuales | Rojo discontinuo = EMV por estaci\u00f3n | Azul = posterior suavizado (banda IC 90%)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/es/station_return_level_curves.png",
       p_diag, width = 16, height = 10, dpi = 200, bg = "white")
cat("  Guardado figures/es/station_return_level_curves.png\n")

cat("\n========================================\n")
cat("Todas las figuras guardadas en figures/es/\n")
cat("========================================\n")
