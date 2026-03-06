# 06_station_diagnostics.R — Return level curves for select stations
#
# Classic GEV return level plot:
#   x-axis: return period (log scale)
#   y-axis: return level (mm)
#   Shows: empirical points, MLE fit, smoothed posterior fit + 90% credible band
#
# Run from project root: Rscript R/06_station_diagnostics.R

library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Station return level curves\n")
cat("========================================\n")

# ---- 1. Load data ----
cat("Loading data...\n")

am <- readRDS("data/annual_maxima_andalucia.rds")
s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

meta <- s1$station_meta
mles <- s1$mles
n_draws <- nrow(s2$psi.selected)

# ---- 2. Select stations ----
# Pick 6 stations spanning provinces and risk levels
selected <- c(
  "6155A",  # Málaga Aeropuerto — 66 yr, must include
  "5783",   # Sevilla Aeropuerto — 66 yr, had 241.6mm event
  "6058I",  # Estepona — wet coast, 211mm event
  "6325O",  # Almería Aeropuerto — dry SE
  "5402",   # Córdoba Aeropuerto — interior
  "5514"    # Granada Base Aérea — interior mountains
)

station_labels <- c(
  "6155A" = "M\u00e1laga Aeropuerto",
  "5783"  = "Sevilla Aeropuerto",
  "6058I" = "Estepona",
  "6325O" = "Almer\u00eda Aeropuerto",
  "5402"  = "C\u00f3rdoba Aeropuerto",
  "5514"  = "Granada Base A\u00e9rea"
)

cat("Selected stations:\n")
for (id in selected) {
  idx <- which(meta$indicativo == id)
  n_yr <- sum(am$indicativo == id)
  cat(sprintf("  %s: %s (%s) — %d years\n",
    id, station_labels[id], meta$provincia[idx], n_yr))
}

# ---- 3. GEV return level function ----
# Given (psi, tau, phi) in model parameterisation, compute return level for period T
rl_from_params <- function(psi, tau, phi, T) {
  mu    <- exp(psi)
  sigma <- exp(psi + tau)
  xi    <- g(phi)
  mu + sigma / xi * ((-log(1 - 1 / T))^(-xi) - 1)
}

# ---- 4. Compute curves ----
cat("Computing return level curves...\n")

return_periods <- exp(seq(log(1.001), log(200), length.out = 150))

plots <- list()

for (i in seq_along(selected)) {
  id <- selected[i]
  idx <- which(meta$indicativo == id)
  label <- station_labels[id]

  # Observed annual maxima
  obs <- am[am$indicativo == id, ]
  obs <- obs[order(obs$max_prec), ]
  n <- nrow(obs)
  # Empirical return periods (Gringorten plotting position)
  obs$rp <- (n + 0.12) / (seq_len(n) - 0.44 + 0.12)
  # Flip so largest has longest return period
  obs$rp <- rev(obs$rp)

  # MLE fit curve
  psi_mle <- mles[idx, "psi"]
  tau_mle <- mles[idx, "tau"]
  phi_mle <- mles[idx, "phi"]

  rl_mle <- sapply(return_periods, function(T)
    rl_from_params(psi_mle, tau_mle, phi_mle, T))

  # Smoothed posterior curves
  psi_draws <- s2$psi.selected[, idx]
  tau_draws <- s2$tau.selected[, idx]
  phi_draws <- s2$phi.selected[, idx]

  rl_draws <- matrix(NA, n_draws, length(return_periods))
  for (j in seq_along(return_periods)) {
    rl_draws[, j] <- rl_from_params(psi_draws, tau_draws, phi_draws, return_periods[j])
  }

  rl_post_mean <- colMeans(rl_draws)
  rl_post_lo   <- apply(rl_draws, 2, quantile, 0.05)
  rl_post_hi   <- apply(rl_draws, 2, quantile, 0.95)

  # Data frames for plotting
  curve_df <- data.frame(
    rp = return_periods,
    mle = rl_mle,
    post_mean = rl_post_mean,
    post_lo = rl_post_lo,
    post_hi = rl_post_hi
  )

  # Province and years info
  prov <- meta$provincia[idx]
  yr_range <- paste0(min(obs$year), "–", max(obs$year))

  plots[[i]] <- ggplot() +
    # 90% credible band
    geom_ribbon(data = curve_df, aes(x = rp, ymin = post_lo, ymax = post_hi),
                fill = "#3366CC", alpha = 0.2) +
    # MLE curve
    geom_line(data = curve_df, aes(x = rp, y = mle, colour = "MLE"),
              linewidth = 0.7, linetype = "dashed") +
    # Posterior mean curve
    geom_line(data = curve_df, aes(x = rp, y = post_mean, colour = "Smoothed"),
              linewidth = 0.9) +
    # Observed data points
    geom_point(data = obs, aes(x = rp, y = max_prec),
               shape = 21, fill = "grey20", colour = "white", size = 2, stroke = 0.4) +
    scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                  labels = c("2", "5", "10", "20", "50", "100", "200")) +
    scale_colour_manual(values = c("MLE" = "#CC3333", "Smoothed" = "#3366CC"),
                        name = NULL) +
    labs(
      title = label,
      subtitle = paste0(prov, " · ", n, " years (", yr_range, ")"),
      x = "Return period (years)",
      y = "Daily rainfall (mm)"
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

cat("  Curves computed for", length(plots), "stations\n")

# ---- 5. Assemble and save ----
cat("Assembling panel plot...\n")

p_all <- (plots[[1]] + plots[[2]] + plots[[3]]) /
         (plots[[4]] + plots[[5]] + plots[[6]]) +
  plot_annotation(
    title = "Return level curves: observed data vs model fits",
    subtitle = "Points = annual maxima | Dashed red = site-only MLE | Blue = spatially smoothed posterior (90% CI band)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("figures/station_return_level_curves.png",
       p_all, width = 16, height = 10, dpi = 200, bg = "white")

cat("Saved figures/station_return_level_curves.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
