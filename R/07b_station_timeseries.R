# 07b_station_timeseries.R — Annual maxima time series for select stations
#
# Bar chart of annual maxima per year with horizontal lines for
# smoothed posterior return levels (RL20, RL50, RL100).
#
# Run from project root: Rscript R/07b_station_timeseries.R

library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Station annual maxima time series\n")
cat("========================================\n")

# ---- 1. Load data ----
cat("Loading data...\n")

am <- readRDS("data/annual_maxima_andalucia.rds")
s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

meta <- s1$station_meta
n_draws <- nrow(s2$psi.selected)

# ---- 2. Select stations ----
selected <- c(
  "6155A",  # Málaga Aeropuerto
  "5783",   # Sevilla Aeropuerto
  "6058I",  # Estepona
  "6325O",  # Almería Aeropuerto
  "5402",   # Córdoba Aeropuerto
  "5514"    # Granada Base Aérea
)

station_labels <- c(
  "6155A" = "M\u00e1laga Aeropuerto",
  "5783"  = "Sevilla Aeropuerto",
  "6058I" = "Estepona",
  "6325O" = "Almer\u00eda Aeropuerto",
  "5402"  = "C\u00f3rdoba Aeropuerto",
  "5514"  = "Granada Base A\u00e9rea"
)

# ---- 3. Compute return levels from smoothed posterior ----
cat("Computing return levels...\n")

rl_from_params <- function(psi, tau, phi, T) {
  mu    <- exp(psi)
  sigma <- exp(psi + tau)
  xi    <- g(phi)
  mu + sigma / xi * ((-log(1 - 1 / T))^(-xi) - 1)
}

return_periods <- c(20, 50, 100)
rl_colours <- c("20-year" = "#2196F3", "50-year" = "#FF9800", "100-year" = "#E91E63")

plots <- list()

for (i in seq_along(selected)) {
  id <- selected[i]
  idx <- which(meta$indicativo == id)
  label <- station_labels[id]
  prov <- meta$provincia[idx]

  # Annual maxima for this station
  obs <- am[am$indicativo == id, ]
  obs <- obs[order(obs$year), ]

  # Smoothed return levels (posterior mean)
  psi_draws <- s2$psi.selected[, idx]
  tau_draws <- s2$tau.selected[, idx]
  phi_draws <- s2$phi.selected[, idx]

  rl_vals <- sapply(return_periods, function(T) {
    mean(rl_from_params(psi_draws, tau_draws, phi_draws, T))
  })
  names(rl_vals) <- paste0(return_periods, "-year")

  # Colour bars by whether they exceed RL thresholds
  obs$colour_cat <- "Normal"
  obs$colour_cat[obs$max_prec >= rl_vals["20-year"]] <- ">RL20"
  obs$colour_cat[obs$max_prec >= rl_vals["50-year"]] <- ">RL50"
  obs$colour_cat[obs$max_prec >= rl_vals["100-year"]] <- ">RL100"
  obs$colour_cat <- factor(obs$colour_cat,
    levels = c("Normal", ">RL20", ">RL50", ">RL100"))

  bar_colours <- c("Normal" = "grey60", ">RL20" = "#2196F3",
                    ">RL50" = "#FF9800", ">RL100" = "#E91E63")

  # RL line data
  rl_df <- data.frame(
    label = paste0(return_periods, "-year"),
    value = rl_vals,
    colour = rl_colours[paste0(return_periods, "-year")]
  )

  yr_range <- paste0(min(obs$year), "\u2013", max(obs$year))

  plots[[i]] <- ggplot(obs, aes(x = year, y = max_prec)) +
    geom_col(aes(fill = colour_cat), width = 0.8) +
    scale_fill_manual(values = bar_colours, name = NULL, drop = FALSE) +
    # RL lines
    geom_hline(yintercept = rl_vals["20-year"],
               colour = "#2196F3", linewidth = 0.6, linetype = "dashed") +
    geom_hline(yintercept = rl_vals["50-year"],
               colour = "#FF9800", linewidth = 0.6, linetype = "dashed") +
    geom_hline(yintercept = rl_vals["100-year"],
               colour = "#E91E63", linewidth = 0.6, linetype = "dashed") +
    # RL labels on right margin
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["20-year"],
             label = "RL20", colour = "#2196F3", size = 2.8, hjust = 0, fontface = "bold") +
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["50-year"],
             label = "RL50", colour = "#FF9800", size = 2.8, hjust = 0, fontface = "bold") +
    annotate("text", x = max(obs$year) + 1.5, y = rl_vals["100-year"],
             label = "RL100", colour = "#E91E63", size = 2.8, hjust = 0, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(
      title = label,
      subtitle = paste0(prov, " \u00b7 ", nrow(obs), " years (", yr_range, ")"),
      x = NULL,
      y = "Annual max daily rainfall (mm)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  cat(sprintf("  %s: RL20=%.0f, RL50=%.0f, RL100=%.0f mm\n",
    label, rl_vals[1], rl_vals[2], rl_vals[3]))
}

# ---- 4. Assemble and save ----
cat("Assembling panel plot...\n")

p_all <- (plots[[1]] + plots[[2]] + plots[[3]]) /
         (plots[[4]] + plots[[5]] + plots[[6]]) +
  plot_annotation(
    title = "Annual maximum daily rainfall with smoothed return level thresholds",
    subtitle = paste0(
      "Bars coloured by exceedance: grey = normal, ",
      "blue = >RL20, orange = >RL50, red = >RL100"),
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  )

ggsave("figures/jesus-figures/station_annual_maxima.png",
       p_all, width = 18, height = 10, dpi = 200, bg = "white")

cat("Saved figures/jesus-figures/station_annual_maxima.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
