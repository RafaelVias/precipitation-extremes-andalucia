# 07c_seasonal_rainfall.R — Seasonal cycle of daily rainfall at select stations
#
# Shows monthly patterns:
#   Top row: mean daily precipitation by month (includes dry days)
#   Bottom row: fraction of days exceeding 20mm (extreme day frequency)
#
# Run from project root: Rscript R/07c_seasonal_rainfall.R

library(dplyr)
library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

cat("========================================\n")
cat("Seasonal rainfall patterns\n")
cat("========================================\n")

# ---- 1. Load data ----
cat("Loading daily data...\n")

daily <- readRDS("data/daily_precip_andalucia_raw.rds")
daily$prec_num <- as.numeric(gsub(",", ".", daily$prec))
daily$month <- as.integer(format(daily$fecha, "%m"))
daily$year  <- as.integer(format(daily$fecha, "%Y"))

# ---- 2. Select stations ----
selected <- c("6155A", "5783", "6058I", "6325O", "5402", "5514")

station_labels <- c(
  "6155A" = "M\u00e1laga Aeropuerto",
  "5783"  = "Sevilla Aeropuerto",
  "6058I" = "Estepona",
  "6325O" = "Almer\u00eda Aeropuerto",
  "5402"  = "C\u00f3rdoba Aeropuerto",
  "5514"  = "Granada Base A\u00e9rea"
)

month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# ---- 3. Compute monthly statistics ----
cat("Computing monthly statistics...\n")

all_stats <- list()

for (id in selected) {
  sub <- daily %>%
    filter(indicativo == id, !is.na(prec_num))

  monthly <- sub %>%
    group_by(month) %>%
    summarise(
      mean_precip = mean(prec_num),
      p90         = quantile(prec_num, 0.90),
      p95         = quantile(prec_num, 0.95),
      p99         = quantile(prec_num, 0.99),
      frac_wet    = mean(prec_num > 1),     # P(rain > 1mm)
      frac_heavy  = mean(prec_num > 20),    # P(rain > 20mm)
      max_precip  = max(prec_num),
      n_days      = n(),
      .groups     = "drop"
    ) %>%
    mutate(station = station_labels[id])

  all_stats[[id]] <- monthly

  cat(sprintf("  %s: driest month = %s (%.1f mm/day), wettest = %s (%.1f mm/day)\n",
    station_labels[id],
    month.abb[monthly$month[which.min(monthly$mean_precip)]],
    min(monthly$mean_precip),
    month.abb[monthly$month[which.max(monthly$mean_precip)]],
    max(monthly$mean_precip)))
}

stats_df <- bind_rows(all_stats)
stats_df$station <- factor(stats_df$station, levels = station_labels)

# ---- 4. Also compute the month of annual maxima ----
cat("Computing month of annual maxima...\n")

am_month <- list()
for (id in selected) {
  sub <- daily %>%
    filter(indicativo == id, !is.na(prec_num)) %>%
    group_by(year) %>%
    slice_max(prec_num, n = 1, with_ties = FALSE) %>%
    ungroup()

  month_counts <- sub %>%
    count(month) %>%
    mutate(
      frac = n / sum(n),
      station = station_labels[id]
    )

  # Fill in missing months with 0
  all_months <- data.frame(month = 1:12)
  month_counts <- merge(all_months, month_counts, by = "month", all.x = TRUE)
  month_counts$n[is.na(month_counts$n)] <- 0
  month_counts$frac[is.na(month_counts$frac)] <- 0
  month_counts$station <- station_labels[id]

  am_month[[id]] <- month_counts
}

am_df <- bind_rows(am_month)
am_df$station <- factor(am_df$station, levels = station_labels)

# ---- 5. Plot ----
cat("Generating plots...\n")

bar_fill <- "#3366CC"
heavy_fill <- "#E91E63"
am_fill <- "#FF9800"

# Top row: mean daily rainfall by month
p_mean <- ggplot(stats_df, aes(x = month, y = mean_precip)) +
  geom_col(fill = bar_fill, alpha = 0.8, width = 0.7) +
  facet_wrap(~ station, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  labs(
    title = "Mean daily precipitation by month",
    x = NULL,
    y = "mm/day"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

# Middle row: fraction of heavy rain days (>20mm)
p_heavy <- ggplot(stats_df, aes(x = month, y = frac_heavy * 100)) +
  geom_col(fill = heavy_fill, alpha = 0.8, width = 0.7) +
  facet_wrap(~ station, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  labs(
    title = "Frequency of heavy rainfall days (> 20 mm/day)",
    x = NULL,
    y = "% of days"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

# Bottom row: month of annual maximum
p_am <- ggplot(am_df, aes(x = month, y = frac * 100)) +
  geom_col(fill = am_fill, alpha = 0.8, width = 0.7) +
  facet_wrap(~ station, nrow = 1, scales = "free_y") +
  scale_x_continuous(breaks = 1:12, labels = month_labels) +
  labs(
    title = "When do annual maxima occur? (month of wettest day each year)",
    x = NULL,
    y = "% of years"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 12)
  )

p_all <- p_mean / p_heavy / p_am +
  plot_annotation(
    title = "Seasonal rainfall patterns across Andaluc\u00eda",
    subtitle = "Daily precipitation data from AEMET (1960\u20132025)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, colour = "grey40")
    )
  )

ggsave("figures/jesus-figures/seasonal_rainfall.png",
       p_all, width = 18, height = 11, dpi = 200, bg = "white")

cat("Saved figures/jesus-figures/seasonal_rainfall.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
