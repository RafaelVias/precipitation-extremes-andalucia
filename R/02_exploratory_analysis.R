# =============================================================================
# 02_exploratory_analysis.R
# Exploratory data analysis of annual precipitation maxima in Andalucia
# =============================================================================
#
# Inputs:
#   - data/annual_maxima_andalucia.rds
#   - data/stations_andalucia.rds
#
# Outputs (figures/):
#   - station_map.pdf              Station locations coloured by record length
#   - annual_maxima_histogram.pdf  Distribution of annual maxima
#   - record_length_dist.pdf       Record length by station
#   - maxima_by_province.pdf       Boxplots by province
#   - maxima_vs_elevation.pdf      Annual maxima vs elevation
#   - maxima_vs_coords.pdf         Annual maxima vs lat/lon
#   - empirical_variogram.pdf      Empirical variogram of median annual max
#   - gev_fits_selected.pdf        At-site GEV fits at long-record stations
# =============================================================================

library(dplyr)
library(ggplot2)
library(sf)
library(gstat)
library(evd)

# ---- Load data -------------------------------------------------------------
annual_max <- readRDS("data/annual_maxima_andalucia.rds")
stations   <- readRDS("data/stations_andalucia.rds")

cat("Annual maxima:", nrow(annual_max), "obs from",
    n_distinct(annual_max$indicativo), "stations\n")

# Station-level summaries
station_summary <- annual_max %>%
  group_by(indicativo, nombre, provincia, altitud, latitud, longitud) %>%
  summarise(
    n_years    = n(),
    median_max = median(max_prec),
    mean_max   = mean(max_prec),
    max_max    = max(max_prec),
    .groups    = "drop"
  )

# ---- Figure 1: Station map -------------------------------------------------
cat("Figure 1: Station map...\n")

# Get Andalucia outline from rnaturalearth if available, otherwise use convex hull
station_sf <- st_as_sf(station_summary,
                        coords = c("longitud", "latitud"), crs = 4326)

p1 <- ggplot() +
  geom_sf(data = station_sf,
          aes(colour = n_years, size = median_max),
          alpha = 0.8) +
  scale_colour_viridis_c(name = "Record\nlength\n(years)", option = "C") +
  scale_size_continuous(name = "Median\nannual\nmax (mm)", range = c(1, 5)) +
  labs(title = "AEMET stations in Andalucia",
       subtitle = sprintf("%d stations, coloured by record length, sized by median annual maximum",
                          nrow(station_summary))) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")

ggsave("figures/station_map.pdf", p1, width = 10, height = 6)

# ---- Figure 2: Distribution of annual maxima --------------------------------
cat("Figure 2: Annual maxima histogram...\n")

p2 <- ggplot(annual_max, aes(x = max_prec)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, fill = "steelblue", alpha = 0.7) +
  geom_density(colour = "darkred", linewidth = 0.8) +
  labs(title = "Distribution of annual maximum daily precipitation",
       subtitle = sprintf("%d station-year observations across Andalucia",
                          nrow(annual_max)),
       x = "Annual maximum daily precipitation (mm)",
       y = "Density") +
  theme_minimal(base_size = 11)

ggsave("figures/annual_maxima_histogram.pdf", p2, width = 8, height = 5)

# ---- Figure 3: Record length distribution -----------------------------------
cat("Figure 3: Record length distribution...\n")

p3 <- ggplot(station_summary, aes(x = reorder(indicativo, -n_years),
                                    y = n_years, fill = provincia)) +
  geom_col() +
  labs(title = "Record length by station",
       subtitle = sprintf("%d stations (no minimum threshold)",
                          nrow(station_summary)),
       x = "Station (sorted by record length)",
       y = "Number of years with valid annual maxima",
       fill = "Province") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave("figures/record_length_dist.pdf", p3, width = 10, height = 5)

# ---- Figure 4: Boxplots by province ----------------------------------------
cat("Figure 4: Maxima by province...\n")

p4 <- ggplot(annual_max, aes(x = reorder(provincia, max_prec, FUN = median),
                               y = max_prec)) +
  geom_boxplot(aes(fill = provincia), alpha = 0.7, outlier.size = 0.8) +
  labs(title = "Annual maximum daily precipitation by province",
       x = "Province",
       y = "Annual max daily precipitation (mm)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  coord_flip()

ggsave("figures/maxima_by_province.pdf", p4, width = 8, height = 5)

# ---- Figure 5: Annual maxima vs elevation -----------------------------------
cat("Figure 5: Maxima vs elevation...\n")

p5 <- ggplot(station_summary, aes(x = altitud, y = median_max)) +
  geom_point(aes(colour = provincia, size = n_years), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
  labs(title = "Median annual maximum vs station elevation",
       x = "Elevation (m a.s.l.)",
       y = "Median annual max precip (mm)",
       colour = "Province",
       size = "Record\nlength") +
  theme_minimal(base_size = 11)

ggsave("figures/maxima_vs_elevation.pdf", p5, width = 8, height = 5)

# ---- Figure 6: Annual maxima vs lat/lon -------------------------------------
cat("Figure 6: Maxima vs coordinates...\n")

p6a <- ggplot(station_summary, aes(x = longitud, y = median_max)) +
  geom_point(aes(colour = provincia), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
  labs(x = "Longitude", y = "Median annual max (mm)",
       title = "East-west gradient") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

p6b <- ggplot(station_summary, aes(x = latitud, y = median_max)) +
  geom_point(aes(colour = provincia), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", linewidth = 0.5) +
  labs(x = "Latitude", y = "Median annual max (mm)",
       title = "North-south gradient") +
  theme_minimal(base_size = 10)

# Combine with patchwork if available, otherwise save separately
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p6 <- p6a + p6b + plot_layout(guides = "collect") +
    plot_annotation(title = "Spatial gradients in median annual maxima")
  ggsave("figures/maxima_vs_coords.pdf", p6, width = 10, height = 5)
} else {
  ggsave("figures/maxima_vs_lon.pdf", p6a, width = 5, height = 5)
  ggsave("figures/maxima_vs_lat.pdf", p6b, width = 5, height = 5)
}

# ---- Figure 7: Empirical variogram -----------------------------------------
cat("Figure 7: Empirical variogram...\n")

# Use median annual max per station as the summary for the variogram
vario_sf <- st_as_sf(station_summary,
                      coords = c("longitud", "latitud"), crs = 4326)
# Project to UTM zone 30N for distance in km
vario_utm <- st_transform(vario_sf, 25830)
vario_sp  <- as(vario_utm, "Spatial")

emp_vario <- variogram(median_max ~ 1, data = vario_sp)

p7 <- ggplot(emp_vario, aes(x = dist / 1000, y = gamma)) +
  geom_point(aes(size = np), alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, colour = "steelblue",
              linewidth = 0.8) +
  scale_size_continuous(name = "# pairs") +
  labs(title = "Empirical variogram of median annual maximum precipitation",
       x = "Distance (km)",
       y = "Semivariance") +
  theme_minimal(base_size = 11)

ggsave("figures/empirical_variogram.pdf", p7, width = 8, height = 5)

# ---- Figure 8: At-site GEV fits at selected stations -----------------------
cat("Figure 8: At-site GEV fits...\n")

# Select top 6 stations by record length
top_stations <- station_summary %>%
  arrange(desc(n_years)) %>%
  slice_head(n = 6)

gev_plots <- list()
gev_results <- list()

for (i in seq_len(nrow(top_stations))) {
  stn <- top_stations$indicativo[i]
  stn_name <- top_stations$nombre[i]
  stn_data <- annual_max %>% filter(indicativo == stn) %>% pull(max_prec)

  fit <- tryCatch(fgev(stn_data), error = function(e) NULL)

  if (!is.null(fit)) {
    gev_results[[stn]] <- list(
      name = stn_name,
      loc  = fit$estimate["loc"],
      scale = fit$estimate["scale"],
      shape = fit$estimate["shape"],
      n = length(stn_data)
    )

    # Return level plot data
    m <- c(2, 5, 10, 20, 50, 100)
    rl <- sapply(m, function(p) {
      mu <- fit$estimate["loc"]
      sig <- fit$estimate["scale"]
      xi <- fit$estimate["shape"]
      yp <- -log(1 - 1/p)
      if (abs(xi) < 1e-6) {
        mu - sig * log(yp)
      } else {
        mu + sig/xi * (yp^(-xi) - 1)
      }
    })

    # Empirical return levels (plotting positions)
    n_obs <- length(stn_data)
    sorted <- sort(stn_data)
    emp_prob <- seq_len(n_obs) / (n_obs + 1)
    emp_rp <- 1 / (1 - emp_prob)

    plot_df <- data.frame(
      return_period = emp_rp,
      observed = sorted
    )
    fitted_df <- data.frame(
      return_period = m,
      fitted = rl
    )

    gev_plots[[i]] <- ggplot() +
      geom_point(data = plot_df,
                 aes(x = return_period, y = observed),
                 alpha = 0.6) +
      geom_line(data = fitted_df,
                aes(x = return_period, y = fitted),
                colour = "red", linewidth = 0.8) +
      scale_x_log10() +
      labs(title = sprintf("%s (%d years)", stn_name, length(stn_data)),
           x = "Return period (years)", y = "Precip (mm)") +
      theme_minimal(base_size = 9)
  }
}

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p8 <- wrap_plots(gev_plots, ncol = 3) +
    plot_annotation(title = "At-site GEV return level plots (top 6 longest-record stations)")
  ggsave("figures/gev_fits_selected.pdf", p8, width = 12, height = 8)
} else {
  for (i in seq_along(gev_plots)) {
    ggsave(sprintf("figures/gev_fit_%d.pdf", i), gev_plots[[i]],
           width = 5, height = 4)
  }
}

# ---- Print GEV parameter summary -------------------------------------------
cat("\n=== At-site GEV parameter estimates (top stations) ===\n")
for (stn in names(gev_results)) {
  r <- gev_results[[stn]]
  cat(sprintf("  %s (%s): mu=%.1f, sigma=%.1f, xi=%.3f (n=%d)\n",
              stn, r$name, r$loc, r$scale, r$shape, r$n))
}

cat("\n=== Exploratory analysis complete ===\n")
