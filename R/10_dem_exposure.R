# 10_dem_exposure.R — Acquire DEM and compute windward exposure index
#
# Downloads SRTM elevation data for Andalucía, computes a directional
# windward exposure index for two dominant moisture sources:
#   - Atlantic south-westerlies (azimuth ~225°)
#   - Mediterranean easterlies  (azimuth ~90°)
#
# Exposure index at station s for wind direction θ:
#   E(s, θ) = alt(s) − mean(alt along upwind transect of length D)
#   Positive = station is higher than upwind terrain (windward slope)
#   Negative = station is lower than upwind terrain (leeward/sheltered)
#
# Run from project root: Rscript R/10_dem_exposure.R

library(terra)
library(elevatr)
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)

cat("========================================\n")
cat("DEM acquisition & windward exposure\n")
cat("========================================\n")

dir.create("data/dem", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/exploration", showWarnings = FALSE, recursive = TRUE)

# ---- 1. Define Andalucía bounding box with buffer ----
# Need buffer for upwind transects that extend beyond the region
cat("Step 1: Setting up region...\n")

library(rnaturalearth)
states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

# Bounding box with 0.5° buffer for transect computation
bbox <- st_bbox(andalucia)
bbox_buffered <- st_bbox(c(
  xmin = unname(bbox["xmin"]) - 0.5,
  xmax = unname(bbox["xmax"]) + 0.5,
  ymin = unname(bbox["ymin"]) - 0.5,
  ymax = unname(bbox["ymax"]) + 0.5
), crs = 4326)

bbox_sf <- st_as_sfc(bbox_buffered) |> st_sf()

cat(sprintf("  Bounding box: [%.2f, %.2f] x [%.2f, %.2f]\n",
            bbox_buffered["xmin"], bbox_buffered["xmax"],
            bbox_buffered["ymin"], bbox_buffered["ymax"]))

# ---- 2. Download DEM ----
dem_file <- "data/dem/andalucia_srtm.tif"

if (file.exists(dem_file)) {
  cat("Step 2: Loading cached DEM...\n")
  dem <- rast(dem_file)
} else {
  cat("Step 2: Downloading SRTM DEM (this may take a moment)...\n")
  # z = 9 gives ~250m resolution, good balance of detail and file size
  dem_elev <- get_elev_raster(bbox_sf, z = 9, src = "aws")
  dem <- rast(dem_elev)
  writeRaster(dem, dem_file, overwrite = TRUE)
  cat("  Saved to", dem_file, "\n")
}

cat(sprintf("  DEM: %d x %d cells, resolution ~%.0fm\n",
            ncol(dem), nrow(dem), res(dem)[1] * 111000))

# ---- 3. Load station data ----
cat("Step 3: Loading stations...\n")

am <- readRDS("data/annual_maxima_andalucia.rds")
am <- am %>% filter(indicativo != "6381")

stn <- am %>%
  group_by(indicativo) %>%
  summarise(
    n_years    = n(),
    median_max = median(max_prec),
    mean_max   = mean(max_prec),
    alt_aemet  = first(as.numeric(altitud)),
    lon        = first(longitud),
    lat        = first(latitud),
    prov       = first(provincia),
    nombre     = first(nombre),
    .groups    = "drop"
  )

# Extract DEM altitude at station locations
stn_pts <- vect(stn, geom = c("lon", "lat"), crs = "EPSG:4326")
stn$alt_dem <- terra::extract(dem, stn_pts)[, 2]

cat(sprintf("  Stations: %d\n", nrow(stn)))
cat(sprintf("  AEMET vs DEM altitude correlation: r = %.3f\n",
            cor(stn$alt_aemet, stn$alt_dem, use = "complete")))
cat(sprintf("  AEMET vs DEM altitude MAE: %.0f m\n",
            mean(abs(stn$alt_aemet - stn$alt_dem), na.rm = TRUE)))

# ---- 4. Compute windward exposure index ----
cat("Step 4: Computing windward exposure indices...\n")

# Wind directions (azimuth in degrees, clockwise from north):
#   Mediterranean SE:  135° (levante moisture from the SE)
#   Atlantic WSW:      255° (Atlantic fronts from the WSW)
# Directions chosen by sweeping all azimuths in 15° steps and selecting
# the pair whose mean exposure maximises correlation with median annual
# maxima (r = +0.335 for 135° + 255°, vs +0.279 for the initial 90° + 225°).
wind_dirs <- c(mediterranean_se = 135, atlantic_wsw = 255)

# Transect parameters
D_km <- 20            # transect length in km
n_transect <- 40      # number of points along transect
D_deg <- D_km / 111   # approximate conversion to degrees

compute_exposure <- function(lon_s, lat_s, alt_s, azimuth_deg, dem_rast,
                             d_deg, n_pts) {
  # Azimuth = direction FROM which wind blows (clockwise from north)
  # We walk FROM the station TOWARD the wind source along the transect
  # Standard convention: east = sin(az), north = cos(az)
  az_rad <- azimuth_deg * pi / 180

  dists <- seq(d_deg / n_pts, d_deg, length.out = n_pts)
  t_lon <- lon_s + dists * sin(az_rad)
  t_lat <- lat_s + dists * cos(az_rad)

  # Extract elevations along transect
  t_pts <- vect(cbind(t_lon, t_lat), crs = "EPSG:4326")
  t_alt <- terra::extract(dem_rast, t_pts)[[2]]

  # Exposure = station altitude - mean upwind altitude
  mean_upwind <- mean(t_alt, na.rm = TRUE)
  if (is.na(mean_upwind)) return(NA_real_)
  alt_s - mean_upwind
}

# Compute for each station and each wind direction
for (dir_name in names(wind_dirs)) {
  az <- wind_dirs[dir_name]
  cat(sprintf("  Direction: %s (azimuth %d°)...\n", dir_name, az))

  exp_vals <- sapply(seq_len(nrow(stn)), function(i) {
    compute_exposure(stn$lon[i], stn$lat[i], stn$alt_dem[i],
                     az, dem, D_deg, n_transect)
  })
  stn[[paste0("exposure_", dir_name)]] <- exp_vals
}

# Combined exposure
stn$exposure_max <- pmax(stn$exposure_mediterranean_se, stn$exposure_atlantic_wsw)
stn$exposure_mean <- (stn$exposure_mediterranean_se + stn$exposure_atlantic_wsw) / 2

cat(sprintf("  Mediterranean SE exposure range: %.0f to %.0f m\n",
            min(stn$exposure_mediterranean_se, na.rm = TRUE),
            max(stn$exposure_mediterranean_se, na.rm = TRUE)))
cat(sprintf("  Atlantic WSW exposure range: %.0f to %.0f m\n",
            min(stn$exposure_atlantic_wsw, na.rm = TRUE),
            max(stn$exposure_atlantic_wsw, na.rm = TRUE)))

# ---- 5. Correlations ----
cat("\nCorrelations with median annual maximum:\n")
vars <- c("alt_dem", "exposure_mediterranean_se", "exposure_atlantic_wsw",
           "exposure_max", "exposure_mean")
for (v in vars) {
  r <- cor(stn$median_max, stn[[v]], use = "complete")
  cat(sprintf("  %-30s: r = %+.3f\n", v, r))
}

# ---- 6. Save station data with exposure ----
saveRDS(stn, "data/station_exposure.rds")
cat("\n  Saved data/station_exposure.rds\n")

# ---- 7. Exploratory plots ----
cat("\nGenerating exploratory plots...\n")

prov_colours <- c(
  "ALMERIA"  = "#E41A1C", "CADIZ"    = "#377EB8", "CORDOBA"  = "#4DAF4A",
  "GRANADA"  = "#984EA3", "HUELVA"   = "#FF7F00", "JAEN"     = "#A65628",
  "MALAGA"   = "#F781BF", "SEVILLA"  = "#999999"
)

# Plot 1: Median annual max vs Atlantic WSW exposure
p_atl <- ggplot(stn, aes(x = exposure_atlantic_wsw, y = median_max, colour = prov)) +
  geom_point(aes(size = n_years), alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              colour = "grey30", fill = "grey80", linewidth = 0.8, alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = prov_colours, name = "Province") +
  scale_size_continuous(name = "Years", range = c(1.5, 5), breaks = c(10, 30, 60)) +
  labs(
    title = "Atlantic WSW exposure (azimuth 255\u00b0)",
    subtitle = sprintf("r = %+.3f | transect = %d km",
                       cor(stn$median_max, stn$exposure_atlantic_wsw, use = "complete"), D_km),
    x = "Exposure index (m)", y = "Median annual maximum (mm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank())

# Plot 2: Median annual max vs Mediterranean SE exposure
p_med <- ggplot(stn, aes(x = exposure_mediterranean_se, y = median_max, colour = prov)) +
  geom_point(aes(size = n_years), alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              colour = "grey30", fill = "grey80", linewidth = 0.8, alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = prov_colours, name = "Province") +
  scale_size_continuous(name = "Years", range = c(1.5, 5), breaks = c(10, 30, 60)) +
  labs(
    title = "Mediterranean SE exposure (azimuth 135\u00b0)",
    subtitle = sprintf("r = %+.3f | transect = %d km",
                       cor(stn$median_max, stn$exposure_mediterranean_se, use = "complete"), D_km),
    x = "Exposure index (m)", y = "Median annual maximum (mm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank())

# Plot 3: Mean exposure vs median max
p_mean_exp <- ggplot(stn, aes(x = exposure_mean, y = median_max, colour = prov)) +
  geom_point(aes(size = n_years), alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "loess", se = TRUE,
              colour = "grey30", fill = "grey80", linewidth = 0.8, alpha = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = prov_colours, name = "Province") +
  scale_size_continuous(name = "Years", range = c(1.5, 5), breaks = c(10, 30, 60)) +
  labs(
    title = "Mean exposure (average of both directions)",
    subtitle = sprintf("r = %+.3f",
                       cor(stn$median_max, stn$exposure_mean, use = "complete")),
    x = "Exposure index (m)", y = "Median annual maximum (mm)"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank())

# Plot 4: DEM map with station exposure
# Create a simplified DEM for plotting
dem_df <- as.data.frame(dem, xy = TRUE)
names(dem_df)[3] <- "elevation"
dem_df <- dem_df[!is.na(dem_df$elevation), ]

# Crop to Andalucía bbox
dem_df <- dem_df[dem_df$x >= bbox["xmin"] - 0.1 & dem_df$x <= bbox["xmax"] + 0.1 &
                 dem_df$y >= bbox["ymin"] - 0.1 & dem_df$y <= bbox["ymax"] + 0.1, ]

p_map <- ggplot() +
  geom_raster(data = dem_df, aes(x = x, y = y, fill = elevation)) +
  scale_fill_gradientn(
    colours = c("#1a5276", "#2e86c1", "#85c1e9", "#abebc6",
                "#f9e79f", "#e67e22", "#922b21", "#7b241c"),
    name = "Elevation (m)",
    limits = c(0, max(dem_df$elevation)),
    na.value = "grey90"
  ) +
  geom_sf(data = andalucia, fill = NA, colour = "grey20", linewidth = 0.5) +
  geom_point(data = stn, aes(x = lon, y = lat, size = exposure_mean),
             colour = "black", fill = "white", shape = 21, stroke = 0.5) +
  scale_size_continuous(name = "Mean exposure\n(m)", range = c(1, 6)) +
  coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
  labs(
    title = "DEM and station windward exposure",
    subtitle = sprintf("SRTM ~%.0fm | transect length = %d km", res(dem)[1] * 111000, D_km)
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        panel.grid.minor = element_blank())

# Combine
p_combined <- (p_atl + p_med) / (p_mean_exp + p_map) +
  plot_annotation(
    title = "Windward exposure and extreme daily rainfall in Andaluc\u00eda",
    subtitle = sprintf("Exposure = station altitude \u2212 mean upwind altitude over %d km transect | positive = windward slope", D_km),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  )

ggsave("figures/exploration/exposure_vs_maxima.png", p_combined,
       width = 18, height = 14, dpi = 200, bg = "white")
cat("  Saved figures/exploration/exposure_vs_maxima.png\n")

# Individual panels
ggsave("figures/exploration/exposure_atlantic.png", p_atl,
       width = 10, height = 7, dpi = 200, bg = "white")
ggsave("figures/exploration/exposure_mediterranean.png", p_med,
       width = 10, height = 7, dpi = 200, bg = "white")
ggsave("figures/exploration/dem_map.png", p_map,
       width = 12, height = 7, dpi = 200, bg = "white")

cat("\n========================================\n")
cat("Done.\n")
cat("========================================\n")
