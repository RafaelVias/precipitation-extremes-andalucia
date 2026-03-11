# 11_grid_covariates.R — Precompute altitude and windward exposure on prediction grid
#
# Creates the same prediction grid used by scripts 04-06, then extracts DEM
# altitude and computes windward exposure at each grid cell.  The result is
# saved to data/grid_covariates.rds for use by prediction scripts.
#
# Wind directions and transect parameters match R/10_dem_exposure.R.
#
# Run from project root: Rscript R/11_grid_covariates.R

library(terra)
library(sf)
library(rnaturalearth)

cat("========================================\n")
cat("Grid-level covariates (altitude + exposure)\n")
cat("========================================\n")

# ---- 1. Create prediction grid (same as scripts 04-06) ----
cat("Step 1: Creating prediction grid...\n")

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

pred_res <- 0.005
bbox <- st_bbox(andalucia)
pred_pts <- expand.grid(
  lon = seq(bbox["xmin"] + pred_res / 2, bbox["xmax"], by = pred_res),
  lat = seq(bbox["ymin"] + pred_res / 2, bbox["ymax"], by = pred_res)
)
pred_sf <- st_as_sf(pred_pts, coords = c("lon", "lat"), crs = 4326)
inside <- st_intersects(pred_sf, andalucia, sparse = FALSE)[, 1]
pred_pts <- pred_pts[inside, ]
n_pred <- nrow(pred_pts)

cat(sprintf("  Grid: %d points (res = %.3f°)\n", n_pred, pred_res))

# ---- 2. Load DEM ----
cat("Step 2: Loading DEM...\n")

dem_file <- "data/dem/andalucia_srtm.tif"
if (!file.exists(dem_file)) {
  stop("DEM not found at ", dem_file, ". Run R/10_dem_exposure.R first.")
}
dem <- rast(dem_file)
cat(sprintf("  DEM: %d x %d cells\n", ncol(dem), nrow(dem)))

# ---- 3. Extract altitude at grid points ----
cat("Step 3: Extracting altitude at grid points...\n")

grid_vect <- vect(as.matrix(pred_pts), crs = "EPSG:4326")
alt_grid <- terra::extract(dem, grid_vect)[[2]]

cat(sprintf("  Altitude range: %.0f - %.0f m\n",
            min(alt_grid, na.rm = TRUE), max(alt_grid, na.rm = TRUE)))
cat(sprintf("  NAs: %d\n", sum(is.na(alt_grid))))

# Set coastal/sea-level NAs to 0
alt_grid[is.na(alt_grid)] <- 0

# ---- 4. Compute exposure at grid points ----
cat("Step 4: Computing windward exposure at grid points...\n")

# Same directions and transect parameters as R/10_dem_exposure.R
wind_dirs <- c(mediterranean_se = 135, atlantic_wsw = 255)
D_km <- 20
n_transect <- 40
D_deg <- D_km / 111

compute_exposure <- function(lon_s, lat_s, alt_s, azimuth_deg, dem_rast,
                             d_deg, n_pts) {
  az_rad <- azimuth_deg * pi / 180
  dists <- seq(d_deg / n_pts, d_deg, length.out = n_pts)
  t_lon <- lon_s + dists * sin(az_rad)
  t_lat <- lat_s + dists * cos(az_rad)
  t_pts <- vect(cbind(t_lon, t_lat), crs = "EPSG:4326")
  t_alt <- terra::extract(dem_rast, t_pts)[[2]]
  mean_upwind <- mean(t_alt, na.rm = TRUE)
  if (is.na(mean_upwind)) return(NA_real_)
  alt_s - mean_upwind
}

exposure_grid <- list()

for (dir_name in names(wind_dirs)) {
  az <- wind_dirs[dir_name]
  cat(sprintf("  Direction: %s (azimuth %d°)...\n", dir_name, az))
  t0 <- proc.time()

  exp_vals <- sapply(seq_len(n_pred), function(i) {
    if (i %% 5000 == 0) cat(sprintf("    %d / %d\n", i, n_pred))
    compute_exposure(pred_pts$lon[i], pred_pts$lat[i], alt_grid[i],
                     az, dem, D_deg, n_transect)
  })

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("    Done in %.0f seconds\n", elapsed))
  exposure_grid[[dir_name]] <- exp_vals
}

exp_mean_grid <- (exposure_grid$mediterranean_se + exposure_grid$atlantic_wsw) / 2

cat(sprintf("  Mean exposure range: %.0f to %.0f m\n",
            min(exp_mean_grid, na.rm = TRUE), max(exp_mean_grid, na.rm = TRUE)))

# Replace NAs with 0 (coastal cells where transect goes over sea)
exposure_grid$mediterranean_se[is.na(exposure_grid$mediterranean_se)] <- 0
exposure_grid$atlantic_wsw[is.na(exposure_grid$atlantic_wsw)] <- 0
exp_mean_grid[is.na(exp_mean_grid)] <- 0

# ---- 5. Save ----
cat("Step 5: Saving grid covariates...\n")

grid_cov <- data.frame(
  lon          = pred_pts$lon,
  lat          = pred_pts$lat,
  alt_dem      = alt_grid,
  exposure_mediterranean_se = exposure_grid$mediterranean_se,
  exposure_atlantic_wsw     = exposure_grid$atlantic_wsw,
  exposure_mean             = exp_mean_grid
)

saveRDS(grid_cov, "data/grid_covariates.rds")
cat(sprintf("  Saved data/grid_covariates.rds (%d grid points)\n", n_pred))

cat("\n========================================\n")
cat("Done. Grid covariates ready.\n")
cat("========================================\n")
