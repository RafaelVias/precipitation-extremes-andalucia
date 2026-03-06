# =============================================================================
# 01_acquire_data.R
# Download daily precipitation data from AEMET for Andalucia
# =============================================================================
#
# Prerequisites:
#   - climaemet package installed: install.packages("climaemet")
#   - AEMET API key registered:
#       climaemet::aemet_api_key("YOUR_KEY", install = TRUE)
#     Get a free key at: https://opendata.aemet.es/centrodedescargas/obtencionAPIKey
#
# Outputs:
#   - data/stations_andalucia.rds      Station inventory for Andalucia
#   - data/stations_andalucia.csv      Same, in CSV
#   - data/daily_precip_andalucia_raw.rds  Raw daily climate data
#   - data/annual_maxima_andalucia.rds  Annual block maxima with QC
#
# NOTE: The AEMET API is rate-limited (50 requests/min). This script uses
#       aemet_daily_period_all() which fetches all stations in one call per
#       year, then filters locally. Full download takes ~5-10 minutes.
# =============================================================================

library(climaemet)
library(dplyr)
library(lubridate)

# ---- Configuration ----------------------------------------------------------
ANDALUCIA_PROVINCES <- c("MALAGA", "SEVILLA", "CORDOBA", "GRANADA",
                          "JAEN", "ALMERIA", "CADIZ", "HUELVA")
START_YEAR <- 1960
END_YEAR   <- as.integer(format(Sys.Date(), "%Y"))
MIN_DAYS   <- 330   # minimum days per year (90% completeness)
API_SLEEP  <- 3     # seconds between API calls

# ---- Step 1: Station inventory ----------------------------------------------
cat("Step 1: Downloading station inventory...\n")

stations_all <- aemet_stations()
stations_and <- stations_all %>%
  filter(provincia %in% ANDALUCIA_PROVINCES)

cat("  Total AEMET stations:", nrow(stations_all), "\n")
cat("  Andalucia stations:", nrow(stations_and), "\n")
cat("  By province:\n")
print(stations_and %>% count(provincia, sort = TRUE), n = 8)

saveRDS(stations_and, "data/stations_andalucia.rds")
write.csv(stations_and, "data/stations_andalucia.csv", row.names = FALSE)

# ---- Step 2: Download daily data (year by year) -----------------------------
cat("\nStep 2: Downloading daily data ", START_YEAR, "-", END_YEAR, "...\n")

and_ids <- stations_and$indicativo
all_data <- list()

for (yr in START_YEAR:END_YEAR) {
  cat(sprintf("  %d ... ", yr))

  chunk <- tryCatch(
    aemet_daily_period_all(start = yr, end = yr),
    error = function(e) {
      cat("ERROR:", conditionMessage(e), " ")
      NULL
    }
  )

  if (!is.null(chunk) && nrow(chunk) > 0) {
    chunk_and <- chunk %>% filter(indicativo %in% and_ids)
    if (nrow(chunk_and) > 0) {
      all_data[[as.character(yr)]] <- chunk_and
      cat(nrow(chunk_and), "rows,",
          n_distinct(chunk_and$indicativo), "stations\n")
    } else {
      cat("no Andalucia data\n")
    }
  } else {
    cat("no data\n")
  }

  Sys.sleep(API_SLEEP)
}

daily_all <- bind_rows(all_data)
cat("\n  Total rows:", nrow(daily_all), "\n")
cat("  Date range:", as.character(min(daily_all$fecha, na.rm = TRUE)), "to",
    as.character(max(daily_all$fecha, na.rm = TRUE)), "\n")
cat("  Unique stations:", n_distinct(daily_all$indicativo), "\n")

saveRDS(daily_all, "data/daily_precip_andalucia_raw.rds")
cat("  Saved to data/daily_precip_andalucia_raw.rds\n")

# ---- Step 3: Compute annual block maxima with QC ----------------------------
cat("\nStep 3: Computing annual block maxima...\n")

annual_max <- daily_all %>%
  mutate(year = year(fecha)) %>%
  group_by(indicativo, year) %>%
  summarise(
    max_prec = max(prec, na.rm = TRUE),
    n_obs    = sum(!is.na(prec)),
    .groups  = "drop"
  ) %>%
  # QC: require 90% completeness and finite max
  filter(n_obs >= MIN_DAYS, is.finite(max_prec), max_prec >= 0)

# No minimum-year threshold: the Bayesian spatial model borrows strength
# across stations, so even short-record stations contribute useful information.

# Join station metadata
annual_max <- annual_max %>%
  left_join(
    stations_and %>% select(indicativo, nombre, provincia, altitud,
                             latitud, longitud),
    by = "indicativo"
  )

cat("  Stations:", n_distinct(annual_max$indicativo), "\n")
cat("  Total station-year observations:", nrow(annual_max), "\n")
cat("  By province:\n")
print(annual_max %>%
        distinct(indicativo, provincia) %>%
        count(provincia, sort = TRUE), n = 8)

saveRDS(annual_max, "data/annual_maxima_andalucia.rds")
cat("  Saved to data/annual_maxima_andalucia.rds\n")

cat("\n=== Acquisition complete ===\n")
