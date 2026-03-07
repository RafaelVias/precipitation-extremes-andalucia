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

# Handle AEMET special codes before computing maxima:
#   "Ip"   = inapreciable (trace, < 0.01mm) — sample from U(0, 0.01)
#   "Acum" = accumulated over multiple days — distribute via Dirichlet(α)
# For Acum sequences: the reporting day after "Acum" days carries the multi-day
# total. We distribute it across all days in the window (Acum days + reporting
# day) using a Dirichlet(0.73, ..., 0.73) draw. α = 0.73 was estimated from
# 41,435 wet spells across 126 stations (bootstrap 95% CI: [0.727, 0.737]).
# With α < 1, the Dirichlet concentrates mass on fewer days — realistic for
# Mediterranean precipitation where rain events are typically concentrated.
daily_all <- daily_all %>%
  arrange(indicativo, fecha)

set.seed(1)

# --- Acum: Dirichlet disaggregation ---
is_acum <- daily_all$prec == "Acum" & !is.na(daily_all$prec)
n_acum <- sum(is_acum)
n_acum_events <- 0L
if (n_acum > 0) {
  ALPHA_DIRICHLET <- 0.73
  acum_idx <- which(is_acum)

  # Group consecutive Acum rows within the same station into events
  # Each event = sequence of Acum days + reporting day (first non-Acum after)
  processed <- rep(FALSE, length(acum_idx))
  for (k in seq_along(acum_idx)) {
    if (processed[k]) next
    # Start of a new Acum sequence
    seq_start <- acum_idx[k]
    stn <- daily_all$indicativo[seq_start]
    seq_end <- seq_start
    # Extend to consecutive Acum rows for the same station
    while (seq_end + 1 <= nrow(daily_all) &&
           daily_all$indicativo[seq_end + 1] == stn &&
           !is.na(daily_all$prec[seq_end + 1]) &&
           daily_all$prec[seq_end + 1] == "Acum") {
      seq_end <- seq_end + 1
    }
    # Mark all Acum rows in this sequence as processed
    processed[acum_idx >= seq_start & acum_idx <= seq_end] <- TRUE

    # Reporting day = next row after the Acum sequence
    report_i <- seq_end + 1
    if (report_i > nrow(daily_all) ||
        daily_all$indicativo[report_i] != stn ||
        is.na(daily_all$prec[report_i]) ||
        daily_all$prec[report_i] == "Acum") next

    total_mm <- as.numeric(gsub(",", ".", daily_all$prec[report_i]))
    if (is.na(total_mm) || total_mm < 0) next

    window_idx <- seq_start:report_i  # Acum days + reporting day
    n_days <- length(window_idx)

    # Dirichlet draw: sample k gamma(α, 1) values and normalise
    gdraws <- rgamma(n_days, shape = ALPHA_DIRICHLET, rate = 1)
    props <- gdraws / sum(gdraws)
    daily_vals <- total_mm * props

    daily_all$prec[window_idx] <- as.character(round(daily_vals, 2))
    n_acum_events <- n_acum_events + 1L

    cat(sprintf("    Acum: %s %s–%s: %.1fmm / %d days -> [%s]\n",
                stn, daily_all$fecha[seq_start], daily_all$fecha[report_i],
                total_mm, n_days,
                paste(round(daily_vals, 1), collapse = ", ")))
  }
}

# --- Ip: sample from U(0, 0.01) ---
ip_mask <- daily_all$prec == "Ip" & !is.na(daily_all$prec)
n_ip <- sum(ip_mask)
daily_all$prec[ip_mask] <- as.character(runif(n_ip, 0, 0.01))

daily_all <- daily_all %>%
  mutate(prec = as.numeric(gsub(",", ".", prec)))

cat(sprintf("  AEMET codes: %d Ip -> U(0,0.01) | %d Acum days (%d events) -> Dirichlet(0.73)\n",
            n_ip, n_acum, n_acum_events))

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
