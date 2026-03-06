# =============================================================================
# 02_stage1_mle.R
# Stage 1: Poisson point process GEV MLEs at each station
# =============================================================================
#
# Adapts the code from:
#   https://github.com/arnabstatswithR/max_and_smooth
#   (Hazra, Huser, Jóhannesson - BLGM Chapter 7)
#
# Inputs:
#   - data/daily_precip_andalucia_raw.rds
#   - data/stations_andalucia.rds
#   - vendor/max_and_smooth/stage1_functions.R
#
# Outputs:
#   - data/stage1_results.rds         Site-wise MLEs and covariance matrices
# =============================================================================

library(dplyr)
library(lubridate)

# ---- Source the M&S functions ------------------------------------------------
source("vendor/max_and_smooth/stage1_functions.R")

# =============================================================================
# STEP 0: Prepare daily data as station x day matrix
# =============================================================================
cat("Step 0: Preparing daily data matrix...\n")

daily_all <- readRDS("data/daily_precip_andalucia_raw.rds")
stations  <- readRDS("data/stations_andalucia.rds")

# Fix European decimal notation
daily_all <- daily_all %>%
  mutate(prec = as.numeric(gsub(",", ".", prec)),
         fecha = as.Date(fecha))

# Keep only stations that appear in our annual maxima dataset (passed QC)
annual_max <- readRDS("data/annual_maxima_andalucia.rds")
valid_stations <- unique(annual_max$indicativo)

daily_clean <- daily_all %>%
  filter(indicativo %in% valid_stations,
         !is.na(prec),
         prec >= 0)

# Get station metadata (coordinates) in order
station_meta <- stations %>%
  filter(indicativo %in% valid_stations) %>%
  arrange(indicativo) %>%
  select(indicativo, nombre, provincia, altitud, latitud, longitud)

station_ids <- station_meta$indicativo
ns <- length(station_ids)

cat("  Stations:", ns, "\n")

# Create the Y matrix: rows = stations, columns = days
# Each station has different date ranges, so we use the union of all dates
all_dates <- sort(unique(daily_clean$fecha))
nt <- length(all_dates)
cat("  Date range:", as.character(min(all_dates)), "to",
    as.character(max(all_dates)), "(", nt, "days)\n")

# Build matrix efficiently
cat("  Building Y matrix (", ns, "x", nt, ")...\n")
Y <- matrix(NA, nrow = ns, ncol = nt)
date_idx <- match(daily_clean$fecha, all_dates)
stn_idx  <- match(daily_clean$indicativo, station_ids)

for (i in seq_len(nrow(daily_clean))) {
  si <- stn_idx[i]
  di <- date_idx[i]
  if (!is.na(si) && !is.na(di)) {
    Y[si, di] <- daily_clean$prec[i]
  }
}

# Keep Y with NAs: actual obs are numeric, missing days are NA
# For each station, the actual observation vector is Y[i, !is.na(Y[i,])]
# This is critical for correct nt in the PPP likelihood

cat("  Y matrix dimensions:", dim(Y), "\n")
cat("  Non-NA entries:", sum(!is.na(Y)), "\n")
cat("  Positive entries:", sum(Y > 0, na.rm = TRUE), "\n")

# Compute station-specific observation counts
nt_per_station <- apply(Y, 1, function(x) sum(!is.na(x)))
cat("  Observation days per station: min=", min(nt_per_station),
    " median=", median(nt_per_station), " max=", max(nt_per_station), "\n")

# Station coordinates matrix
loc <- as.matrix(station_meta[, c("longitud", "latitud")])
colnames(loc) <- c("lon", "lat")

# ---- Merge co-located station pairs ---
# Pairs at < 0.01 deg apart are the same physical station with different IDs
# (typically renumbered across time periods). Merge their daily records to get
# longer time series. On overlapping days, prefer the station with more data.
cat("\n  Merging co-located station pairs (< 0.01 deg)...\n")
dists <- as.matrix(dist(loc))
diag(dists) <- Inf
merged_away <- c()  # indices of secondary stations (to remove after merging)
for (i in 1:(ns - 1)) {
  if (i %in% merged_away) next
  for (j in (i + 1):ns) {
    if (j %in% merged_away) next
    if (dists[i, j] < 0.01) {
      # Determine primary (more data) and secondary
      if (nt_per_station[i] >= nt_per_station[j]) {
        primary <- i; secondary <- j
      } else {
        primary <- j; secondary <- i
      }

      # Merge: fill NAs in primary with secondary's values
      na_in_primary <- is.na(Y[primary, ]) & !is.na(Y[secondary, ])
      n_filled <- sum(na_in_primary)
      Y[primary, na_in_primary] <- Y[secondary, na_in_primary]
      nt_new <- sum(!is.na(Y[primary, ]))

      cat(sprintf("    %s (%s) + %s (%s) -> %s (%d + %d new = %d days)\n",
          station_meta$nombre[primary], station_ids[primary],
          station_meta$nombre[secondary], station_ids[secondary],
          station_ids[primary], nt_per_station[primary], n_filled, nt_new))

      nt_per_station[primary] <- nt_new
      merged_away <- c(merged_away, secondary)
    }
  }
}

if (length(merged_away) > 0) {
  keep_idx <- setdiff(1:ns, merged_away)
  Y <- Y[keep_idx, ]
  loc <- loc[keep_idx, ]
  station_meta <- station_meta[keep_idx, ]
  station_ids <- station_ids[keep_idx]
  nt_per_station <- nt_per_station[keep_idx]
  ns <- length(keep_idx)
  cat("  After merging:", ns, "stations\n")
}

# ---- Remove non-mainland stations ---
# Alborán (6381): island 50km off mainland, between Spain and Morocco
problem_ids <- c("6381")
problem_idx <- which(station_ids %in% problem_ids)
if (length(problem_idx) > 0) {
  cat("\n  Removing non-mainland stations:\n")
  for (pi in problem_idx) {
    cat(sprintf("    %s (%s)\n", station_meta$nombre[pi], station_ids[pi]))
  }
  keep_idx <- setdiff(1:ns, problem_idx)
  Y <- Y[keep_idx, ]
  loc <- loc[keep_idx, ]
  station_meta <- station_meta[keep_idx, ]
  station_ids <- station_ids[keep_idx]
  nt_per_station <- nt_per_station[keep_idx]
  ns <- length(keep_idx)
  cat("  After removals:", ns, "stations\n")
}

# =============================================================================
# STAGE 1: Poisson point process MLEs at each station
# =============================================================================
cat("\n========================================\n")
cat("STAGE 1: Computing site-wise MLEs...\n")
cat("========================================\n")

# Site-specific thresholds: 75th percentile of positive precipitation
# (following Hazra et al. 2023)
# Use only actual observations (not NA-padded zeros!)
thresholds <- apply(Y, 1, function(x) {
  obs <- x[!is.na(x)]
  pos <- obs[obs > 0]
  if (length(pos) > 10) quantile(pos, probs = 0.75) else NA
})

cat("  Threshold range:", round(min(thresholds, na.rm = TRUE), 1), "to",
    round(max(thresholds, na.rm = TRUE), 1), "mm\n")

# Number of observations per year for npy parameter
npy <- 365.25

# Fit point process MLEs at each station
# CRITICAL: pass only actual observation days per station (not zero-padded row)
cat("  Fitting MLEs at", ns, "stations...\n")

ppfit.mles <- list()
failed <- c()

for (i in 1:ns) {
  if (i %% 20 == 0) cat("    Station", i, "/", ns, "\n")

  # Extract only actual observations for this station
  obs_i <- Y[i, !is.na(Y[i, ])]

  result <- tryCatch({
    pp.fit.fast(obs_i, threshold = thresholds[i], npy = npy)
  }, error = function(e) {
    NULL
  })

  if (!is.null(result) && result$conv == 0) {
    ppfit.mles[[i]] <- result
  } else {
    failed <- c(failed, i)
    ppfit.mles[[i]] <- NULL
  }
}

n_failed <- length(failed)
cat("  Successful fits:", ns - n_failed, "/", ns, "\n")
if (n_failed > 0) {
  cat("  Failed stations:", paste(station_ids[failed], collapse = ", "), "\n")
}

# Extract MLEs
mles <- t(sapply(seq_along(ppfit.mles), function(i) {
  if (!is.null(ppfit.mles[[i]])) ppfit.mles[[i]]$mle else rep(NA, 3)
}))
colnames(mles) <- c("psi", "tau", "phi")

# Transform to original scale for display
mles.original <- mles
mles.original[, 1] <- exp(mles[, 1])        # mu
mles.original[, 2] <- exp(mles[, 1] + mles[, 2])  # sigma
mles.original[, 3] <- g(mles[, 3])          # xi
colnames(mles.original) <- c("mu", "sigma", "xi")

cat("\n  MLE summary (original scale):\n")
print(summary(mles.original))

# Remove failed stations from all objects
if (n_failed > 0) {
  keep <- setdiff(1:ns, failed)
  Y <- Y[keep, ]
  loc <- loc[keep, ]
  station_meta <- station_meta[keep, ]
  thresholds <- thresholds[keep]
  nt_per_station <- nt_per_station[keep]
  ppfit.mles <- ppfit.mles[keep]
  mles <- mles[keep, ]
  mles.original <- mles.original[keep, ]
  station_ids <- station_ids[keep]
  ns <- length(keep)
  cat("  After removing failed stations:", ns, "stations\n")
}

# Parametric bootstrap for covariance matrices
cat("\n  Computing bootstrap covariance matrices (1000 samples per station)...\n")
cat("  This will take several minutes...\n")

covmats <- list()
for (i in 1:ns) {
  if (i %% 20 == 0) cat("    Station", i, "/", ns, "\n")

  mu.station <- mles.original[i, 1]
  sigma.station <- mles.original[i, 2]
  xi.station <- mles.original[i, 3]
  threshold.station <- thresholds[i]
  tau.station <- sigma.station + xi.station * (threshold.station - mu.station)
  # Use station-specific nt (actual observation days)
  nt.station <- nt_per_station[i]
  lambda.station <- nt.station / npy *
    (1 + xi.station * (threshold.station - mu.station) / sigma.station)^(-1/xi.station)

  log.mu.init <- mles[i, 1]
  log.sigma.by.mu.init <- mles[i, 2]
  phi.init <- mles[i, 3]
  init <- c(log.mu.init, log.sigma.by.mu.init, phi.init)

  est.PPP <- function(rep.no) {
    set.seed(rep.no)
    n.exceed <- rpois(1, lambda.station)
    if (n.exceed < 3) return(rep(NA, 3))

    sample.PPP <- threshold.station + tau.station / xi.station *
      (runif(n.exceed)^{-xi.station} - 1)

    pp.lik <- function(a) {
      mu <- exp(a[1])
      sc <- exp(a[1] + a[2])
      xi <- g(a[3])

      if ((1 + xi * (threshold.station - mu)/sc) < 0) { l <- 10^6 } else {
        y <- (sample.PPP - mu)/sc
        y <- 1 + xi * y
        if (min(y) <= 0) { l <- 10^6 } else {
          ll <- log(sc) * n.exceed + sum(log(y)) * (1/xi + 1) +
            nt.station/npy * (1 + xi * (threshold.station - mu)/sc)^(-1/xi)
          pp <- (4 - c.phi) * log(xi - xi_lo) + (4 - 1) * log(xi_hi - xi) +
            (a[3] - a.phi) / b.phi - exp((a[3] - a.phi) / b.phi)
          l <- ll - pp
        }
      }
      l
    }

    x <- tryCatch(
      optim(init, pp.lik, control = list(maxit = 1000)),
      error = function(e) list(convergence = 1)
    )
    if (x$convergence == 0) return(x$par) else return(rep(NA, 3))
  }

  ests <- t(sapply(1:1000, est.PPP))
  ests <- ests[complete.cases(ests), ]

  if (nrow(ests) > 50) {
    covmats[[i]] <- cov(ests)
  } else {
    # Fall back to Hessian-based covariance if bootstrap fails
    covmats[[i]] <- ppfit.mles[[i]]$cov
  }
}

# Package Stage 1 results
mles.covmats <- lapply(1:ns, function(i) {
  list(mle = ppfit.mles[[i]]$mle, covmat = covmats[[i]])
})

saveRDS(list(mles.covmats = mles.covmats, loc = loc,
             station_meta = station_meta, thresholds = thresholds,
             mles = mles, mles.original = mles.original,
             nt_per_station = nt_per_station),
        "data/stage1_results.rds")
cat("  Stage 1 results saved to data/stage1_results.rds\n")

cat("\n=== Stage 1 complete ===\n")
cat("Run R/03_stage2_smooth.R for spatial smoothing.\n")
