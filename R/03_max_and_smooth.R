# =============================================================================
# 03_max_and_smooth.R
# Spatial extreme value model using Max-and-Smooth (Hazra et al. 2023)
# Poisson point process likelihood with SPDE spatial random effects
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
#   - vendor/max_and_smooth/stage2_functions.R
#
# Outputs:
#   - data/daily_matrix.rds           Station x day precipitation matrix
#   - data/stage1_results.rds         Site-wise MLEs and covariance matrices
#   - data/stage2_mcmc_results.rds    MCMC posterior chains
#   - figures/mesh_andalucia.pdf      SPDE triangular mesh
#   - figures/mle_maps.pdf            Stage 1 MLE parameter maps
#   - figures/return_level_maps.pdf   Posterior return level maps
# =============================================================================

library(dplyr)
library(lubridate)
library(INLA)
library(Matrix)
library(spam)
library(fields)
library(sf)
library(ggplot2)

# ---- Source the M&S functions ------------------------------------------------
source("vendor/max_and_smooth/stage1_functions.R")
source("vendor/max_and_smooth/stage2_functions.R")

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

saveRDS(list(Y = Y, loc = loc, station_meta = station_meta,
             all_dates = all_dates, nt_per_station = nt_per_station),
        "data/daily_matrix.rds")

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

# =============================================================================
# STAGE 1 DIAGNOSTIC PLOTS
# =============================================================================
cat("\nPlotting Stage 1 diagnostics...\n")

plot_df <- data.frame(station_meta, mles.original, mles)

# MLE parameter maps
station_sf <- st_as_sf(plot_df, coords = c("longitud", "latitud"), crs = 4326)

p_mu <- ggplot(station_sf) +
  geom_sf(aes(colour = mu), size = 2) +
  scale_colour_viridis_c(option = "C") +
  labs(title = expression(hat(mu) ~ "(location)"), colour = "mm") +
  theme_minimal()

p_sigma <- ggplot(station_sf) +
  geom_sf(aes(colour = sigma), size = 2) +
  scale_colour_viridis_c(option = "C") +
  labs(title = expression(hat(sigma) ~ "(scale)"), colour = "mm") +
  theme_minimal()

p_xi <- ggplot(station_sf) +
  geom_sf(aes(colour = xi), size = 2) +
  scale_colour_viridis_c(option = "D") +
  labs(title = expression(hat(xi) ~ "(shape)"), colour = "") +
  theme_minimal()

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_mle <- p_mu + p_sigma + p_xi +
    plot_annotation(title = "Stage 1: Site-wise Poisson point process MLEs")
  ggsave("figures/mle_maps.pdf", p_mle, width = 14, height = 5)
}

# =============================================================================
# STAGE 2: MCMC Smoothing with SPDE spatial effects
# =============================================================================
cat("\n========================================\n")
cat("STAGE 2: Building SPDE mesh and running MCMC...\n")
cat("========================================\n")

# Build INLA mesh over Andalucia
# Using station locations with boundary buffer
cat("  Building triangular mesh...\n")

mesh.domain <- inla.mesh.2d(
  loc = loc,
  max.edge = c(0.3, 1.0),   # inner/outer max triangle edge (degrees)
  offset = c(0.2, 0.5),     # inner/outer buffer
  cutoff = 0.15             # min distance between points
)

cat("  Mesh: ", mesh.domain$n, "nodes,", nrow(mesh.domain$graph$tv), "triangles\n")

# Save mesh plot
pdf("figures/mesh_andalucia.pdf", width = 10, height = 8)
plot(mesh.domain, main = sprintf("SPDE mesh for Andalucia (%d nodes)", mesh.domain$n))
points(loc, pch = 16, col = "red", cex = 0.8)
dev.off()

# Build INLA matrices
A <- inla.spde.make.A(mesh = mesh.domain, loc = loc)
fem.mesh <- inla.mesh.fem(mesh.domain, order = 2)
inla.mats <- list(
  c.mat = fem.mesh$c0,
  g1.mat = fem.mesh$g1,
  g2.mat = fem.mesh$g2,
  A = A
)

# ---- Fix xi to a constant across all stations ----
# The data cannot support spatially varying xi (SNR = 1.2, and
# cor(tau, phi) ~ 0.98 in the bootstrap covariances creates
# near-collinearity that cripples the MCMC).
# Fix xi = 0.1 (the MLE median/mean) and pin phi in the surrogate.
xi_fixed <- 0.1
phi_fixed <- h(xi_fixed)
cat(sprintf("  Fixing xi = %.2f (phi = %.4f) for all stations\n", xi_fixed, phi_fixed))

for (i in seq_along(mles.covmats)) {
  mles.covmats[[i]]$mle[3] <- phi_fixed
  # Zero out cross-correlations with phi but keep phi's own variance
  # (it's the cor(tau,phi)~0.98 that cripples the sampler, not the variance)
  orig_phi_var <- mles.covmats[[i]]$covmat[3, 3]
  mles.covmats[[i]]$covmat[3, ] <- 0
  mles.covmats[[i]]$covmat[, 3] <- 0
  mles.covmats[[i]]$covmat[3, 3] <- orig_phi_var
}

# Run MCMC
cat("  Starting MCMC (1,000 iterations, 200 burn-in, thin=1)...\n")
cat("  Expected time: ~1 minute...\n")

fit.mcmc <- mcmc.maxNsmooth(
  mles.covmats = mles.covmats,
  loc = loc,
  inla.mats = inla.mats,
  alpha = 2,
  # Fix phi: provide non-zero init since var(phi)=0 when all phi MLEs are identical
  sigmaSq_phi.init = 1e-4,
  # priors (PC priors via exponential rates)
  sd_beta_psi = 1e2,
  sd_beta_tau = 1e2,
  sd_beta_phi = 1e2,
  lambda_sigma_psi = 0.1,
  lambda_sigma_tau = 0.1,
  lambda_sigma_phi = 0.1,
  lambda_s_psi = 0.1,
  lambda_s_tau = 0.1,
  lambda_rho_psi = 0.1,
  lambda_rho_tau = 0.1,
  # MCMC settings (quick diagnostic run)
  iters = 1000,
  burn = 200,
  thin = 1
)

cat("  MCMC completed in", round(fit.mcmc$minutes, 1), "minutes\n")

saveRDS(fit.mcmc, "data/stage2_mcmc_results_v2.rds")
cat("  Stage 2 results saved to data/stage2_mcmc_results_v2.rds\n")

# =============================================================================
# RETURN LEVEL MAPS
# =============================================================================
cat("\nComputing return level maps...\n")

# Transform MCMC chains back to original parameters
mu.chain <- exp(fit.mcmc$psi.selected)
sigma.chain <- exp(fit.mcmc$psi.selected + fit.mcmc$tau.selected)
xi.chain <- g(fit.mcmc$phi.selected)

# Return levels for M = 20, 50, 100 years
for (M in c(20, 50, 100)) {
  rl.chain <- mu.chain + sigma.chain / xi.chain *
    ((-log(1 - 1/M))^(-xi.chain) - 1)
  rl.posmean <- apply(rl.chain, 2, mean)
  rl.possd <- apply(rl.chain, 2, sd)

  assign(paste0("rl", M, ".mean"), rl.posmean)
  assign(paste0("rl", M, ".sd"), rl.possd)

  cat(sprintf("  %d-year return level: range %.1f - %.1f mm (mean), SD %.1f - %.1f\n",
              M, min(rl.posmean), max(rl.posmean), min(rl.possd), max(rl.possd)))
}

# Smoothed parameter posterior means
mu.posmean <- exp(fit.mcmc$psi.posmean)
sigma.posmean <- exp(fit.mcmc$psi.posmean + fit.mcmc$tau.posmean)
xi.posmean <- g(fit.mcmc$phi.posmean)

cat("\n  Smoothed parameter summaries:\n")
cat(sprintf("    mu:    %.1f - %.1f\n", min(mu.posmean), max(mu.posmean)))
cat(sprintf("    sigma: %.1f - %.1f\n", min(sigma.posmean), max(sigma.posmean)))
cat(sprintf("    xi:    %.3f - %.3f\n", min(xi.posmean), max(xi.posmean)))

# Plot return level maps
rl_df <- data.frame(
  station_meta,
  mu_smooth = mu.posmean,
  sigma_smooth = sigma.posmean,
  xi_smooth = xi.posmean,
  rl20_mean = rl20.mean, rl20_sd = rl20.sd,
  rl50_mean = rl50.mean, rl50_sd = rl50.sd,
  rl100_mean = rl100.mean, rl100_sd = rl100.sd
)

rl_sf <- st_as_sf(rl_df, coords = c("longitud", "latitud"), crs = 4326)

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)

  p20 <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl20_mean), size = 2.5) +
    scale_colour_viridis_c(option = "B") +
    labs(title = "20-year return level", colour = "mm") +
    theme_minimal()

  p50 <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl50_mean), size = 2.5) +
    scale_colour_viridis_c(option = "B") +
    labs(title = "50-year return level", colour = "mm") +
    theme_minimal()

  p100 <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl100_mean), size = 2.5) +
    scale_colour_viridis_c(option = "B") +
    labs(title = "100-year return level", colour = "mm") +
    theme_minimal()

  p20sd <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl20_sd), size = 2.5) +
    scale_colour_viridis_c(option = "E") +
    labs(title = "20-year SD", colour = "mm") +
    theme_minimal()

  p50sd <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl50_sd), size = 2.5) +
    scale_colour_viridis_c(option = "E") +
    labs(title = "50-year SD", colour = "mm") +
    theme_minimal()

  p100sd <- ggplot(rl_sf) +
    geom_sf(aes(colour = rl100_sd), size = 2.5) +
    scale_colour_viridis_c(option = "E") +
    labs(title = "100-year SD", colour = "mm") +
    theme_minimal()

  p_rl <- (p20 + p50 + p100) / (p20sd + p50sd + p100sd) +
    plot_annotation(
      title = "Max-and-Smooth: Return level estimates (top) and posterior SD (bottom)",
      subtitle = sprintf("%d stations in Andalucia, Poisson point process + SPDE", ns)
    )
  ggsave("figures/return_level_maps.pdf", p_rl, width = 16, height = 10)
}

# Save results table
saveRDS(rl_df, "data/return_levels_andalucia.rds")
write.csv(rl_df, "data/return_levels_andalucia.csv", row.names = FALSE)

# Hyperparameter summaries
cat("\n=== Hyperparameter posterior summaries ===\n")
hyper_names <- c("beta_psi", "beta_tau", "beta_phi",
                 "sigma_psi", "sigma_tau", "sigma_phi",
                 "s_psi", "s_tau", "rho_psi", "rho_tau")
for (h in hyper_names) {
  chain <- fit.mcmc[[h]]
  cat(sprintf("  %s: mean=%.4f, sd=%.4f, 95%%CI=(%.4f, %.4f)\n",
              h, mean(chain), sd(chain),
              quantile(chain, 0.025), quantile(chain, 0.975)))
}

# ESS diagnostics
cat("\n=== Effective sample sizes ===\n")
library(coda)
for (h in hyper_names) {
  chain <- fit.mcmc[[h]]
  ess <- effectiveSize(chain)
  cat(sprintf("  %-15s ESS = %6.0f / %d\n", h, ess, length(chain)))
}

cat("\n=== Max-and-Smooth analysis complete ===\n")
