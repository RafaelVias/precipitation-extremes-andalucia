# 03_stage2_smooth.R — Stage 2: Matern(5/2) GP with PC priors + covariates
#
# Spatial smoothing of Stage 1 MLEs using a Matern(5/2) Gaussian process
# with penalised complexity priors (Fuglstad et al. 2019), fitted in Stan.
#
# Model structure (following Hazra, Huser & Johannesson, BLGM 2023 Ch.7):
#   eta_psi(s) = X(s)' beta_psi + f_psi(s)   (covariates + GP)
#   eta_tau(s) = mu_tau + f_tau(s)             (intercept + GP)
#   eta_phi(s) = mu_phi + f_phi(s)             (intercept + GP)
#
# All three parameters receive Matern(5/2) GPs. The phi GP uses tighter
# PC priors reflecting the smaller spatial variation in the shape parameter.
#
# Covariates on psi: intercept, altitude (std), exposure (std), alt x exposure
#
# Run from project root: Rscript R/03_stage2_smooth.R

library(cmdstanr)
library(Matrix)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("STAGE 2: Matern(5/2) GP + PC priors + covariates\n")
cat("  All 3 parameters get GPs (phi with tighter priors)\n")
cat("========================================\n")

# ---- Load Stage 1 results ----
s1 <- readRDS("data/stage1_results.rds")
mles.covmats <- s1$mles.covmats
loc          <- s1$loc
meta         <- s1$station_meta
ns           <- nrow(loc)

cat("  Stations:", ns, "\n")

# ---- 1. Load covariates ----
cat("  Loading covariates...\n")

stn_exp <- readRDS("data/station_exposure.rds")

# Match station order to Stage 1 results
s1_ids <- meta$indicativo
exp_idx <- match(s1_ids, stn_exp$indicativo)

if (any(is.na(exp_idx))) {
  stop("Station mismatch: ", sum(is.na(exp_idx)), " Stage 1 stations not found in exposure data")
}

alt_raw <- stn_exp$alt_dem[exp_idx]
exp_raw <- stn_exp$exposure_mean[exp_idx]

# Standardise covariates (save params for prediction)
alt_mean <- mean(alt_raw)
alt_sd   <- sd(alt_raw)
exp_mean <- mean(exp_raw)
exp_sd   <- sd(exp_raw)

alt_std <- (alt_raw - alt_mean) / alt_sd
exp_std <- (exp_raw - exp_mean) / exp_sd

# Design matrix: intercept, alt, exposure, alt x exposure
X_psi <- cbind(1, alt_std, exp_std, alt_std * exp_std)
n_cov_psi <- ncol(X_psi)

cat(sprintf("  Covariates: altitude (mean=%.0f, sd=%.0f), exposure (mean=%.0f, sd=%.0f)\n",
            alt_mean, alt_sd, exp_mean, exp_sd))
cat(sprintf("  Design matrix: %d x %d (intercept + alt + exposure + alt x exposure)\n",
            nrow(X_psi), n_cov_psi))

# ---- 2. Compute pairwise distance matrix ----
cat("  Computing pairwise distance matrix...\n")
dist_mat <- as.matrix(dist(loc))

cat("    Distance range:", round(min(dist_mat[upper.tri(dist_mat)]), 3), "-",
    round(max(dist_mat[upper.tri(dist_mat)]), 3), "degrees\n")

# ---- 3. Build block-diagonal precision matrix ----
cat("  Building precision matrix from bootstrap covariances...\n")

psi_hat <- sapply(mles.covmats, function(x) x$mle[1])
tau_hat <- sapply(mles.covmats, function(x) x$mle[2])
phi_hat <- sapply(mles.covmats, function(x) x$mle[3])
eta_hat <- c(psi_hat, tau_hat, phi_hat)

n_total <- 3 * ns
rows <- cols <- vals <- c()
n_singular <- 0

for (i in seq_len(ns)) {
  Sigma_i <- mles.covmats[[i]]$covmat
  ev <- eigen(Sigma_i, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) / max(ev) < 1e-10 || min(ev) <= 0) {
    n_singular <- n_singular + 1
    Sigma_i <- Sigma_i + diag(3) * max(ev) * 1e-6
  }
  Q_i <- solve(Sigma_i)
  for (r in 1:3) {
    for (cc in 1:3) {
      rows <- c(rows, i + (r - 1) * ns)
      cols <- c(cols, i + (cc - 1) * ns)
      vals <- c(vals, Q_i[r, cc])
    }
  }
}

if (n_singular > 0) cat("    Warning:", n_singular, "stations needed ridge\n")

Q_full <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_total, n_total))
L <- t(chol(Q_full))
n_values <- diff(L@p)
index    <- L@i + 1L
value    <- L@x
log_det_Q <- sum(log(diag(L)))

cat("    Precision matrix:", n_total, "x", n_total, "\n")
cat("    Non-zeros in L:", length(value), "\n")

# ---- 4. PC prior rate parameters ----
# All 3 GPs: psi, tau, phi (phi gets tighter priors)
n_gp <- 3L

# sigma_gp: psi,tau: P(sigma > 1) = 0.05 => lambda = 3.0
#           phi:     P(sigma > 0.3) = 0.05 => lambda = 10.0
lambda_s <- c(3.0, 3.0, 10.0)

# phi_gp (range): psi,tau: P(rho < 0.1 deg) = 0.05 => lambda = 0.30
#                 phi:     P(rho < 0.5 deg) = 0.05 => lambda = 1.50
lambda_rho <- c(0.30, 0.30, 1.50)

# nugget: psi,tau: P(nugget > 0.5) = 0.05 => lambda = 6.0
#         phi:     P(nugget > 0.3) = 0.05 => lambda = 10.0
lambda_nugget <- c(6.0, 6.0, 10.0)

cat("\n  PC prior rates (psi, tau, phi GPs):\n")
cat(sprintf("    lambda_s     = [%.1f, %.1f, %.1f]\n",
            lambda_s[1], lambda_s[2], lambda_s[3]))
cat(sprintf("    lambda_rho   = [%.2f, %.2f, %.2f]\n",
            lambda_rho[1], lambda_rho[2], lambda_rho[3]))
cat(sprintf("    lambda_nugget= [%.1f, %.1f, %.1f]\n",
            lambda_nugget[1], lambda_nugget[2], lambda_nugget[3]))

# ---- 5. Prepare Stan data ----
cat("  Preparing Stan data...\n")

stan_data <- list(
  n_stations       = ns,
  n_param          = 3L,
  n_gp             = n_gp,
  eta_hat          = eta_hat,
  dist_mat         = dist_mat,
  n_nonzero_chol_Q = length(value),
  n_values         = as.integer(n_values),
  index            = as.integer(index),
  value            = value,
  log_det_Q        = log_det_Q,
  n_cov_psi        = n_cov_psi,
  X_psi            = X_psi,
  prior_sd_beta    = 10.0,
  lambda_s         = lambda_s,
  lambda_rho       = lambda_rho,
  lambda_nugget    = lambda_nugget
)

# ---- 6. Prepare initial values ----
# OLS estimate for beta_psi as starting point
beta_psi_init <- as.numeric(solve(t(X_psi) %*% X_psi, t(X_psi) %*% psi_hat))

mu_tau <- mean(tau_hat)
mu_phi <- mean(phi_hat)
sd_psi <- sd(psi_hat - X_psi %*% beta_psi_init)
sd_tau <- sd(tau_hat)
sd_phi <- sd(phi_hat)

psi_resid <- psi_hat - X_psi %*% beta_psi_init
psi_raw <- as.numeric(psi_resid / sd_psi)
tau_raw <- (tau_hat - mu_tau) / sd_tau
phi_raw <- (phi_hat - mu_phi) / max(sd_phi, 1e-4)

inits <- list(
  beta_psi   = beta_psi_init,
  mu_tau     = mu_tau,
  mu_phi     = mu_phi,
  sigma_gp   = c(sd_psi, sd_tau, max(sd_phi, 0.01)),
  phi_gp     = c(0.3, 0.3, 0.5),
  nugget     = c(0.01, 0.01, 0.01),
  eta_raw    = cbind(psi_raw, tau_raw, phi_raw)
)

cat(sprintf("  OLS init for beta_psi: [%.3f, %.3f, %.3f, %.3f]\n",
            beta_psi_init[1], beta_psi_init[2], beta_psi_init[3], beta_psi_init[4]))

# ---- 7. Compile and run Stan model ----
cat("  Compiling Stan model...\n")
model <- cmdstan_model("Stan/smooth_matern_pc_cov.stan")

cat("  Running NUTS (4 chains x 2000 iterations)...\n")
t0 <- proc.time()

fit <- model$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 200,
  adapt_delta     = 0.9,
  max_treedepth   = 12,
  init            = rep(list(inits), 4)
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Stan completed in %.1f seconds (%.1f min)\n", elapsed, elapsed / 60))

# ---- 8. Diagnostics ----
cat("\n  === Stan Diagnostics ===\n")
cat("  Beta_psi (covariate effects on GEV location):\n")
print(fit$summary(c("beta_psi")))
cat("\n  Intercepts and GP hyperparameters:\n")
print(fit$summary(c("mu_tau", "mu_phi", "sigma_gp", "phi_gp", "nugget")))

diag_summary <- fit$diagnostic_summary()
cat("\n  Divergences per chain:", diag_summary$num_divergent, "\n")
cat("  Max treedepth per chain:", diag_summary$num_max_treedepth, "\n")

# ---- 9. Convert to output format ----
cat("\n  Converting output format...\n")

draws <- fit$draws(format = "draws_matrix")

psi_cols <- paste0("eta[", 1:ns, ",1]")
tau_cols <- paste0("eta[", 1:ns, ",2]")
phi_cols <- paste0("eta[", 1:ns, ",3]")

psi_draws <- as.matrix(draws[, psi_cols])
tau_draws <- as.matrix(draws[, tau_cols])
phi_draws <- as.matrix(draws[, phi_cols])

beta_psi_cols <- paste0("beta_psi[", 1:n_cov_psi, "]")
beta_psi_draws <- as.matrix(draws[, beta_psi_cols])

mu_tau_draws <- as.numeric(draws[, "mu_tau"])
mu_phi_draws <- as.numeric(draws[, "mu_phi"])

sigma_gp_draws <- as.matrix(draws[, c("sigma_gp[1]", "sigma_gp[2]", "sigma_gp[3]")])
phi_gp_draws   <- as.matrix(draws[, c("phi_gp[1]", "phi_gp[2]", "phi_gp[3]")])
nugget_draws   <- as.matrix(draws[, c("nugget[1]", "nugget[2]", "nugget[3]")])

result <- list(
  psi.selected = psi_draws,
  tau.selected = tau_draws,
  phi.selected = phi_draws,

  psi.posmean = colMeans(psi_draws),
  tau.posmean = colMeans(tau_draws),
  phi.posmean = colMeans(phi_draws),
  psi.posvar  = apply(psi_draws, 2, var),
  tau.posvar  = apply(tau_draws, 2, var),
  phi.posvar  = apply(phi_draws, 2, var),

  # Covariate coefficients for psi (n_draws x n_cov_psi matrix)
  beta_psi_draws = beta_psi_draws,
  # Intercepts for tau, phi (vectors of length n_draws)
  beta_tau = mu_tau_draws,
  beta_phi = mu_phi_draws,

  # GP hyperparameters (psi, tau, phi)
  sigma_gp_psi = sigma_gp_draws[, 1],
  sigma_gp_tau = sigma_gp_draws[, 2],
  sigma_gp_phi = sigma_gp_draws[, 3],
  phi_gp_psi   = phi_gp_draws[, 1],
  phi_gp_tau   = phi_gp_draws[, 2],
  phi_gp_phi   = phi_gp_draws[, 3],
  nugget_psi   = nugget_draws[, 1],
  nugget_tau   = nugget_draws[, 2],
  nugget_phi   = nugget_draws[, 3],

  # Flag: phi has spatial GP
  phi_has_gp = TRUE,

  # Covariate standardisation parameters (needed for prediction)
  cov_standardisation = list(
    alt_mean = alt_mean, alt_sd = alt_sd,
    exp_mean = exp_mean, exp_sd = exp_sd
  ),

  method   = "Matern52_GP_PC_cov_Stan",
  minutes  = elapsed / 60,
  pc_prior_rates = list(lambda_s = lambda_s, lambda_rho = lambda_rho,
                        lambda_nugget = lambda_nugget),
  stan_fit = fit
)

saveRDS(result, "data/stage2_matern_pc_results.rds")
cat("  Results saved to data/stage2_matern_pc_results.rds\n")

cat("\n========================================\n")
cat("Done. Stage 2 complete.\n")
cat("========================================\n")
