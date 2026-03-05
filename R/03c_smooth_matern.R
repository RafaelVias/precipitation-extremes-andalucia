# 03c_smooth_matern.R — Stage 2: Matérn(5/2) GP spatial smoothing via Stan/NUTS
#
# Replaces the BYM2 (ICAR + iid) prior from 03b_smooth_stan.R with a proper
# Matérn(5/2) Gaussian Process. Each GEV parameter field (psi, tau, phi) gets
# its own GP hyperparameters: marginal SD (sigma_gp), range (phi_gp), nugget.
#
# The Stage 1 likelihood is unchanged — same block-diagonal precision matrix
# from bootstrap covariances.
#
# Run from project root: Rscript R/03c_smooth_matern.R

library(cmdstanr)
library(Matrix)

# Load Stage 1 transform functions (g, h)
source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("STAGE 2 (Stan): Matern(5/2) GP smoothing\n")
cat("========================================\n")

# ---- Load Stage 1 results ----
s1 <- readRDS("data/stage1_results.rds")
mles.covmats <- s1$mles.covmats
loc          <- s1$loc
meta         <- s1$station_meta
ns           <- nrow(loc)

cat("  Stations:", ns, "\n")

# ---- 1. Compute pairwise distance matrix ----
cat("  Computing pairwise distance matrix...\n")

dist_mat <- as.matrix(dist(loc))

cat("    Dimension:", nrow(dist_mat), "x", ncol(dist_mat), "\n")
cat("    Distance range:", round(min(dist_mat[upper.tri(dist_mat)]), 3), "-",
    round(max(dist_mat[upper.tri(dist_mat)]), 3), "degrees\n")

# ---- 2. Build block-diagonal precision matrix (parameter-major order) ----
cat("  Building precision matrix from bootstrap covariances...\n")

# Stack eta_hat in parameter-major order: [psi_1..psi_P, tau_1..tau_P, phi_1..phi_P]
psi_hat <- sapply(mles.covmats, function(x) x$mle[1])
tau_hat <- sapply(mles.covmats, function(x) x$mle[2])
phi_hat <- sapply(mles.covmats, function(x) x$mle[3])
eta_hat <- c(psi_hat, tau_hat, phi_hat)

# Build the 3P x 3P precision matrix in parameter-major order.
# Each station contributes a 3x3 block Q_i = solve(Sigma_i),
# placed at positions (i + (r-1)*P, i + (c-1)*P) for r,c in 1:3.
n_total <- 3 * ns
rows <- cols <- vals <- c()

n_singular <- 0
for (i in seq_len(ns)) {
  Sigma_i <- mles.covmats[[i]]$covmat

  # Check condition number — use pseudoinverse if near-singular
  ev <- eigen(Sigma_i, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) / max(ev) < 1e-10 || min(ev) <= 0) {
    n_singular <- n_singular + 1
    # Add small ridge to make invertible
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

if (n_singular > 0) {
  cat("    Warning:", n_singular, "stations needed ridge regularization\n")
}

Q_full <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_total, n_total))

# Sparse Cholesky
L <- t(chol(Q_full))

# Extract CSC format for Stan
n_values <- diff(L@p)  # non-zeros per column
index    <- L@i + 1L    # row indices (1-indexed for Stan)
value    <- L@x         # values
log_det_Q <- sum(log(diag(L)))

cat("    Precision matrix:", n_total, "x", n_total, "\n")
cat("    Non-zeros in L:", length(value), "\n")
cat("    log|L|:", round(log_det_Q, 2), "\n")

# ---- 3. Prepare Stan data ----
cat("  Preparing Stan data...\n")

stan_data <- list(
  n_stations       = ns,
  n_param          = 3L,
  eta_hat          = eta_hat,
  dist_mat         = dist_mat,
  n_nonzero_chol_Q = length(value),
  n_values         = as.integer(n_values),
  index            = as.integer(index),
  value            = value,
  log_det_Q        = log_det_Q
)

# ---- 4. Prepare initial values ----
mu_psi <- mean(psi_hat)
mu_tau <- mean(tau_hat)
mu_phi <- mean(phi_hat)
sd_psi <- sd(psi_hat)
sd_tau <- sd(tau_hat)
sd_phi <- sd(phi_hat)

psi_raw <- (psi_hat - mu_psi) / sd_psi
tau_raw <- (tau_hat - mu_tau) / sd_tau
phi_raw <- (phi_hat - mu_phi) / max(sd_phi, 1e-4)

inits <- list(
  mu       = c(mu_psi, mu_tau, mu_phi),
  sigma_gp = c(sd_psi, sd_tau, max(sd_phi, 0.01)),
  phi_gp   = c(0.3, 0.3, 0.5),
  nugget   = c(0.01, 0.01, 0.01),
  eta_raw  = cbind(psi_raw, tau_raw, phi_raw)
)

# ---- 5. Compile and run Stan model ----
cat("  Compiling Stan model...\n")
model <- cmdstan_model("Stan/smooth_matern.stan")

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
cat(sprintf("  Stan completed in %.1f seconds\n", elapsed))

# ---- 6. Diagnostics ----
cat("\n  === Stan Diagnostics ===\n")
cat("  Hyperparameter summary:\n")
print(fit$summary(c("mu", "sigma_gp", "phi_gp", "nugget")))

diag_summary <- fit$diagnostic_summary()
cat("\n  Divergences per chain:", diag_summary$num_divergent, "\n")
cat("  Max treedepth per chain:", diag_summary$num_max_treedepth, "\n")

# ---- 7. Convert to output format ----
cat("\n  Converting output format...\n")

draws <- fit$draws(format = "draws_matrix")

# Extract eta draws: eta[station, param] — param 1=psi, 2=tau, 3=phi
psi_cols <- paste0("eta[", 1:ns, ",1]")
tau_cols <- paste0("eta[", 1:ns, ",2]")
phi_cols <- paste0("eta[", 1:ns, ",3]")

psi_draws <- as.matrix(draws[, psi_cols])  # [n_draws x ns]
tau_draws <- as.matrix(draws[, tau_cols])
phi_draws <- as.matrix(draws[, phi_cols])

# Hyperparameter draws
mu_draws       <- as.matrix(draws[, c("mu[1]", "mu[2]", "mu[3]")])
sigma_gp_draws <- as.matrix(draws[, c("sigma_gp[1]", "sigma_gp[2]", "sigma_gp[3]")])
phi_gp_draws   <- as.matrix(draws[, c("phi_gp[1]", "phi_gp[2]", "phi_gp[3]")])
nugget_draws   <- as.matrix(draws[, c("nugget[1]", "nugget[2]", "nugget[3]")])

# Build output list
result <- list(
  # Parameter chains [n_draws x ns]
  psi.selected = psi_draws,
  tau.selected = tau_draws,
  phi.selected = phi_draws,

  # Posterior summaries
  psi.posmean = colMeans(psi_draws),
  tau.posmean = colMeans(tau_draws),
  phi.posmean = colMeans(phi_draws),
  psi.posvar  = apply(psi_draws, 2, var),
  tau.posvar  = apply(tau_draws, 2, var),
  phi.posvar  = apply(phi_draws, 2, var),

  # Mean hyperparameter chains (compatible with 04_update_figures.R)
  beta_psi = mu_draws[, 1],
  beta_tau = mu_draws[, 2],
  beta_phi = mu_draws[, 3],

  # GP hyperparameter chains
  sigma_gp_psi = sigma_gp_draws[, 1],
  sigma_gp_tau = sigma_gp_draws[, 2],
  sigma_gp_phi = sigma_gp_draws[, 3],
  phi_gp_psi   = phi_gp_draws[, 1],
  phi_gp_tau   = phi_gp_draws[, 2],
  phi_gp_phi   = phi_gp_draws[, 3],
  nugget_psi   = nugget_draws[, 1],
  nugget_tau   = nugget_draws[, 2],
  nugget_phi   = nugget_draws[, 3],

  # Metadata
  method   = "Matern52_GP_Stan",
  minutes  = elapsed / 60,
  stan_fit = fit
)

saveRDS(result, "data/stage2_matern_results.rds")
cat("  Results saved to data/stage2_matern_results.rds\n")

cat("\n========================================\n")
cat("Done. Matern GP smooth step complete.\n")
cat("========================================\n")
