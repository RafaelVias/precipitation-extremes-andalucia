# 03_stage2_smooth.R — Stage 2: Matérn(5/2) GP with PC priors
#
# Spatial smoothing of Stage 1 MLEs using a Matérn(5/2) Gaussian process
# with penalised complexity priors (Fuglstad et al. 2019), fitted in Stan.
#
# Priors:
#   - mu ~ N(0, 100)
#   - sigma_gp ~ Exp(lambda_s)
#   - phi_gp ~ PC prior: lambda_rho * rho^{-2} * exp(-lambda_rho / rho)
#   - nugget ~ Exp(lambda_nugget)
#
# Run from project root: Rscript R/03_stage2_smooth.R

library(cmdstanr)
library(Matrix)

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("STAGE 2 (Stan): Matern(5/2) GP + PC priors\n")
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

cat("    Distance range:", round(min(dist_mat[upper.tri(dist_mat)]), 3), "-",
    round(max(dist_mat[upper.tri(dist_mat)]), 3), "degrees\n")

# ---- 2. Build block-diagonal precision matrix ----
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

# ---- 3. PC prior rate parameters ----
# PC prior for range: P(rho < rho_0) = alpha  =>  lambda_rho = -rho_0 * log(alpha)
# PC prior for sigma: P(sigma > U) = alpha    =>  lambda_s   = -log(alpha) / U

# sigma_gp: P(sigma > 1) = 0.05  =>  lambda_s = 3.0
lambda_s <- c(3.0, 3.0, 3.0)

# phi_gp (range): P(rho < 0.1 deg) = 0.05  =>  lambda_rho = -0.1 * log(0.05) = 0.30
lambda_rho <- c(0.30, 0.30, 0.30)

# nugget: P(nugget > 0.5) = 0.05  =>  lambda_nugget = -log(0.05)/0.5 = 6.0
lambda_nugget <- c(6.0, 6.0, 6.0)

cat("\n  PC prior rates:\n")
cat(sprintf("    lambda_s     = [%.1f, %.1f, %.1f]  (P(sigma>1)=0.05)\n",
            lambda_s[1], lambda_s[2], lambda_s[3]))
cat(sprintf("    lambda_rho   = [%.2f, %.2f, %.2f]  (P(rho<0.1)=0.05)\n",
            lambda_rho[1], lambda_rho[2], lambda_rho[3]))
cat(sprintf("    lambda_nugget= [%.1f, %.1f, %.1f]  (P(nugget>0.5)=0.05)\n",
            lambda_nugget[1], lambda_nugget[2], lambda_nugget[3]))

# ---- 4. Prepare Stan data ----
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
  log_det_Q        = log_det_Q,
  lambda_s         = lambda_s,
  lambda_rho       = lambda_rho,
  lambda_nugget    = lambda_nugget
)

# ---- 5. Prepare initial values ----
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

# ---- 6. Compile and run Stan model ----
cat("  Compiling Stan model...\n")
model <- cmdstan_model("Stan/smooth_matern_pc.stan")

cat("  Running NUTS (4 chains x 2000 iterations)...\n")
t0 <- proc.time()

fit <- model$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 100,
  adapt_delta     = 0.9,
  max_treedepth   = 12,
  init            = rep(list(inits), 4)
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Stan completed in %.1f seconds (%.1f min)\n", elapsed, elapsed / 60))

# ---- 7. Diagnostics ----
cat("\n  === Stan Diagnostics ===\n")
cat("  Hyperparameter summary:\n")
print(fit$summary(c("mu", "sigma_gp", "phi_gp", "nugget")))

diag_summary <- fit$diagnostic_summary()
cat("\n  Divergences per chain:", diag_summary$num_divergent, "\n")
cat("  Max treedepth per chain:", diag_summary$num_max_treedepth, "\n")

# ---- 8. Convert to output format ----
cat("\n  Converting output format...\n")

draws <- fit$draws(format = "draws_matrix")

psi_cols <- paste0("eta[", 1:ns, ",1]")
tau_cols <- paste0("eta[", 1:ns, ",2]")
phi_cols <- paste0("eta[", 1:ns, ",3]")

psi_draws <- as.matrix(draws[, psi_cols])
tau_draws <- as.matrix(draws[, tau_cols])
phi_draws <- as.matrix(draws[, phi_cols])

mu_draws       <- as.matrix(draws[, c("mu[1]", "mu[2]", "mu[3]")])
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

  beta_psi = mu_draws[, 1],
  beta_tau = mu_draws[, 2],
  beta_phi = mu_draws[, 3],

  sigma_gp_psi = sigma_gp_draws[, 1],
  sigma_gp_tau = sigma_gp_draws[, 2],
  sigma_gp_phi = sigma_gp_draws[, 3],
  phi_gp_psi   = phi_gp_draws[, 1],
  phi_gp_tau   = phi_gp_draws[, 2],
  phi_gp_phi   = phi_gp_draws[, 3],
  nugget_psi   = nugget_draws[, 1],
  nugget_tau   = nugget_draws[, 2],
  nugget_phi   = nugget_draws[, 3],

  method   = "Matern52_GP_PC_Stan",
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
