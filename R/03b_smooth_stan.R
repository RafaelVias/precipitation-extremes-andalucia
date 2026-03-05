# 03b_smooth_stan.R — Stage 2: BYM2 spatial smoothing via Stan/NUTS
#
# Replaces the custom MCMC SPDE sampler (03_max_and_smooth.R Stage 2) with
# Brynjolfur Gauti's BYM2 + Stan approach. Uses the same Stage 1 MLEs.
#
# Run from project root: Rscript R/03b_smooth_stan.R

library(cmdstanr)
library(Matrix)
library(deldir)
library(INLA)

# Load Stage 1 transform functions (g, h)
source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("STAGE 2 (Stan): BYM2 spatial smoothing\n")
cat("========================================\n")

# ---- Load Stage 1 results ----
s1 <- readRDS("data/stage1_results.rds")
mles.covmats <- s1$mles.covmats
loc          <- s1$loc
meta         <- s1$station_meta
ns           <- nrow(loc)

cat("  Stations:", ns, "\n")

# ---- 1. Build neighborhood graph via Delaunay triangulation ----
cat("  Building Delaunay neighborhood graph...\n")

tri <- deldir(loc[, "lon"], loc[, "lat"])

# Extract edges from Delaunay triangulation
# deldir returns edges in $delsgs with columns (x1, y1, x2, y2, ind1, ind2)
del_edges <- tri$delsgs[, c("ind1", "ind2")]

# Ensure both directions are present and node1 < node2 for ICAR
edges_raw <- rbind(
  data.frame(station = del_edges$ind1, neighbor = del_edges$ind2),
  data.frame(station = del_edges$ind2, neighbor = del_edges$ind1)
)

# Keep only edges where station < neighbor (for ICAR prior)
edges_upper <- edges_raw[edges_raw$station < edges_raw$neighbor, ]
edges_upper <- unique(edges_upper)

# But the Stan model needs ALL directed edges (both directions)
edges_all <- rbind(
  edges_upper,
  data.frame(station = edges_upper$neighbor, neighbor = edges_upper$station)
)

cat("    Edges:", nrow(edges_upper), "(undirected),",
    nrow(edges_all), "(directed)\n")
cat("    Mean neighbors per station:",
    round(nrow(edges_all) / ns, 1), "\n")

# ---- 2. Compute BYM2 scaling factor ----
cat("  Computing BYM2 scaling factor...\n")

get_scaling_factor <- function(edges_upper, N) {
  # Build adjacency matrix
  adj_mat <- sparseMatrix(
    i = edges_upper$station,
    j = edges_upper$neighbor,
    x = 1,
    dims = c(N, N),
    symmetric = TRUE
  )

  # ICAR precision matrix (singular — graph Laplacian)
  Q <- Diagonal(N, rowSums(adj_mat)) - adj_mat

  # Small jitter for numerical stability
  Q_pert <- Q + Diagonal(N) * max(diag(Q)) * sqrt(.Machine$double.eps)

  # Constrained inverse (sum-to-zero)
  Q_inv <- inla.qinv(Q_pert, constr = list(A = matrix(1, 1, N), e = 0))

  # Geometric mean of marginal variances
  exp(mean(log(diag(Q_inv))))
}

scaling_factor <- get_scaling_factor(edges_upper, ns)
cat("    Scaling factor:", round(scaling_factor, 4), "\n")

# ---- 3. Build block-diagonal precision matrix (parameter-major order) ----
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

# ---- 4. Prepare Stan data ----
cat("  Preparing Stan data...\n")

stan_data <- list(
  n_stations       = ns,
  n_obs            = ns,  # placeholder (not used in likelihood)
  n_param          = 3L,
  eta_hat          = eta_hat,
  n_edges          = nrow(edges_all),
  node1            = edges_all$station,
  node2            = edges_all$neighbor,
  scaling_factor   = scaling_factor,
  n_nonzero_chol_Q = length(value),
  n_values         = as.integer(n_values),
  index            = as.integer(index),
  value            = value,
  log_det_Q        = log_det_Q
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
phi_raw <- (phi_hat - mu_phi) / max(sd_phi, 1e-4)  # guard against zero SD

eta_raw <- cbind(psi_raw, tau_raw, phi_raw)

inits <- list(
  mu          = c(mu_psi, mu_tau, mu_phi),
  sigma       = c(sd_psi, sd_tau, max(sd_phi, 0.01)),
  rho         = c(0.5, 0.5, 0.5),
  eta_spatial = eta_raw,
  eta_random  = matrix(0, nrow = ns, ncol = 3)
)

# ---- 6. Compile and run Stan model ----
cat("  Compiling Stan model...\n")
model <- cmdstan_model("Stan/smooth_bym2.stan")

cat("  Running NUTS (4 chains x 2000 iterations)...\n")
t0 <- proc.time()

fit <- model$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  refresh         = 200,
  init            = rep(list(inits), 4)
)

elapsed <- (proc.time() - t0)[3]
cat(sprintf("  Stan completed in %.1f seconds\n", elapsed))

# ---- 7. Diagnostics ----
cat("\n  === Stan Diagnostics ===\n")
cat("  Hyperparameter summary:\n")
print(fit$summary(c("mu", "sigma", "rho")))

diag_summary <- fit$diagnostic_summary()
cat("\n  Divergences per chain:", diag_summary$num_divergent, "\n")
cat("  Max treedepth per chain:", diag_summary$num_max_treedepth, "\n")

# ---- 8. Convert to format compatible with 04_update_figures.R ----
cat("\n  Converting output format...\n")

draws <- fit$draws(format = "draws_matrix")

# Extract eta draws: eta[station, param] — param 1=psi, 2=tau, 3=phi
# Stan returns eta[i,j] as "eta[i,j]" columns
psi_cols <- paste0("eta[", 1:ns, ",1]")
tau_cols <- paste0("eta[", 1:ns, ",2]")
phi_cols <- paste0("eta[", 1:ns, ",3]")

psi_draws <- as.matrix(draws[, psi_cols])  # [n_draws x ns]
tau_draws <- as.matrix(draws[, tau_cols])
phi_draws <- as.matrix(draws[, phi_cols])

# Hyperparameter draws
mu_draws    <- as.matrix(draws[, c("mu[1]", "mu[2]", "mu[3]")])
sigma_draws <- as.matrix(draws[, c("sigma[1]", "sigma[2]", "sigma[3]")])
rho_draws   <- as.matrix(draws[, c("rho[1]", "rho[2]", "rho[3]")])

# Build output list matching 04_update_figures.R expectations
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

  # Hyperparameter chains (for trace plots)
  beta_psi  = mu_draws[, 1],
  beta_tau  = mu_draws[, 2],
  beta_phi  = mu_draws[, 3],
  sigma_psi = sigma_draws[, 1],
  sigma_tau = sigma_draws[, 2],
  sigma_phi = sigma_draws[, 3],
  rho_psi   = rho_draws[, 1],
  rho_tau   = rho_draws[, 2],
  rho_phi   = rho_draws[, 3],

  # Metadata
  method  = "BYM2_Stan",
  minutes = elapsed / 60,

  # Store the full Stan fit for detailed diagnostics
  stan_fit = fit
)

saveRDS(result, "data/stage2_stan_results.rds")
cat("  Results saved to data/stage2_stan_results.rds\n")

cat("\n========================================\n")
cat("Done. Stan BYM2 smooth step complete.\n")
cat("========================================\n")
