// Matern(5/2) GP spatial smoothing with PC priors and covariates on psi
//
// Model structure (following Hazra, Huser & Johannesson, BLGM 2023 Ch.7):
//   eta_psi(s) = X(s)' beta_psi + f_psi(s)   (covariates + GP)
//   eta_tau(s) = mu_tau + f_tau(s)             (intercept + GP)
//   eta_phi(s) = mu_phi + nugget_phi * z(s)    (intercept + iid noise, NO GP)
//
// The shape parameter phi has no spatial GP because empirical variograms
// show negligible spatial structure (50-75% unstructured noise).
// This eliminates divergences from the phi GP and is more parsimonious.
//
// Covariates on psi: intercept, altitude (std), exposure (std), alt x exposure
//
// Priors:
//   - beta_psi ~ N(0, prior_sd_beta)
//   - mu_tau, mu_phi ~ N(0, 100)
//   - sigma_gp[1:2] ~ Exp(lambda_s)     (psi, tau only)
//   - phi_gp[1:2] ~ PC prior            (psi, tau only)
//   - nugget[1:2] ~ Exp(lambda_nugget)   (psi, tau only)
//   - nugget_phi ~ Exp(lambda_nugget_phi)

functions {
  real normal_prec_chol_lpdf(vector y, vector x, array[] int n_values, array[] int index, vector values, real log_det) {
    int N = num_elements(x);
    int counter = 1;
    vector[N] q = rep_vector(0, N);

    for (i in 1:N) {
      for (j in 1:n_values[i]) {
        q[i] += values[counter] * (y[index[counter]] - x[index[counter]]);
        counter += 1;
      }
    }

    return log_det - dot_self(q) / 2;
  }
}

data {
  int<lower = 1> n_stations;
  int<lower = 1> n_param;           // still 3 (psi, tau, phi)
  int<lower = 1> n_gp;              // 2 (only psi and tau get GPs)

  vector[n_stations * n_param] eta_hat;

  matrix[n_stations, n_stations] dist_mat;

  int<lower = 1> n_nonzero_chol_Q;
  array[n_param * n_stations] int n_values;
  array[n_nonzero_chol_Q] int index;
  vector[n_nonzero_chol_Q] value;
  real<lower = 0> log_det_Q;

  // Covariates for psi
  int<lower = 1> n_cov_psi;
  matrix[n_stations, n_cov_psi] X_psi;
  real<lower = 0> prior_sd_beta;

  // PC prior rate parameters (for psi and tau GPs only)
  vector<lower = 0>[n_gp] lambda_s;
  vector<lower = 0>[n_gp] lambda_rho;
  vector<lower = 0>[n_gp] lambda_nugget;
  real<lower = 0> lambda_nugget_phi;   // separate nugget rate for phi
}

parameters {
  vector[n_cov_psi] beta_psi;
  real mu_tau;
  real mu_phi;
  vector<lower = 0>[n_gp] sigma_gp;    // GP marginal SD (psi, tau)
  vector<lower = 0>[n_gp] phi_gp;      // GP range (psi, tau)
  vector<lower = 0>[n_gp] nugget;      // GP nugget (psi, tau)
  real<lower = 0> nugget_phi;           // iid noise SD for phi
  matrix[n_stations, n_gp] eta_raw_gp; // non-centered GP draws (psi, tau)
  vector[n_stations] phi_raw;           // non-centered iid draws for phi
}

transformed parameters {
  matrix[n_stations, n_param] eta;

  // --- psi and tau: covariates/intercept + Matern GP ---
  for (p in 1:n_gp) {
    matrix[n_stations, n_stations] Sigma;
    real sig2 = square(sigma_gp[p]);
    real nug2 = square(nugget[p]);

    for (i in 1:n_stations) {
      Sigma[i, i] = sig2 + nug2;
      for (j in (i + 1):n_stations) {
        real s5 = sqrt(5.0) * dist_mat[i, j] / phi_gp[p];
        real cov_val = sig2 * (1.0 + s5 + square(s5) / 3.0) * exp(-s5);
        Sigma[i, j] = cov_val;
        Sigma[j, i] = cov_val;
      }
    }

    if (p == 1)
      eta[, p] = X_psi * beta_psi + cholesky_decompose(Sigma) * eta_raw_gp[, p];
    else
      eta[, p] = mu_tau + cholesky_decompose(Sigma) * eta_raw_gp[, p];
  }

  // --- phi: intercept + iid noise (no spatial GP) ---
  eta[, 3] = mu_phi + nugget_phi * phi_raw;
}

model {
  vector[n_param * n_stations] eta_vec;
  for (p in 1:n_param) {
    eta_vec[((p - 1) * n_stations + 1):(p * n_stations)] = eta[, p];
  }

  // Stage 1 likelihood
  target += normal_prec_chol_lpdf(eta_hat | eta_vec, n_values, index, value, log_det_Q);

  // Non-centered priors on GP latent fields (psi, tau)
  target += std_normal_lpdf(to_vector(eta_raw_gp));

  // Non-centered prior on phi iid noise
  target += std_normal_lpdf(phi_raw);

  // === Priors ===

  // Regression coefficients for psi
  target += normal_lpdf(beta_psi | 0, prior_sd_beta);

  // Intercepts for tau, phi
  target += normal_lpdf(mu_tau | 0, 100);
  target += normal_lpdf(mu_phi | 0, 100);

  // GP marginal SD: Exp(lambda_s) — psi and tau only
  for (p in 1:n_gp)
    target += exponential_lpdf(sigma_gp[p] | lambda_s[p]);

  // GP range: PC prior — psi and tau only
  for (p in 1:n_gp)
    target += log(lambda_rho[p]) - 2 * log(phi_gp[p]) - lambda_rho[p] / phi_gp[p];

  // GP nugget SD: Exp(lambda_nugget) — psi and tau only
  for (p in 1:n_gp)
    target += exponential_lpdf(nugget[p] | lambda_nugget[p]);

  // Phi nugget SD: Exp(lambda_nugget_phi)
  target += exponential_lpdf(nugget_phi | lambda_nugget_phi);
}
