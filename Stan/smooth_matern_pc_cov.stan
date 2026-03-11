// Matern(5/2) GP spatial smoothing with PC priors and covariates on psi
//
// Model structure (following Hazra, Huser & Johannesson, BLGM 2023 Ch.7):
//   eta_psi(s) = X(s)' beta_psi + f_psi(s)   (covariates + GP)
//   eta_tau(s) = mu_tau + f_tau(s)             (intercept + GP)
//   eta_phi(s) = mu_phi + f_phi(s)             (intercept + GP)
//
// All three parameters receive Matern(5/2) GPs. The phi GP uses tighter
// PC priors reflecting the smaller spatial variation in the shape parameter.
//
// Covariates on psi: intercept, altitude (std), exposure (std), alt x exposure
//
// Priors:
//   - beta_psi ~ N(0, prior_sd_beta)
//   - mu_tau, mu_phi ~ N(0, 100)
//   - sigma_gp[1:3] ~ Exp(lambda_s)     (separate rates per parameter)
//   - phi_gp[1:3] ~ PC prior            (separate rates per parameter)
//   - nugget[1:3] ~ Exp(lambda_nugget)   (separate rates per parameter)

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
  int<lower = 1> n_param;           // 3 (psi, tau, phi)
  int<lower = 1> n_gp;              // 3 (all get GPs)

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

  // PC prior rate parameters (length 3: psi, tau, phi)
  vector<lower = 0>[n_gp] lambda_s;
  vector<lower = 0>[n_gp] lambda_rho;
  vector<lower = 0>[n_gp] lambda_nugget;
}

parameters {
  vector[n_cov_psi] beta_psi;
  real mu_tau;
  real mu_phi;
  vector<lower = 0>[n_gp] sigma_gp;    // GP marginal SD
  vector<lower = 0>[n_gp] phi_gp;      // GP range
  vector<lower = 0>[n_gp] nugget;      // GP nugget
  matrix[n_stations, n_gp] eta_raw;    // non-centered GP draws
}

transformed parameters {
  matrix[n_stations, n_param] eta;

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
      eta[, p] = X_psi * beta_psi + cholesky_decompose(Sigma) * eta_raw[, p];
    else if (p == 2)
      eta[, p] = mu_tau + cholesky_decompose(Sigma) * eta_raw[, p];
    else
      eta[, p] = mu_phi + cholesky_decompose(Sigma) * eta_raw[, p];
  }
}

model {
  vector[n_param * n_stations] eta_vec;
  for (p in 1:n_param) {
    eta_vec[((p - 1) * n_stations + 1):(p * n_stations)] = eta[, p];
  }

  // Stage 1 likelihood
  target += normal_prec_chol_lpdf(eta_hat | eta_vec, n_values, index, value, log_det_Q);

  // Non-centered priors on GP latent fields
  target += std_normal_lpdf(to_vector(eta_raw));

  // === Priors ===

  // Regression coefficients for psi
  target += normal_lpdf(beta_psi | 0, prior_sd_beta);

  // Intercepts for tau, phi
  target += normal_lpdf(mu_tau | 0, 100);
  target += normal_lpdf(mu_phi | 0, 100);

  // GP marginal SD: Exp(lambda_s)
  for (p in 1:n_gp)
    target += exponential_lpdf(sigma_gp[p] | lambda_s[p]);

  // GP range: PC prior
  for (p in 1:n_gp)
    target += log(lambda_rho[p]) - 2 * log(phi_gp[p]) - lambda_rho[p] / phi_gp[p];

  // GP nugget SD: Exp(lambda_nugget)
  for (p in 1:n_gp)
    target += exponential_lpdf(nugget[p] | lambda_nugget[p]);
}
