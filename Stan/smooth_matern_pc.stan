// Matérn(5/2) GP spatial smoothing with PC priors (Fuglstad et al. 2019)
//
// Same model structure as smooth_matern.stan but with penalized complexity
// priors matching the book chapter (Hazra, Huser & Jóhannesson, 2023):
//   - mu ~ N(0, 100)           (vague intercept)
//   - sigma_gp ~ Exp(lambda_s) (PC prior on marginal SD)
//   - phi_gp ~ PC prior:  pi(rho) = lambda_rho * rho^{-2} * exp(-lambda_rho / rho)
//   - nugget ~ Exp(lambda_nugget)  (PC prior on nugget SD)

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
  int<lower = 1> n_param;

  vector[n_stations * n_param] eta_hat;

  matrix[n_stations, n_stations] dist_mat;

  int<lower = 1> n_nonzero_chol_Q;
  array[n_param * n_stations] int n_values;
  array[n_nonzero_chol_Q] int index;
  vector[n_nonzero_chol_Q] value;
  real<lower = 0> log_det_Q;

  // PC prior rate parameters
  vector<lower = 0>[n_param] lambda_s;       // rate for sigma_gp
  vector<lower = 0>[n_param] lambda_rho;     // rate for phi_gp (range)
  vector<lower = 0>[n_param] lambda_nugget;  // rate for nugget
}

parameters {
  vector[n_param] mu;
  vector<lower = 0>[n_param] sigma_gp;
  vector<lower = 0>[n_param] phi_gp;
  vector<lower = 0>[n_param] nugget;
  matrix[n_stations, n_param] eta_raw;
}

transformed parameters {
  matrix[n_stations, n_param] eta;

  for (p in 1:n_param) {
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

    eta[, p] = mu[p] + cholesky_decompose(Sigma) * eta_raw[, p];
  }
}

model {
  vector[n_param * n_stations] eta_vec;
  for (p in 1:n_param) {
    eta_vec[((p - 1) * n_stations + 1):(p * n_stations)] = eta[, p];
  }

  // Stage 1 likelihood (unchanged)
  target += normal_prec_chol_lpdf(eta_hat | eta_vec, n_values, index, value, log_det_Q);

  // Non-centered prior
  target += std_normal_lpdf(to_vector(eta_raw));

  // === PC Priors (Fuglstad et al. 2019) ===

  // Intercept: vague Gaussian
  target += normal_lpdf(mu | 0, 100);

  // Marginal SD: Exp(lambda_s)
  for (p in 1:n_param)
    target += exponential_lpdf(sigma_gp[p] | lambda_s[p]);

  // Range: PC prior  pi(rho) = lambda * rho^{-2} * exp(-lambda / rho)
  // log pi(rho) = log(lambda) - 2*log(rho) - lambda/rho
  for (p in 1:n_param)
    target += log(lambda_rho[p]) - 2 * log(phi_gp[p]) - lambda_rho[p] / phi_gp[p];

  // Nugget SD: Exp(lambda_nugget)
  for (p in 1:n_param)
    target += exponential_lpdf(nugget[p] | lambda_nugget[p]);
}
