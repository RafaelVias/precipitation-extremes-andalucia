// Matérn(5/2) GP spatial smoothing model for Max-and-Smooth
//
// Replaces BYM2 (ICAR + iid) with a proper Matérn(5/2) Gaussian Process.
// Each GEV parameter field gets independent GP hyperparameters (sigma, phi, nugget).
// Non-centered parameterization: eta = mu + L * eta_raw, eta_raw ~ N(0, I).
//
// The Stage 1 likelihood (normal_prec_chol_lpdf) is unchanged from smooth_bym2.stan.

functions {
  /*
  Log-density of N(y | x, Q^{-1}) using sparse Cholesky factor L where Q = LL^T.
  Stores L in compressed sparse column (CSC) format:
    n_values[i] = number of non-zeros in column i
    index[k]    = row index of k-th non-zero
    values[k]   = value of k-th non-zero
  Computes: log|L| - 0.5 * ||L^T(y - x)||^2
  */
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

  // Pairwise distance matrix (precomputed in R)
  matrix[n_stations, n_stations] dist_mat;

  // Stage 1 precision matrix in sparse Cholesky CSC format
  int<lower = 1> n_nonzero_chol_Q;
  array[n_param * n_stations] int n_values;
  array[n_nonzero_chol_Q] int index;
  vector[n_nonzero_chol_Q] value;
  real<lower = 0> log_det_Q;
}

parameters {
  vector[n_param] mu;
  vector<lower = 0>[n_param] sigma_gp;   // marginal SD (partial sill)
  vector<lower = 0>[n_param] phi_gp;     // range parameter
  vector<lower = 0>[n_param] nugget;     // nugget SD
  matrix[n_stations, n_param] eta_raw;   // non-centered latent field
}

transformed parameters {
  matrix[n_stations, n_param] eta;

  for (p in 1:n_param) {
    // Build Matérn(5/2) covariance matrix
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

    // Non-centered: eta = mu + L * eta_raw
    eta[, p] = mu[p] + cholesky_decompose(Sigma) * eta_raw[, p];
  }
}

model {
  // Flatten eta to parameter-major vector for likelihood
  vector[n_param * n_stations] eta_vec;
  for (p in 1:n_param) {
    eta_vec[((p - 1) * n_stations + 1):(p * n_stations)] = eta[, p];
  }

  // Stage 1 likelihood (unchanged)
  target += normal_prec_chol_lpdf(eta_hat | eta_vec, n_values, index, value, log_det_Q);

  // Non-centered prior
  target += std_normal_lpdf(to_vector(eta_raw));

  // Hyperpriors
  target += normal_lpdf(mu | 0, 10);
  target += exponential_lpdf(sigma_gp | 1);
  target += inv_gamma_lpdf(phi_gp | 3.0, 1.0);
  target += exponential_lpdf(nugget | 2);
}
