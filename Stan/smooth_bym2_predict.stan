// BYM2 spatial smoothing model with prediction at unobserved locations
//
// Extension of smooth_bym2.stan: adds n_pred prediction nodes to the ICAR
// graph. These nodes have no data (no likelihood contribution) but the ICAR
// prior propagates spatial information from observed stations to grid points.
//
// Indices 1..n_stations = observed stations (contribute to likelihood)
// Indices (n_stations+1)..n_total = prediction grid points (prior only)

functions {
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
    return - 0.5 * dot_self((phi[node1] - phi[node2])) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }

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
  int<lower = 0> n_pred;
  int<lower = 1> n_param;

  // Stage 1 MLEs (stations only)
  vector[n_stations * n_param] eta_hat;

  // Extended ICAR graph (stations + prediction grid)
  int<lower = 1> n_edges;
  array[n_edges] int<lower = 1, upper = n_stations + n_pred> node1;
  array[n_edges] int<lower = 1, upper = n_stations + n_pred> node2;
  real<lower = 0> scaling_factor;

  // Precision matrix from bootstrap covariances (stations only)
  int<lower = 1> n_nonzero_chol_Q;
  array[n_param * n_stations] int n_values;
  array[n_nonzero_chol_Q] int index;
  vector[n_nonzero_chol_Q] value;
  real<lower = 0> log_det_Q;
}

transformed data {
  int n_total = n_stations + n_pred;
}

parameters {
  matrix[n_total, n_param] eta_spatial;
  matrix[n_total, n_param] eta_random;
  vector[n_param] mu;
  vector<lower = 0>[n_param] sigma;
  vector<lower = 0, upper = 1>[n_param] rho;
}

model {
  // Construct eta at OBSERVED stations only (for likelihood)
  vector[n_param * n_stations] eta_obs;

  for (p in 1:n_param) {
    int start = ((p - 1) * n_stations + 1);
    int end = (p * n_stations);
    eta_obs[start:end] = mu[p] + sigma[p] *
      (sqrt(rho[p] / scaling_factor) * eta_spatial[1:n_stations, p] +
       sqrt(1 - rho[p]) * eta_random[1:n_stations, p]);

    // ICAR prior on the FULL graph (stations + grid)
    target += icar_normal_lpdf(eta_spatial[, p] | n_total, node1, node2);
  }

  target += std_normal_lpdf(to_vector(eta_random));
  target += exponential_lpdf(sigma | 1);
  target += beta_lpdf(rho | 1, 1);

  // Likelihood: only at observed stations
  target += normal_prec_chol_lpdf(eta_hat | eta_obs, n_values, index, value, log_det_Q);
}

generated quantities {
  // Full eta at observed stations
  matrix[n_stations, n_param] eta;
  // Predicted eta at grid points
  matrix[n_pred, n_param] eta_pred;

  for (p in 1:n_param) {
    eta[, p] = mu[p] + sigma[p] *
      (sqrt(rho[p] / scaling_factor) * eta_spatial[1:n_stations, p] +
       sqrt(1 - rho[p]) * eta_random[1:n_stations, p]);
    eta_pred[, p] = mu[p] + sigma[p] *
      (sqrt(rho[p] / scaling_factor) * eta_spatial[(n_stations + 1):n_total, p] +
       sqrt(1 - rho[p]) * eta_random[(n_stations + 1):n_total, p]);
  }
}
