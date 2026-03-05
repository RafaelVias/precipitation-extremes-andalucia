// BYM2 spatial smoothing model for Max-and-Smooth
// Adapted from Brynjolfur Gauti Jonsson's maxandsmooth R package
// https://github.com/bgautijonsson/maxandsmooth
//
// Replaces the custom MCMC SPDE sampler with Stan/NUTS + BYM2 (ICAR + iid)

functions {
  /*
  ICAR prior on spatial coefficients with soft sum-to-zero constraint
  */
  real icar_normal_lpdf(vector phi, int N, array[] int node1, array[] int node2) {
    return - 0.5 * dot_self((phi[node1] - phi[node2])) +
      normal_lpdf(sum(phi) | 0, 0.001 * N);
  }

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
  int<lower = 1> n_obs;
  int<lower = 1> n_param;

  vector[n_stations * n_param] eta_hat;

  int<lower = 1> n_edges;
  array[n_edges] int<lower = 1, upper = n_stations> node1;
  array[n_edges] int<lower = 1, upper = n_stations> node2;
  real<lower = 0> scaling_factor;

  int<lower = 1> n_nonzero_chol_Q;
  array[n_param * n_stations] int n_values;
  array[n_nonzero_chol_Q] int index;
  vector[n_nonzero_chol_Q] value;
  real<lower = 0> log_det_Q;
}

transformed data {
  vector[n_stations] psi_hat = eta_hat[1:n_stations];
  vector[n_stations] tau_hat = eta_hat[(n_stations + 1):(2 * n_stations)];
  vector[n_stations] phi_hat = eta_hat[(2 * n_stations + 1):(3 * n_stations)];
}

parameters {
  matrix[n_stations, n_param] eta_spatial;
  matrix[n_stations, n_param] eta_random;
  vector[n_param] mu;
  vector<lower = 0>[n_param] sigma;
  vector<lower = 0, upper = 1>[n_param] rho;
}

model {
  vector[n_param * n_stations] eta;

  for (p in 1:n_param) {
    int start = ((p - 1) * n_stations + 1);
    int end = (p * n_stations);
    eta[start:end] = mu[p] + sigma[p] *
      (sqrt(rho[p] / scaling_factor) * eta_spatial[, p] +
       sqrt(1 - rho[p]) * eta_random[, p]);
    target += icar_normal_lpdf(eta_spatial[, p] | n_edges, node1, node2);
  }

  target += std_normal_lpdf(to_vector(eta_random));
  target += exponential_lpdf(sigma | 1);
  target += beta_lpdf(rho | 1, 1);

  target += normal_prec_chol_lpdf(eta_hat | eta, n_values, index, value, log_det_Q);
}

generated quantities {
  matrix[n_stations, n_param] eta;
  for (p in 1:n_param) {
    eta[, p] = mu[p] + sigma[p] * (sqrt(rho[p] / scaling_factor) * eta_spatial[, p] + sqrt(1 - rho[p]) * eta_random[, p]);
  }
}
