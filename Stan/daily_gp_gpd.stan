// daily_gp_gpd.stan — Single-station daily rainfall model with periodic GP
//
// Hurdle model with GP seasonal structure:
//   - Occurrence: logistic with periodic GP over day-of-year
//   - Bulk: Gamma with periodic GP modulating intensity
//   - Tail: GPD above threshold u (constant parameters)
//
// The GP is defined on M knot points (~36, every 10 days) and
// interpolated to all 365 days via precomputed weights.
// This avoids expensive 365x365 Cholesky — only M×M is needed.

functions {
  // GPD log-density
  real gpd_lpdf(real y, real sigma, real xi) {
    if (abs(xi) < 1e-8) {
      return -log(sigma) - y / sigma;
    } else {
      real z = 1.0 + xi * y / sigma;
      if (z <= 0) return negative_infinity();
      return -log(sigma) - (1.0 / xi + 1.0) * log(z);
    }
  }
}

data {
  int<lower=1> N;                // total observed days
  array[N] int<lower=1, upper=365> doy;  // day-of-year for each obs
  array[N] int<lower=0, upper=1> is_wet; // 1 if precip > 0

  int<lower=1> N_wet;            // number of wet days
  vector<lower=0>[N_wet] y_wet;  // precipitation on wet days
  array[N_wet] int<lower=1, upper=365> doy_wet; // day-of-year for wet days
  real<lower=0> u;               // GPD threshold
  array[N_wet] int<lower=0, upper=1> above_u;
  int<lower=0> N_above;          // number of exceedances

  // GP knot configuration (precomputed in R)
  int<lower=2> M;                       // number of GP knots
  matrix[M, M] sin2_mat;                // sin²(π(k_i - k_j)/365.25) between knots
  matrix[365, M] interp_weights;        // linear interpolation from M knots to 365 days
}

parameters {
  // GP hyperparameters — occurrence
  real<lower=0> sigma_occ;
  real<lower=0> ell_occ;

  // GP hyperparameters — intensity
  real<lower=0> sigma_int;
  real<lower=0> ell_int;

  // GP latent values at knots (non-centered)
  vector[M] eta_occ;
  vector[M] eta_int;

  // Gamma bulk shape
  real<lower=0> gamma_shape;

  // GPD tail
  real<lower=0> gpd_sigma;
  real<lower=-0.5, upper=1.0> gpd_xi;
}

transformed parameters {
  // GP seasonal curves for all 365 days
  vector[365] f_occ;
  vector[365] f_int;

  {
    // Build M×M periodic kernel matrices
    real sigma2_occ = square(sigma_occ);
    real sigma2_int = square(sigma_int);
    real inv_ell2_occ = 1.0 / square(ell_occ);
    real inv_ell2_int = 1.0 / square(ell_int);

    matrix[M, M] K_occ;
    matrix[M, M] K_int;

    for (i in 1:M) {
      K_occ[i, i] = sigma2_occ + 1e-6;
      K_int[i, i] = sigma2_int + 1e-6;
      for (j in (i+1):M) {
        K_occ[i, j] = sigma2_occ * exp(-2.0 * sin2_mat[i, j] * inv_ell2_occ);
        K_occ[j, i] = K_occ[i, j];
        K_int[i, j] = sigma2_int * exp(-2.0 * sin2_mat[i, j] * inv_ell2_int);
        K_int[j, i] = K_int[i, j];
      }
    }

    // GP at knots (non-centered)
    vector[M] f_knots_occ = cholesky_decompose(K_occ) * eta_occ;
    vector[M] f_knots_int = cholesky_decompose(K_int) * eta_int;

    // Interpolate to all 365 days
    f_occ = interp_weights * f_knots_occ;
    f_int = interp_weights * f_knots_int;
  }
}

model {
  // --- Priors ---

  // GP hyperparameters
  sigma_occ ~ normal(0, 2);
  ell_occ ~ normal(0, 3);
  sigma_int ~ normal(0, 2);
  ell_int ~ normal(0, 3);

  // Non-centered GP
  eta_occ ~ std_normal();
  eta_int ~ std_normal();

  // Bulk
  gamma_shape ~ exponential(0.2);

  // GPD tail
  gpd_sigma ~ exponential(0.05);
  gpd_xi ~ normal(0.1, 0.3);

  // --- Likelihood ---

  // Occurrence: logistic(f_occ[doy])
  for (i in 1:N) {
    is_wet[i] ~ bernoulli_logit(f_occ[doy[i]]);
  }

  // Wet-day amounts
  for (i in 1:N_wet) {
    real rate_i = exp(-f_int[doy_wet[i]]);

    if (above_u[i] == 0) {
      target += gamma_lpdf(y_wet[i] | gamma_shape, rate_i);
      target += -gamma_lcdf(u | gamma_shape, rate_i);
    } else {
      real excess = y_wet[i] - u;
      target += gpd_lpdf(excess | gpd_sigma, gpd_xi);
    }
  }
}

generated quantities {
  // Implied GEV via Pickands-Balkema-de Haan
  real lambda_annual = N_above * 365.25 / (N * 1.0);
  real gev_xi = gpd_xi;
  real gev_sigma;
  real gev_mu;

  if (abs(gpd_xi) < 1e-8) {
    gev_sigma = gpd_sigma;
    gev_mu = u + gpd_sigma * log(lambda_annual);
  } else {
    gev_sigma = gpd_sigma * (lambda_annual ^ gpd_xi);
    gev_mu = u + gpd_sigma / gpd_xi * ((lambda_annual ^ gpd_xi) - 1);
  }

  real rl20  = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/20))^(-gev_xi) - 1);
  real rl50  = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/50))^(-gev_xi) - 1);
  real rl100 = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/100))^(-gev_xi) - 1);
}
