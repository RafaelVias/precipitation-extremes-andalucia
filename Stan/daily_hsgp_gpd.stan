// daily_hsgp_gpd.stan — HSGP approximation for periodic GP
//
// Same model as daily_gp_gpd.stan but replaces the knot-based GP
// with a Hilbert Space GP (Fourier basis + spectral density weights).
//
// For the periodic kernel k(τ) = σ²exp(-2sin²(πτ/P)/ℓ²):
//   eigenvalues:  λ_m = σ² exp(-κ) I_m(κ),  κ = 2/ℓ²
//   eigenfunctions: Fourier basis {1, cos(2πjd/P), sin(2πjd/P), ...}
//
// No Cholesky needed — f = PHI * (diagSPD .* z)

functions {
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
  int<lower=1> N;
  array[N] int<lower=1, upper=365> doy;
  array[N] int<lower=0, upper=1> is_wet;

  int<lower=1> N_wet;
  vector<lower=0>[N_wet] y_wet;
  array[N_wet] int<lower=1, upper=365> doy_wet;
  real<lower=0> u;
  array[N_wet] int<lower=0, upper=1> above_u;
  int<lower=0> N_above;

  int<lower=1> J;  // number of Fourier harmonics
}

transformed data {
  int B = 1 + 2 * J;   // total basis functions: 1 constant + J cos/sin pairs
  real P = 365.25;

  // Precompute Fourier basis matrix: 365 × B
  matrix[365, B] PHI;
  for (d in 1:365) {
    PHI[d, 1] = 1.0;
    for (j in 1:J) {
      PHI[d, 2 * j]     = cos(2.0 * pi() * j * d / P);
      PHI[d, 2 * j + 1] = sin(2.0 * pi() * j * d / P);
    }
  }
}

parameters {
  real<lower=0> sigma_occ;
  real<lower=0> ell_occ;
  real<lower=0> sigma_int;
  real<lower=0> ell_int;

  vector[B] z_occ;
  vector[B] z_int;

  real<lower=0> gamma_shape;
  real<lower=0> gpd_sigma;
  real<lower=-0.5, upper=1.0> gpd_xi;
}

transformed parameters {
  vector[365] f_occ;
  vector[365] f_int;

  {
    // Spectral density for periodic kernel:
    //   λ_m = σ² exp(-κ) I_m(κ),  κ = 2/ℓ²
    // diagSPD = sqrt(λ_m) for each basis function

    real kappa_occ = 2.0 / square(ell_occ);
    real kappa_int = 2.0 / square(ell_int);

    vector[B] spd_occ;
    vector[B] spd_int;

    // j=0: constant basis
    spd_occ[1] = sigma_occ * sqrt(exp(-kappa_occ) *
                  modified_bessel_first_kind(0, kappa_occ));
    spd_int[1] = sigma_int * sqrt(exp(-kappa_int) *
                  modified_bessel_first_kind(0, kappa_int));

    // j=1..J: cos/sin pairs share the same spectral weight
    for (j in 1:J) {
      real w_occ = sigma_occ * sqrt(2.0 * exp(-kappa_occ) *
                    modified_bessel_first_kind(j, kappa_occ));
      real w_int = sigma_int * sqrt(2.0 * exp(-kappa_int) *
                    modified_bessel_first_kind(j, kappa_int));
      spd_occ[2 * j]     = w_occ;
      spd_occ[2 * j + 1] = w_occ;
      spd_int[2 * j]     = w_int;
      spd_int[2 * j + 1] = w_int;
    }

    // GP approximation: f = PHI * (diagSPD .* z)
    f_occ = PHI * (spd_occ .* z_occ);
    f_int = PHI * (spd_int .* z_int);
  }
}

model {
  // --- Priors (same as knot-based version) ---
  sigma_occ ~ normal(0, 2);
  ell_occ ~ normal(0, 3);
  sigma_int ~ normal(0, 2);
  ell_int ~ normal(0, 3);

  z_occ ~ std_normal();
  z_int ~ std_normal();

  gamma_shape ~ exponential(0.2);
  gpd_sigma ~ exponential(0.05);
  gpd_xi ~ normal(0.1, 0.3);

  // --- Likelihood ---
  for (i in 1:N) {
    is_wet[i] ~ bernoulli_logit(f_occ[doy[i]]);
  }

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
