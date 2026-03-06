// daily_hurdle_gpd.stan — Single-station daily rainfall model
//
// Hurdle model:
//   Component 1: P(rain > 0) — logistic with seasonal Fourier terms
//   Component 2: Y | Y > 0 — bulk (Gamma) + tail (GPD above threshold u)
//
// Seasonal structure via Fourier harmonics of day-of-year (K harmonics).
// The GPD shape xi is constant (no seasonal variation for the tail shape).

functions {
  // Log-PDF of Generalized Pareto Distribution
  // y > 0, sigma > 0, xi real
  real gpd_lpdf(real y, real sigma, real xi) {
    if (abs(xi) < 1e-8) {
      // Exponential limit
      return -log(sigma) - y / sigma;
    } else {
      real z = 1.0 + xi * y / sigma;
      if (z <= 0) return negative_infinity();
      return -log(sigma) - (1.0 / xi + 1.0) * log(z);
    }
  }

  // Log-CDF of GPD (for computing survival function)
  real gpd_lcdf(real y, real sigma, real xi) {
    if (abs(xi) < 1e-8) {
      return log1m_exp(-y / sigma);
    } else {
      real z = 1.0 + xi * y / sigma;
      if (z <= 0) return 0;  // CDF = 1
      return log1m(z ^ (-1.0 / xi));
    }
  }
}

data {
  int<lower=1> N;               // total observed days
  array[N] int<lower=0, upper=1> is_wet; // 1 if precip > 0
  int<lower=1> N_wet;           // number of wet days
  vector<lower=0>[N_wet] y_wet; // precipitation on wet days (mm)
  real<lower=0> u;              // GPD threshold (mm, on wet-day amounts)
  array[N_wet] int<lower=0, upper=1> above_u; // 1 if y_wet > u
  int<lower=0> N_above;         // number of exceedances

  // Seasonal covariates: sin/cos pairs for each day
  int<lower=1> K;               // number of Fourier harmonics
  matrix[N, 2*K] X_season;      // Fourier basis for all days
  matrix[N_wet, 2*K] X_season_wet; // Fourier basis for wet days only
}

parameters {
  // Occurrence model (logistic)
  real alpha_0;                  // intercept
  vector[2*K] alpha_season;     // Fourier coefficients

  // Bulk model (Gamma for wet days below threshold)
  real<lower=0> gamma_shape;    // Gamma shape
  real log_gamma_rate_0;        // log(rate) intercept
  vector[2*K] gamma_rate_season; // seasonal modulation of rate

  // GPD tail model (exceedances above u)
  real<lower=0> gpd_sigma;      // GPD scale
  real<lower=-0.5, upper=1.0> gpd_xi; // GPD shape
}

model {
  // --- Priors ---

  // Occurrence
  alpha_0 ~ normal(-2, 2);       // ~84% dry days -> logit(0.16) ≈ -1.7
  alpha_season ~ normal(0, 2);

  // Gamma bulk
  gamma_shape ~ exponential(0.2);   // weakly informative
  log_gamma_rate_0 ~ normal(0, 2);
  gamma_rate_season ~ normal(0, 1);

  // GPD tail
  gpd_sigma ~ exponential(0.05);  // weakly informative, allows up to ~50mm scale
  gpd_xi ~ normal(0.1, 0.3);     // centered on light positive tail

  // --- Likelihood ---

  // Component 1: Occurrence (logistic regression)
  {
    vector[N] logit_p = alpha_0 + X_season * alpha_season;
    is_wet ~ bernoulli_logit(logit_p);
  }

  // Component 2: Wet-day amounts
  {
    vector[N_wet] log_rate = log_gamma_rate_0 + X_season_wet * gamma_rate_season;

    for (i in 1:N_wet) {
      if (above_u[i] == 0) {
        // Below threshold: Gamma likelihood, truncated to (0, u]
        // y_wet[i] ~ Gamma(shape, rate) truncated to (0, u]
        real rate_i = exp(log_rate[i]);
        target += gamma_lpdf(y_wet[i] | gamma_shape, rate_i);
        target += -gamma_lcdf(u | gamma_shape, rate_i);  // truncation adjustment
      } else {
        // Above threshold: GPD likelihood for exceedance (y - u)
        real excess = y_wet[i] - u;
        target += gpd_lpdf(excess | gpd_sigma, gpd_xi);
      }
    }
  }
}

generated quantities {
  // Implied GEV parameters via Pickands-Balkema-de Haan
  // lambda = expected number of threshold exceedances per year (~365 days)
  // Approximate: lambda ≈ N_above / (N / 365.25)
  real lambda_annual = N_above * 365.25 / (N * 1.0);

  // GEV parameters derived from GPD + exceedance rate
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

  // Return levels
  real rl20  = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/20))^(-gev_xi) - 1);
  real rl50  = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/50))^(-gev_xi) - 1);
  real rl100 = gev_mu + gev_sigma / gev_xi * ((-log(1 - 1.0/100))^(-gev_xi) - 1);
}
