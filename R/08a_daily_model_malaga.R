# 08a_daily_model_malaga.R — Phase 1: Single-station daily hurdle+GPD model
#
# Fits a hurdle model to Málaga Aeropuerto daily rainfall:
#   - Occurrence: logistic with Fourier seasonal terms
#   - Bulk: Gamma (below threshold)
#   - Tail: GPD (above threshold)
#
# Derives implied GEV parameters and compares with existing Max-and-Smooth estimates.
#
# Run from project root: Rscript R/08a_daily_model_malaga.R

library(dplyr)
library(cmdstanr)
library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Phase 1: Daily hurdle+GPD model — M\u00e1laga Aeropuerto\n")
cat("========================================\n")

# ---- 1. Prepare daily data ----
cat("Preparing data...\n")

daily <- readRDS("data/daily_precip_andalucia_raw.rds")
daily$prec_num <- as.numeric(gsub(",", ".", daily$prec))

mal <- daily %>%
  filter(indicativo == "6155A", !is.na(prec_num)) %>%
  mutate(
    doy  = as.integer(format(fecha, "%j")),
    year = as.integer(format(fecha, "%Y"))
  ) %>%
  select(fecha, year, doy, prec_num) %>%
  arrange(fecha)

cat("  Total days:", nrow(mal), "\n")
cat("  Wet days:", sum(mal$prec_num > 0), "\n")

# ---- 2. Fourier seasonal basis ----
K <- 3  # 3 harmonics = 6 coefficients (captures annual + semi-annual + quarter)

make_fourier <- function(doy, K) {
  X <- matrix(0, length(doy), 2 * K)
  for (k in 1:K) {
    X[, 2*k - 1] <- sin(2 * pi * k * doy / 365.25)
    X[, 2*k]     <- cos(2 * pi * k * doy / 365.25)
  }
  X
}

X_all <- make_fourier(mal$doy, K)

# ---- 3. Threshold selection ----
wet <- mal %>% filter(prec_num > 0)
u <- quantile(wet$prec_num, 0.90)  # 90th percentile of wet-day amounts
cat(sprintf("  GPD threshold: %.1f mm (90th pctile of wet days)\n", u))
cat(sprintf("  Exceedances: %d (%.1f/year)\n",
  sum(wet$prec_num > u), sum(wet$prec_num > u) / length(unique(mal$year))))

# ---- 4. Assemble Stan data ----
cat("Assembling Stan data...\n")

is_wet <- as.integer(mal$prec_num > 0)
y_wet  <- wet$prec_num
X_wet  <- make_fourier(wet$doy, K)
above_u <- as.integer(y_wet > u)

stan_data <- list(
  N       = nrow(mal),
  is_wet  = is_wet,
  N_wet   = nrow(wet),
  y_wet   = y_wet,
  u       = as.numeric(u),
  above_u = above_u,
  N_above = sum(above_u),
  K       = K,
  X_season     = X_all,
  X_season_wet = X_wet
)

cat(sprintf("  N=%d, N_wet=%d, N_above=%d, K=%d\n",
  stan_data$N, stan_data$N_wet, stan_data$N_above, stan_data$K))

# ---- 5. Compile and fit ----
cat("Compiling Stan model...\n")
model <- cmdstan_model("Stan/daily_hurdle_gpd.stan")

cat("Fitting model (4 chains, 1000+1000)...\n")
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 42,
  refresh = 200
)

cat("\n---- Diagnostics ----\n")
fit$cmdstan_diagnose()

# ---- 6. Summarise key parameters ----
cat("\n---- Parameter summaries ----\n")

params <- c("alpha_0", "gamma_shape", "log_gamma_rate_0",
            "gpd_sigma", "gpd_xi",
            "gev_mu", "gev_sigma", "gev_xi",
            "rl20", "rl50", "rl100", "lambda_annual")

print(fit$summary(params))

# ---- 7. Compare with existing Max-and-Smooth estimates ----
cat("\n---- Comparison with Max-and-Smooth ----\n")

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")

mal_idx <- which(s1$station_meta$indicativo == "6155A")

# MLE GEV parameters
psi_mle <- s1$mles[mal_idx, "psi"]
tau_mle <- s1$mles[mal_idx, "tau"]
phi_mle <- s1$mles[mal_idx, "phi"]
mu_mle    <- exp(psi_mle)
sigma_mle <- exp(psi_mle + tau_mle)
xi_mle    <- g(phi_mle)

cat(sprintf("  MLE:       mu=%.1f, sigma=%.1f, xi=%.3f\n",
  mu_mle, sigma_mle, xi_mle))

# Smoothed posterior (M&S PC prior)
psi_draws <- s2$psi.selected[, mal_idx]
tau_draws <- s2$tau.selected[, mal_idx]
phi_draws <- s2$phi.selected[, mal_idx]

mu_ms    <- mean(exp(psi_draws))
sigma_ms <- mean(exp(psi_draws + tau_draws))
xi_ms    <- mean(g(phi_draws))

rl100_ms_draws <- exp(psi_draws) +
  exp(psi_draws + tau_draws) / g(phi_draws) *
  ((-log(1 - 1/100))^(-g(phi_draws)) - 1)

cat(sprintf("  M&S smooth: mu=%.1f, sigma=%.1f, xi=%.3f\n",
  mu_ms, sigma_ms, xi_ms))
cat(sprintf("  M&S RL100:  %.1f (%.1f, %.1f)\n",
  mean(rl100_ms_draws),
  quantile(rl100_ms_draws, 0.05),
  quantile(rl100_ms_draws, 0.95)))

# Daily model implied GEV
draws <- fit$draws(format = "df")
cat(sprintf("  Daily GPD:  mu=%.1f (%.1f, %.1f), sigma=%.1f (%.1f, %.1f), xi=%.3f (%.3f, %.3f)\n",
  mean(draws$gev_mu), quantile(draws$gev_mu, 0.05), quantile(draws$gev_mu, 0.95),
  mean(draws$gev_sigma), quantile(draws$gev_sigma, 0.05), quantile(draws$gev_sigma, 0.95),
  mean(draws$gev_xi), quantile(draws$gev_xi, 0.05), quantile(draws$gev_xi, 0.95)))
cat(sprintf("  Daily RL100: %.1f (%.1f, %.1f)\n",
  mean(draws$rl100),
  quantile(draws$rl100, 0.05),
  quantile(draws$rl100, 0.95)))

# ---- 8. Diagnostic plots ----
cat("\nGenerating diagnostic plots...\n")

# Plot 1: Seasonal occurrence probability
doy_grid <- 1:365
X_grid <- make_fourier(doy_grid, K)

alpha_0_draws <- draws$alpha_0
alpha_season_draws <- as.matrix(draws[, paste0("alpha_season[", 1:(2*K), "]")])

logit_p <- sweep(X_grid %*% t(alpha_season_draws), 2, alpha_0_draws, "+")
p_rain <- 1 / (1 + exp(-logit_p))

p_df <- data.frame(
  doy = doy_grid,
  mean = rowMeans(p_rain),
  lo = apply(p_rain, 1, quantile, 0.05),
  hi = apply(p_rain, 1, quantile, 0.95)
)

# Empirical monthly wet-day fraction for comparison
emp_monthly <- mal %>%
  mutate(month = as.integer(format(fecha, "%m"))) %>%
  group_by(month) %>%
  summarise(frac_wet = mean(prec_num > 0), .groups = "drop") %>%
  mutate(doy_mid = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349))

p1 <- ggplot(p_df, aes(x = doy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#3366CC", alpha = 0.2) +
  geom_line(aes(y = mean), colour = "#3366CC", linewidth = 0.8) +
  geom_point(data = emp_monthly, aes(x = doy_mid, y = frac_wet),
             colour = "red", size = 2.5) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  labs(title = "Daily rainfall probability (seasonal cycle)",
       subtitle = "Blue = fitted Fourier model (90% CI) | Red = empirical monthly fraction",
       x = NULL, y = "P(rain > 0)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Plot 2: Return level comparison
rp_grid <- exp(seq(log(2), log(200), length.out = 100))

# Daily model implied RL
rl_daily <- matrix(NA, nrow(draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  rl_daily[, j] <- draws$gev_mu + draws$gev_sigma / draws$gev_xi *
    ((-log(1 - 1/T))^(-draws$gev_xi) - 1)
}

# M&S RL
rl_ms <- matrix(NA, length(psi_draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  mu_d <- exp(psi_draws)
  sig_d <- exp(psi_draws + tau_draws)
  xi_d <- g(phi_draws)
  rl_ms[, j] <- mu_d + sig_d / xi_d * ((-log(1 - 1/T))^(-xi_d) - 1)
}

rl_comp <- data.frame(
  rp = rep(rp_grid, 2),
  mean = c(colMeans(rl_daily), colMeans(rl_ms)),
  lo = c(apply(rl_daily, 2, quantile, 0.05), apply(rl_ms, 2, quantile, 0.05)),
  hi = c(apply(rl_daily, 2, quantile, 0.95), apply(rl_ms, 2, quantile, 0.95)),
  model = rep(c("Daily GPD", "Max-and-Smooth"), each = length(rp_grid))
)

# Empirical return levels from annual maxima
am <- readRDS("data/annual_maxima_andalucia.rds")
mal_am <- am %>% filter(indicativo == "6155A") %>% arrange(max_prec)
n_am <- nrow(mal_am)
mal_am$rp <- rev((n_am + 0.12) / (seq_len(n_am) - 0.44 + 0.12))

p2 <- ggplot(rl_comp, aes(x = rp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15) +
  geom_line(aes(y = mean, colour = model), linewidth = 0.8) +
  geom_point(data = mal_am, aes(x = rp, y = max_prec),
             shape = 21, fill = "grey20", colour = "white", size = 2, stroke = 0.4) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                labels = c("2", "5", "10", "20", "50", "100", "200")) +
  scale_colour_manual(values = c("Daily GPD" = "#E91E63", "Max-and-Smooth" = "#3366CC")) +
  scale_fill_manual(values = c("Daily GPD" = "#E91E63", "Max-and-Smooth" = "#3366CC")) +
  labs(title = "Return level comparison: Daily GPD vs Max-and-Smooth",
       subtitle = "M\u00e1laga Aeropuerto | Points = observed annual maxima | Bands = 90% CI",
       x = "Return period (years)", y = "Daily rainfall (mm)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom")

p_all <- p1 / p2
ggsave("figures/jesus-figures/daily_model_malaga.png",
       p_all, width = 12, height = 10, dpi = 200, bg = "white")

cat("Saved figures/jesus-figures/daily_model_malaga.png\n")

# Save results
saveRDS(list(fit = fit, stan_data = stan_data, threshold = u),
        "data/daily_model_malaga.rds")
cat("Saved data/daily_model_malaga.rds\n")

cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
