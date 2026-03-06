# 08b_daily_gp_malaga.R — Phase 1: Single-station daily model with periodic GP
#
# Fits a hurdle model to Málaga Aeropuerto daily rainfall:
#   - Occurrence: logistic with periodic GP over day-of-year
#   - Bulk: Gamma with periodic GP modulating intensity
#   - Tail: GPD (above threshold)
#
# The GP is defined on M=36 knot points (~every 10 days) and
# interpolated to all 365 days. This makes the Cholesky 36×36
# instead of 365×365 — roughly 1000× faster.
#
# Run from project root: Rscript R/08b_daily_gp_malaga.R

library(dplyr)
library(cmdstanr)
library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("Phase 1: Daily GP model \u2014 M\u00e1laga Aeropuerto\n")
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

# Cap doy at 365 (leap year day 366 -> 365)
mal$doy <- pmin(mal$doy, 365L)

cat("  Total days:", nrow(mal), "\n")
cat("  Wet days:", sum(mal$prec_num > 0), "\n")

# ---- 2. Threshold selection ----
wet <- mal %>% filter(prec_num > 0)
u <- as.numeric(quantile(wet$prec_num, 0.90))
cat(sprintf("  GPD threshold: %.1f mm\n", u))
cat(sprintf("  Exceedances: %d (%.1f/year)\n",
  sum(wet$prec_num > u), sum(wet$prec_num > u) / length(unique(mal$year))))

# ---- 3. GP knot configuration ----
cat("Setting up GP knots...\n")

M <- 36L
knot_days <- seq(5, 360, length.out = M)
cat(sprintf("  M=%d knots, spacing ~%.1f days\n", M, diff(knot_days)[1]))

# Precompute sin² matrix between knots (for periodic kernel)
sin2_mat <- matrix(0, M, M)
for (i in 1:M) {
  for (j in 1:M) {
    sin_val <- sin(pi * (knot_days[i] - knot_days[j]) / 365.25)
    sin2_mat[i, j] <- sin_val^2
  }
}

# Build 365 × M linear interpolation matrix (circular)
interp_weights <- matrix(0, 365, M)
for (d in 1:365) {
  # Circular distances to all knots
  cdist <- abs(d - knot_days)
  cdist <- pmin(cdist, 365 - cdist)

  # Find two nearest knots
  ord <- order(cdist)
  i1 <- ord[1]; i2 <- ord[2]
  d1 <- cdist[i1]; d2 <- cdist[i2]

  if (d1 < 1e-10) {
    interp_weights[d, i1] <- 1.0
  } else {
    interp_weights[d, i1] <- d2 / (d1 + d2)
    interp_weights[d, i2] <- d1 / (d1 + d2)
  }
}

# Verify: each row sums to 1
stopifnot(all(abs(rowSums(interp_weights) - 1) < 1e-10))

# ---- 4. Assemble Stan data ----
cat("Assembling Stan data...\n")

is_wet <- as.integer(mal$prec_num > 0)
y_wet  <- wet$prec_num
doy_wet <- pmin(as.integer(format(wet$fecha, "%j")), 365L)
above_u <- as.integer(y_wet > u)

stan_data <- list(
  N       = nrow(mal),
  doy     = mal$doy,
  is_wet  = is_wet,
  N_wet   = nrow(wet),
  y_wet   = y_wet,
  doy_wet = doy_wet,
  u       = u,
  above_u = above_u,
  N_above = sum(above_u),
  M       = M,
  sin2_mat = sin2_mat,
  interp_weights = interp_weights
)

cat(sprintf("  N=%d, N_wet=%d, N_above=%d\n",
  stan_data$N, stan_data$N_wet, stan_data$N_above))

# ---- 5. Compile and fit ----
cat("Compiling Stan model...\n")
model <- cmdstan_model("Stan/daily_gp_gpd.stan")

# Initial values (now M-dimensional GP vectors)
inits <- list(
  sigma_occ = 1.0, ell_occ = 1.5,
  sigma_int = 1.0, ell_int = 1.5,
  eta_occ = rep(0, M), eta_int = rep(0, M),
  gamma_shape = 0.7,
  gpd_sigma = 18.0, gpd_xi = 0.1
)

cat("Fitting model (4 chains, 500+500)...\n")
fit <- model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 42,
  refresh = 100,
  init = rep(list(inits), 4)
)

cat("\n---- Diagnostics ----\n")
fit$cmdstan_diagnose()

# ---- 6. Summarise key parameters ----
cat("\n---- Parameter summaries ----\n")

params <- c("sigma_occ", "ell_occ", "sigma_int", "ell_int",
            "gamma_shape", "gpd_sigma", "gpd_xi",
            "gev_mu", "gev_sigma", "gev_xi",
            "rl20", "rl50", "rl100")

print(fit$summary(params))

# ---- 7. Compare with Max-and-Smooth ----
cat("\n---- Comparison with Max-and-Smooth ----\n")

s1 <- readRDS("data/stage1_results.rds")
s2 <- readRDS("data/stage2_matern_pc_results.rds")
mal_idx <- which(s1$station_meta$indicativo == "6155A")

psi_draws <- s2$psi.selected[, mal_idx]
tau_draws <- s2$tau.selected[, mal_idx]
phi_draws <- s2$phi.selected[, mal_idx]

rl100_ms <- exp(psi_draws) +
  exp(psi_draws + tau_draws) / g(phi_draws) *
  ((-log(1 - 1/100))^(-g(phi_draws)) - 1)

draws <- fit$draws(format = "df")

cat(sprintf("  M&S:       RL100 = %.1f (%.1f, %.1f)\n",
  mean(rl100_ms), quantile(rl100_ms, 0.05), quantile(rl100_ms, 0.95)))
cat(sprintf("  Daily GP:  RL100 = %.1f (%.1f, %.1f)\n",
  mean(draws$rl100), quantile(draws$rl100, 0.05), quantile(draws$rl100, 0.95)))
cat(sprintf("  Daily GP:  xi = %.3f (%.3f, %.3f)\n",
  mean(draws$gpd_xi), quantile(draws$gpd_xi, 0.05), quantile(draws$gpd_xi, 0.95)))

# ---- 8. Diagnostic plots ----
cat("\nGenerating plots...\n")

# Extract seasonal curves
f_occ_draws <- fit$draws("f_occ", format = "matrix")  # 2000 x 365
f_int_draws <- fit$draws("f_int", format = "matrix")

# Occurrence probability curve
p_rain <- 1 / (1 + exp(-f_occ_draws))

occ_df <- data.frame(
  doy = 1:365,
  mean = colMeans(p_rain),
  lo = apply(p_rain, 2, quantile, 0.05),
  hi = apply(p_rain, 2, quantile, 0.95)
)

# Empirical for comparison
emp_monthly <- mal %>%
  mutate(month = as.integer(format(fecha, "%m"))) %>%
  group_by(month) %>%
  summarise(frac_wet = mean(prec_num > 0), .groups = "drop") %>%
  mutate(doy_mid = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349))

p1 <- ggplot(occ_df, aes(x = doy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#3366CC", alpha = 0.2) +
  geom_line(aes(y = mean), colour = "#3366CC", linewidth = 0.8) +
  geom_point(data = emp_monthly, aes(x = doy_mid, y = frac_wet),
             colour = "red", size = 2.5) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  labs(title = "Daily rainfall probability \u2014 periodic GP seasonal curve",
       subtitle = "Blue = GP posterior (90% CI) | Red = empirical monthly fraction",
       x = NULL, y = "P(rain > 0)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Intensity curve (mean rainfall given rain)
gamma_shape_draws <- draws$gamma_shape
gamma_mean <- sweep(exp(f_int_draws), 1, gamma_shape_draws, "*")

int_df <- data.frame(
  doy = 1:365,
  mean = colMeans(gamma_mean),
  lo = apply(gamma_mean, 2, quantile, 0.05),
  hi = apply(gamma_mean, 2, quantile, 0.95)
)

# Empirical mean wet-day amount by month
emp_int <- mal %>%
  filter(prec_num > 0) %>%
  mutate(month = as.integer(format(fecha, "%m"))) %>%
  group_by(month) %>%
  summarise(mean_wet = mean(prec_num), .groups = "drop") %>%
  mutate(doy_mid = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349))

p2 <- ggplot(int_df, aes(x = doy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#FF9800", alpha = 0.2) +
  geom_line(aes(y = mean), colour = "#FF9800", linewidth = 0.8) +
  geom_point(data = emp_int, aes(x = doy_mid, y = mean_wet),
             colour = "red", size = 2.5) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  labs(title = "Mean wet-day rainfall \u2014 periodic GP seasonal curve",
       subtitle = "Orange = GP posterior (90% CI) | Red = empirical monthly mean",
       x = NULL, y = "Mean rainfall (mm) given rain") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

# Return level comparison
rp_grid <- exp(seq(log(2), log(200), length.out = 100))

rl_daily <- matrix(NA, nrow(draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  rl_daily[, j] <- draws$gev_mu + draws$gev_sigma / draws$gev_xi *
    ((-log(1 - 1/T))^(-draws$gev_xi) - 1)
}

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
  model = rep(c("Daily GP", "Max-and-Smooth"), each = length(rp_grid))
)

am <- readRDS("data/annual_maxima_andalucia.rds")
mal_am <- am %>% filter(indicativo == "6155A") %>% arrange(max_prec)
n_am <- nrow(mal_am)
mal_am$rp <- rev((n_am + 0.12) / (seq_len(n_am) - 0.44 + 0.12))

p3 <- ggplot(rl_comp, aes(x = rp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15) +
  geom_line(aes(y = mean, colour = model), linewidth = 0.8) +
  geom_point(data = mal_am, aes(x = rp, y = max_prec),
             shape = 21, fill = "grey20", colour = "white", size = 2, stroke = 0.4) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                labels = c("2", "5", "10", "20", "50", "100", "200")) +
  scale_colour_manual(values = c("Daily GP" = "#E91E63", "Max-and-Smooth" = "#3366CC")) +
  scale_fill_manual(values = c("Daily GP" = "#E91E63", "Max-and-Smooth" = "#3366CC")) +
  labs(title = "Return level comparison: Daily GP vs Max-and-Smooth",
       subtitle = "M\u00e1laga Aeropuerto | Points = observed annual maxima",
       x = "Return period (years)", y = "Daily rainfall (mm)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

p_all <- p1 / p2 / p3
ggsave("figures/jesus-figures/daily_gp_malaga.png",
       p_all, width = 12, height = 13, dpi = 200, bg = "white")

cat("Saved figures/jesus-figures/daily_gp_malaga.png\n")

saveRDS(list(fit = fit, stan_data = stan_data, threshold = u,
             knot_days = knot_days, interp_weights = interp_weights),
        "data/daily_gp_malaga.rds")
cat("Saved data/daily_gp_malaga.rds\n")

cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
