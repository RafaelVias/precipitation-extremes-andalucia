# 08c_daily_hsgp_malaga.R — HSGP version of the daily GP model
#
# Same model as 08b but using Hilbert Space GP approximation:
#   - Fourier basis functions precomputed in transformed data
#   - Spectral density weights (modified Bessel) instead of Cholesky
#   - No knot matrices needed
#
# Run from project root: Rscript R/08c_daily_hsgp_malaga.R

library(dplyr)
library(cmdstanr)
library(ggplot2)
library(patchwork)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

source("vendor/max_and_smooth/stage1_functions.R")

cat("========================================\n")
cat("HSGP daily model \u2014 M\u00e1laga Aeropuerto\n")
cat("========================================\n")

# ---- 1. Prepare daily data (identical to 08b) ----
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

mal$doy <- pmin(mal$doy, 365L)

cat("  Total days:", nrow(mal), "\n")
cat("  Wet days:", sum(mal$prec_num > 0), "\n")

# ---- 2. Threshold selection ----
wet <- mal %>% filter(prec_num > 0)
u <- as.numeric(quantile(wet$prec_num, 0.90))
cat(sprintf("  GPD threshold: %.1f mm\n", u))
cat(sprintf("  Exceedances: %d (%.1f/year)\n",
  sum(wet$prec_num > u), sum(wet$prec_num > u) / length(unique(mal$year))))

# ---- 3. HSGP: just set J ----
J <- 20L  # 20 Fourier harmonics → 41 basis functions
B <- 1L + 2L * J
cat(sprintf("  J=%d harmonics \u2192 %d basis functions (no knots, no Cholesky)\n", J, B))

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
  J       = J
)

cat(sprintf("  N=%d, N_wet=%d, N_above=%d\n",
  stan_data$N, stan_data$N_wet, stan_data$N_above))

# ---- 5. Compile and fit ----
cat("Compiling Stan model...\n")
model <- cmdstan_model("Stan/daily_hsgp_gpd.stan")

inits <- list(
  sigma_occ = 1.0, ell_occ = 1.5,
  sigma_int = 1.0, ell_int = 1.5,
  z_occ = rep(0, B), z_int = rep(0, B),
  gamma_shape = 0.7,
  gpd_sigma = 18.0, gpd_xi = 0.1
)

cat("Fitting HSGP model (4 chains, 500+500)...\n")
t0 <- proc.time()

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

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n  Wall time: %.1f seconds (%.1f minutes)\n", elapsed, elapsed / 60))

cat("\n---- Diagnostics ----\n")
fit$cmdstan_diagnose()

# ---- 6. Summarise key parameters ----
cat("\n---- Parameter summaries ----\n")

params <- c("sigma_occ", "ell_occ", "sigma_int", "ell_int",
            "gamma_shape", "gpd_sigma", "gpd_xi",
            "gev_mu", "gev_sigma", "gev_xi",
            "rl20", "rl50", "rl100")

print(fit$summary(params))

# ---- 7. Compare with knot-based GP ----
cat("\n---- Comparison: HSGP vs knot-based GP vs Max-and-Smooth ----\n")

# Load knot-based results
knot_res <- readRDS("data/daily_gp_malaga.rds")
knot_draws <- knot_res$fit$draws(format = "df")

# Load M&S results
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

fmt <- "  %-12s RL100 = %6.1f (%5.1f, %5.1f) | xi = %.3f (%.3f, %.3f)\n"
cat(sprintf(fmt, "M&S:",
  mean(rl100_ms), quantile(rl100_ms, 0.05), quantile(rl100_ms, 0.95),
  mean(g(phi_draws)), quantile(g(phi_draws), 0.05), quantile(g(phi_draws), 0.95)))
cat(sprintf(fmt, "Knot GP:",
  mean(knot_draws$rl100), quantile(knot_draws$rl100, 0.05), quantile(knot_draws$rl100, 0.95),
  mean(knot_draws$gpd_xi), quantile(knot_draws$gpd_xi, 0.05), quantile(knot_draws$gpd_xi, 0.95)))
cat(sprintf(fmt, "HSGP:",
  mean(draws$rl100), quantile(draws$rl100, 0.05), quantile(draws$rl100, 0.95),
  mean(draws$gpd_xi), quantile(draws$gpd_xi, 0.05), quantile(draws$gpd_xi, 0.95)))

# ---- 8. Three-way return level plot ----
cat("\nGenerating comparison plot...\n")

rp_grid <- exp(seq(log(2), log(200), length.out = 100))

# HSGP return levels
rl_hsgp <- matrix(NA, nrow(draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  rl_hsgp[, j] <- draws$gev_mu + draws$gev_sigma / draws$gev_xi *
    ((-log(1 - 1/T))^(-draws$gev_xi) - 1)
}

# Knot GP return levels
rl_knot <- matrix(NA, nrow(knot_draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  rl_knot[, j] <- knot_draws$gev_mu + knot_draws$gev_sigma / knot_draws$gev_xi *
    ((-log(1 - 1/T))^(-knot_draws$gev_xi) - 1)
}

# M&S return levels
rl_ms <- matrix(NA, length(psi_draws), length(rp_grid))
for (j in seq_along(rp_grid)) {
  T <- rp_grid[j]
  mu_d <- exp(psi_draws)
  sig_d <- exp(psi_draws + tau_draws)
  xi_d <- g(phi_draws)
  rl_ms[, j] <- mu_d + sig_d / xi_d * ((-log(1 - 1/T))^(-xi_d) - 1)
}

rl_comp <- data.frame(
  rp = rep(rp_grid, 3),
  mean = c(colMeans(rl_hsgp), colMeans(rl_knot), colMeans(rl_ms)),
  lo = c(apply(rl_hsgp, 2, quantile, 0.05),
         apply(rl_knot, 2, quantile, 0.05),
         apply(rl_ms, 2, quantile, 0.05)),
  hi = c(apply(rl_hsgp, 2, quantile, 0.95),
         apply(rl_knot, 2, quantile, 0.95),
         apply(rl_ms, 2, quantile, 0.95)),
  model = rep(c("HSGP", "Knot GP", "Max-and-Smooth"), each = length(rp_grid))
)

am <- readRDS("data/annual_maxima_andalucia.rds")
mal_am <- am %>% filter(indicativo == "6155A") %>% arrange(max_prec)
n_am <- nrow(mal_am)
mal_am$rp <- rev((n_am + 0.12) / (seq_len(n_am) - 0.44 + 0.12))

p_rl <- ggplot(rl_comp, aes(x = rp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.12) +
  geom_line(aes(y = mean, colour = model), linewidth = 0.8) +
  geom_point(data = mal_am, aes(x = rp, y = max_prec),
             shape = 21, fill = "grey20", colour = "white", size = 2, stroke = 0.4) +
  scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100, 200),
                labels = c("2", "5", "10", "20", "50", "100", "200")) +
  scale_colour_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63",
                                  "Max-and-Smooth" = "#3366CC")) +
  scale_fill_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63",
                                "Max-and-Smooth" = "#3366CC")) +
  labs(title = "Return level comparison: HSGP vs Knot GP vs Max-and-Smooth",
       subtitle = sprintf("M\u00e1laga Aeropuerto | HSGP: J=%d harmonics, wall time: %.0fs", J, elapsed),
       x = "Return period (years)", y = "Daily rainfall (mm)",
       colour = NULL, fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# ---- 9. Seasonal curve comparison ----
f_occ_hsgp <- fit$draws("f_occ", format = "matrix")
f_occ_knot <- knot_res$fit$draws("f_occ", format = "matrix")

p_rain_hsgp <- 1 / (1 + exp(-f_occ_hsgp))
p_rain_knot <- 1 / (1 + exp(-f_occ_knot))

seasonal_df <- data.frame(
  doy = rep(1:365, 2),
  mean = c(colMeans(p_rain_hsgp), colMeans(p_rain_knot)),
  lo = c(apply(p_rain_hsgp, 2, quantile, 0.05),
         apply(p_rain_knot, 2, quantile, 0.05)),
  hi = c(apply(p_rain_hsgp, 2, quantile, 0.95),
         apply(p_rain_knot, 2, quantile, 0.95)),
  model = rep(c("HSGP", "Knot GP"), each = 365)
)

emp_monthly <- mal %>%
  mutate(month = as.integer(format(fecha, "%m"))) %>%
  group_by(month) %>%
  summarise(frac_wet = mean(prec_num > 0), .groups = "drop") %>%
  mutate(doy_mid = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349))

p_seas <- ggplot(seasonal_df, aes(x = doy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15) +
  geom_line(aes(y = mean, colour = model, linetype = model), linewidth = 0.8) +
  geom_point(data = emp_monthly, aes(x = doy_mid, y = frac_wet),
             colour = "red", size = 2.5) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  scale_colour_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63")) +
  scale_fill_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63")) +
  scale_linetype_manual(values = c("HSGP" = "solid", "Knot GP" = "dashed")) +
  labs(title = "Rainfall occurrence: HSGP vs Knot GP seasonal curves",
       subtitle = "Red dots = empirical monthly fractions",
       x = NULL, y = "P(rain > 0)", colour = NULL, fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# ---- 10. Intensity curve comparison ----
f_int_hsgp <- fit$draws("f_int", format = "matrix")
f_int_knot <- knot_res$fit$draws("f_int", format = "matrix")

gamma_hsgp <- draws$gamma_shape
gamma_knot <- knot_draws$gamma_shape

mean_rain_hsgp <- sweep(exp(f_int_hsgp), 1, gamma_hsgp, "*")
mean_rain_knot <- sweep(exp(f_int_knot), 1, gamma_knot, "*")

int_df <- data.frame(
  doy = rep(1:365, 2),
  mean = c(colMeans(mean_rain_hsgp), colMeans(mean_rain_knot)),
  lo = c(apply(mean_rain_hsgp, 2, quantile, 0.05),
         apply(mean_rain_knot, 2, quantile, 0.05)),
  hi = c(apply(mean_rain_hsgp, 2, quantile, 0.95),
         apply(mean_rain_knot, 2, quantile, 0.95)),
  model = rep(c("HSGP", "Knot GP"), each = 365)
)

emp_int <- mal %>%
  filter(prec_num > 0) %>%
  mutate(month = as.integer(format(fecha, "%m"))) %>%
  group_by(month) %>%
  summarise(mean_wet = mean(prec_num), .groups = "drop") %>%
  mutate(doy_mid = c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349))

p_int <- ggplot(int_df, aes(x = doy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = model), alpha = 0.15) +
  geom_line(aes(y = mean, colour = model, linetype = model), linewidth = 0.8) +
  geom_point(data = emp_int, aes(x = doy_mid, y = mean_wet),
             colour = "red", size = 2.5) +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  scale_colour_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63")) +
  scale_fill_manual(values = c("HSGP" = "#4CAF50", "Knot GP" = "#E91E63")) +
  scale_linetype_manual(values = c("HSGP" = "solid", "Knot GP" = "dashed")) +
  labs(title = "Mean wet-day rainfall: HSGP vs Knot GP seasonal curves",
       subtitle = "Red dots = empirical monthly mean",
       x = NULL, y = "Mean rainfall (mm) given rain", colour = NULL, fill = NULL, linetype = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

p_all <- p_seas / p_int / p_rl
ggsave("figures/presentation/daily_hsgp_comparison.png",
       p_all, width = 12, height = 14, dpi = 200, bg = "white")

cat("Saved figures/presentation/daily_hsgp_comparison.png\n")

saveRDS(list(fit = fit, stan_data = stan_data, threshold = u, J = J, elapsed = elapsed),
        "data/daily_hsgp_malaga.rds")
cat("Saved data/daily_hsgp_malaga.rds\n")

cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
