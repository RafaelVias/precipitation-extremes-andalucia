# 05b_plot_coarse.R — Quick plot at model resolution (no IDW)
# Just re-runs 05_interpolation_maps.R but skips the fine-grid IDW step
# and plots directly at the 0.1° model grid resolution.

# Source the main script up to the prediction extraction (steps 1-9),
# then plot at coarse resolution.
source("vendor/max_and_smooth/stage1_functions.R")

library(cmdstanr)
library(sf)
library(ggplot2)
library(patchwork)
library(rnaturalearth)
library(Matrix)
library(deldir)
library(INLA)
library(FNN)

s1           <- readRDS("data/stage1_results.rds")
mles.covmats <- s1$mles.covmats
loc          <- s1$loc
meta         <- s1$station_meta
ns           <- nrow(loc)

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

pred_res <- 0.1
bbox <- st_bbox(andalucia)
pred_pts <- expand.grid(
  lon = seq(bbox["xmin"] + pred_res / 2, bbox["xmax"], by = pred_res),
  lat = seq(bbox["ymin"] + pred_res / 2, bbox["ymax"], by = pred_res)
)
pred_sf <- st_as_sf(pred_pts, coords = c("lon", "lat"), crs = 4326)
inside <- st_intersects(pred_sf, andalucia, sparse = FALSE)[, 1]
pred_pts <- pred_pts[inside, ]
n_pred <- nrow(pred_pts)

all_locs <- rbind(
  data.frame(lon = loc[, "lon"], lat = loc[, "lat"]),
  pred_pts
)
n_total <- nrow(all_locs)

cat("Nodes:", n_total, "(", ns, "stations +", n_pred, "grid)\n")

# Build graph
tri <- deldir(all_locs$lon, all_locs$lat)
del_edges <- tri$delsgs[, c("ind1", "ind2")]
edges_upper <- unique(data.frame(
  node1 = pmin(del_edges$ind1, del_edges$ind2),
  node2 = pmax(del_edges$ind1, del_edges$ind2)
))
edges_all <- rbind(edges_upper, data.frame(node1 = edges_upper$node2, node2 = edges_upper$node1))

adj_mat <- sparseMatrix(i = edges_upper$node1, j = edges_upper$node2,
                        x = 1, dims = c(n_total, n_total), symmetric = TRUE)
Q <- Diagonal(n_total, rowSums(adj_mat)) - adj_mat
Q_pert <- Q + Diagonal(n_total) * max(diag(Q)) * sqrt(.Machine$double.eps)
Q_inv <- inla.qinv(Q_pert, constr = list(A = matrix(1, 1, n_total), e = 0))
scaling_factor <- exp(mean(log(diag(Q_inv))))

# Precision matrix (stations only)
psi_hat <- sapply(mles.covmats, function(x) x$mle[1])
tau_hat <- sapply(mles.covmats, function(x) x$mle[2])
phi_hat <- sapply(mles.covmats, function(x) x$mle[3])
eta_hat <- c(psi_hat, tau_hat, phi_hat)

n_lik <- 3 * ns
rows <- cols <- vals <- c()
for (i in seq_len(ns)) {
  Sigma_i <- mles.covmats[[i]]$covmat
  ev <- eigen(Sigma_i, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) / max(ev) < 1e-10 || min(ev) <= 0)
    Sigma_i <- Sigma_i + diag(3) * max(ev) * 1e-6
  Q_i <- solve(Sigma_i)
  for (r in 1:3) for (cc in 1:3) {
    rows <- c(rows, i + (r - 1) * ns)
    cols <- c(cols, i + (cc - 1) * ns)
    vals <- c(vals, Q_i[r, cc])
  }
}
Q_full <- sparseMatrix(i = rows, j = cols, x = vals, dims = c(n_lik, n_lik))
L <- t(chol(Q_full))

stan_data <- list(
  n_stations = ns, n_pred = n_pred, n_param = 3L, eta_hat = eta_hat,
  n_edges = nrow(edges_all), node1 = edges_all$node1, node2 = edges_all$node2,
  scaling_factor = scaling_factor,
  n_nonzero_chol_Q = length(L@x), n_values = as.integer(diff(L@p)),
  index = as.integer(L@i + 1L), value = L@x, log_det_Q = sum(log(diag(L)))
)

# Inits
mu_psi <- mean(psi_hat); mu_tau <- mean(tau_hat); mu_phi <- mean(phi_hat)
sd_psi <- sd(psi_hat); sd_tau <- sd(tau_hat); sd_phi <- max(sd(phi_hat), 1e-4)
psi_raw <- (psi_hat - mu_psi) / sd_psi
tau_raw <- (tau_hat - mu_tau) / sd_tau
phi_raw <- (phi_hat - mu_phi) / sd_phi

nn_init <- get.knnx(as.matrix(loc), as.matrix(pred_pts), k = 3)
eta_spatial_init <- rbind(
  cbind(psi_raw, tau_raw, phi_raw),
  cbind(rowMeans(matrix(psi_raw[nn_init$nn.index], nrow = n_pred)),
        rowMeans(matrix(tau_raw[nn_init$nn.index], nrow = n_pred)),
        rowMeans(matrix(phi_raw[nn_init$nn.index], nrow = n_pred)))
)
inits <- list(
  mu = c(mu_psi, mu_tau, mu_phi),
  sigma = c(sd_psi, sd_tau, max(sd_phi, 0.01)),
  rho = c(0.5, 0.5, 0.5),
  eta_spatial = eta_spatial_init,
  eta_random = matrix(0, n_total, 3)
)

# Run Stan
model <- cmdstan_model("Stan/smooth_bym2_predict.stan")
cat("Running Stan...\n")
fit <- model$sample(
  data = stan_data, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 1000, refresh = 500,
  init = rep(list(inits), 4)
)

# Extract
draws <- fit$draws(format = "draws_matrix")
n_draws <- nrow(draws)
M <- 100

psi_stn <- as.matrix(draws[, paste0("eta[", 1:ns, ",1]")])
tau_stn <- as.matrix(draws[, paste0("eta[", 1:ns, ",2]")])
phi_stn <- as.matrix(draws[, paste0("eta[", 1:ns, ",3]")])

psi_grid <- as.matrix(draws[, paste0("eta_pred[", 1:n_pred, ",1]")])
tau_grid <- as.matrix(draws[, paste0("eta_pred[", 1:n_pred, ",2]")])
phi_grid <- as.matrix(draws[, paste0("eta_pred[", 1:n_pred, ",3]")])

plot_df <- rbind(
  data.frame(
    lon = loc[, "lon"], lat = loc[, "lat"],
    mu_mean = colMeans(exp(psi_stn)),
    sigma_mean = colMeans(exp(psi_stn + tau_stn)),
    xi_mean = colMeans(g(phi_stn)),
    rl100_mean = colMeans(exp(psi_stn) + exp(psi_stn + tau_stn) / g(phi_stn) *
                            ((-log(1 - 1 / M))^(-g(phi_stn)) - 1))
  ),
  data.frame(
    lon = pred_pts$lon, lat = pred_pts$lat,
    mu_mean = colMeans(exp(psi_grid)),
    sigma_mean = colMeans(exp(psi_grid + tau_grid)),
    xi_mean = colMeans(g(phi_grid)),
    rl100_mean = colMeans(exp(psi_grid) + exp(psi_grid + tau_grid) / g(phi_grid) *
                            ((-log(1 - 1 / M))^(-g(phi_grid)) - 1))
  )
)

# Plot at coarse resolution
cat("Plotting at", pred_res, "deg tiles (no IDW)...\n")

stn_sf <- st_as_sf(meta, coords = c("longitud", "latitud"), crs = 4326)
neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

tile_res <- pred_res
base_theme <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

make_panel <- function(data, fill_var, fill_label, title, option = "C") {
  ggplot() +
    geom_sf(data = port_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = mor_crop, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_sf(data = neighbours, fill = "grey90", colour = "grey70", linewidth = 0.2) +
    geom_tile(data = data, aes(x = lon, y = lat, fill = .data[[fill_var]]),
              width = tile_res, height = tile_res) +
    scale_fill_viridis_c(option = option, name = fill_label) +
    geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.4) +
    geom_sf(data = stn_sf, size = 0.5, colour = "white", shape = 16) +
    labs(title = title) +
    coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
    base_theme
}

p1 <- make_panel(plot_df, "mu_mean", "mm", expression(hat(mu) ~ "(location)"), "C")
p2 <- make_panel(plot_df, "sigma_mean", "mm", expression(hat(sigma) ~ "(scale)"), "C")
p3 <- make_panel(plot_df, "xi_mean", expression(xi), expression(hat(xi) ~ "(shape)"), "D")
p4 <- make_panel(plot_df, "rl100_mean", "mm", "100-year return level", "B")

p_all <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Posterior predictive: extreme precipitation (model grid, no IDW)",
    subtitle = sprintf("BYM2 + ICAR | %d stations + %d grid = %d nodes | %.2f deg tiles",
                       ns, n_pred, n_total, pred_res),
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

ggsave("figures/interpolated_maps.pdf", p_all, width = 14, height = 10)
ggsave("figures/interpolated_maps.png", p_all, width = 14, height = 10, dpi = 200)
cat("Saved (no IDW, coarse tiles)\n")
