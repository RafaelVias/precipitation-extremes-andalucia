# 00_station_map.R — Station network overview map
#
# Map of the 127 AEMET stations used in the analysis, coloured by observed
# maximum precipitation and sized by record length. Province boundaries shown.
#
# Run from project root: Rscript R/00_station_map.R

library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)

Sys.setlocale("LC_ALL", "en_US.UTF-8")

cat("========================================\n")
cat("Station network map\n")
cat("========================================\n")

# ---- 1. Load data ----
cat("Loading data...\n")

am <- readRDS("data/annual_maxima_andalucia.rds")

# Exclude non-mainland stations (Alborán island, removed in 02_stage1_mle.R)
am <- am %>% filter(indicativo != "6381")

# Per-station summaries
stn_summary <- am %>%
  group_by(indicativo) %>%
  summarise(
    n_years  = n(),
    max_prec = max(max_prec, na.rm = TRUE),
    nombre   = first(nombre),
    provincia = first(provincia),
    lon      = first(longitud),
    lat      = first(latitud),
    .groups  = "drop"
  )

cat("  Stations:", nrow(stn_summary), "\n")
cat("  Record lengths:", min(stn_summary$n_years), "-",
    max(stn_summary$n_years), "years\n")
cat("  Max observed:", round(max(stn_summary$max_prec), 1), "mm\n")

# ---- 2. Geographic context ----
cat("Loading boundaries...\n")

states <- ne_states(country = "Spain", returnclass = "sf")
andalucia_provs <- states[grep("Andaluc", states$region), ]
andalucia <- st_union(andalucia_provs)

neighbours <- states[states$name %in% c("Murcia", "Albacete", "Ciudad Real", "Badajoz"), ]
portugal <- ne_countries(country = "Portugal", scale = "medium", returnclass = "sf")
morocco  <- ne_countries(country = "Morocco", scale = "medium", returnclass = "sf")
bbox_plot <- st_bbox(c(xmin = -8.0, xmax = -1.3, ymin = 35.8, ymax = 38.8), crs = 4326)
port_crop <- st_crop(portugal, bbox_plot)
mor_crop  <- st_crop(morocco, bbox_plot)

# Province label positions (centroids, manually nudged)
prov_labels <- data.frame(
  name = c("Huelva", "Sevilla", "C\u00f3rdoba", "Ja\u00e9n",
           "C\u00e1diz", "M\u00e1laga", "Granada", "Almer\u00eda"),
  lon  = c(-6.95, -5.85, -4.78, -3.75,
           -5.85, -4.55, -3.45, -2.40),
  lat  = c(37.65, 37.50, 37.95, 37.85,
           36.55, 36.80, 37.25, 37.05)
)

# ---- 3. Plot ----
cat("Generating map...\n")

p <- ggplot() +
  # Background countries/regions
  geom_sf(data = port_crop, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = mor_crop, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_sf(data = neighbours, fill = "grey92", colour = "grey70", linewidth = 0.2) +
  # Province boundaries
  geom_sf(data = andalucia_provs, fill = "grey98", colour = "grey60", linewidth = 0.3) +
  # Province labels
  geom_text(data = prov_labels, aes(x = lon, y = lat, label = name),
            size = 3.0, colour = "grey50", fontface = "italic") +
  # Station points
  geom_point(data = stn_summary,
             aes(x = lon, y = lat, fill = max_prec, size = n_years),
             shape = 21, colour = "grey30", stroke = 0.3) +
  scale_fill_viridis_c(option = "B", name = "Observed\nmax (mm)",
                       breaks = c(50, 100, 150, 200, 250)) +
  scale_size_continuous(name = "Record\nlength (yr)",
                        range = c(1.5, 5),
                        breaks = c(10, 20, 30, 40, 50, 60, 70)) +
  # Andalucía outline
  geom_sf(data = andalucia, fill = NA, colour = "grey30", linewidth = 0.5) +
  coord_sf(xlim = c(-7.8, -1.4), ylim = c(35.9, 38.8)) +
  labs(
    title = "AEMET station network across Andaluc\u00eda",
    subtitle = sprintf("%d stations | %d\u2013%d years of record | QC: \u226590%% daily completeness",
                       nrow(stn_summary), min(am$year), max(am$year))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave("figures/station_map.png", p, width = 12, height = 7, dpi = 200, bg = "white")

cat("Saved figures/station_map.png\n")
cat("========================================\n")
cat("Done.\n")
cat("========================================\n")
