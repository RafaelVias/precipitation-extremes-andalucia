# Spatial Extreme Precipitation in Andalucía

Bayesian spatial extreme value analysis of daily precipitation across Andalucía, Spain, using the **Max-and-Smooth** two-stage framework with Matérn(5/2) Gaussian process spatial smoothing and penalised complexity (PC) priors.

## Method

The analysis follows the Max-and-Smooth approach of Hrafnkelsson et al. (2021) as presented in Hazra, Huser & Jóhannesson (2023, Ch. 7):

1. **Stage 1** — Fit a Poisson point process (PPP) GEV model independently at each station via maximum likelihood. Parametric bootstrap provides per-station covariance matrices.

2. **Stage 2** — Smooth the Stage 1 estimates spatially using a Matérn(5/2) Gaussian process in Stan. The GEV parameters (ψ = log μ, τ = log σ/μ, φ = h(ξ)) are modelled as independent GPs with PC priors (Fuglstad et al., 2019) on the range and marginal variance.

Predictions at unobserved locations are obtained by kriging from the posterior GP.

## Data

- **Source**: AEMET (Agencia Estatal de Meteorología) OpenData API
- **Coverage**: 127 stations across all 8 provinces of Andalucía
- **Period**: 1950–2024 (varies by station)
- **QC**: Station-years require ≥90% daily completeness (≥330 days/year)

## Repository structure

```
R/
├── 01_acquire_aemet_data.R        # Download daily data from AEMET API
├── 02_exploratory_analysis.R      # EDA and station selection
├── 03_max_and_smooth.R            # Full M&S pipeline (Stage 1 + INLA Stage 2)
├── 03d_smooth_matern_pc.R         # Stage 2: Matérn GP + PC priors (Stan)
├── 05b_interpolation_maps_pc.R    # Spatial prediction maps (μ, σ, ξ, RL100)
├── 06a_return_level_panels.R      # Return level panels (mean + uncertainty)
├── 06b_exceedance_panels.R        # Exceedance probability maps
├── 06c_return_level_1x3.R         # 2×2 return level maps with AEMET thresholds
├── 07a_station_return_level_curves.R  # Station-level return level curves
├── 07b_station_timeseries.R       # Annual maxima time series
├── 07c_seasonal_rainfall.R        # Seasonal rainfall patterns
├── 08a_daily_model_malaga.R       # Daily hurdle + GPD model (Málaga)
├── 08b_daily_gp_malaga.R         # Daily GP model (Málaga)
├── 08c_daily_hsgp_malaga.R       # Daily HSGP model (Málaga)
└── 09_spanish_figures.R           # Spanish-language figure variants

Stan/
├── smooth_matern_pc.stan          # Stage 2: Matérn(5/2) GP + PC priors
├── daily_hsgp_gpd.stan            # Hilbert-space GP daily model
├── daily_gp_gpd.stan              # Full GP daily model
└── daily_hurdle_gpd.stan          # Hurdle + GPD daily model

data/
├── daily_precip_andalucia_raw.rds # Raw daily precipitation (AEMET)
├── annual_maxima_andalucia.rds    # Annual block maxima (QC-filtered)
├── stations_andalucia.rds/.csv    # Station metadata and coordinates
├── stage1_results.rds             # Per-station PPP GEV MLEs + covariances
└── stage2_matern_pc_results.rds   # Posterior draws from spatial GP

figures/
└── jesus-figures/                 # Publication-quality figures
```

## Reproducing the analysis

### Prerequisites

```r
install.packages(c("dplyr", "lubridate", "sf", "ggplot2", "patchwork",
                    "rnaturalearth", "rnaturalearthdata", "cmdstanr"))
```

[CmdStan](https://mc-stan.org/cmdstanr/) must be installed separately.

### Pipeline

Run scripts in numerical order from the project root:

```bash
Rscript R/01_acquire_aemet_data.R      # Requires AEMET API key
Rscript R/02_exploratory_analysis.R
Rscript R/03_max_and_smooth.R          # Stage 1 (~10 min)
Rscript R/03d_smooth_matern_pc.R       # Stage 2 (~18 min, 4 chains)
Rscript R/05b_interpolation_maps_pc.R  # Spatial predictions (~5 min)
Rscript R/06a_return_level_panels.R
Rscript R/06b_exceedance_panels.R
Rscript R/06c_return_level_1x3.R
Rscript R/07a_station_return_level_curves.R
```

The `data/` directory includes pre-computed results so that figure scripts (05–09) can be run without re-fitting the models.

## Key results

- **100-year return levels** range from ~80 mm/day in arid eastern Almería to ~200 mm/day on the western Málaga coast.
- The spatial GP with PC priors produces smoother fields than the BYM2 alternative, with longer estimated correlation ranges (ρ ≈ 0.7–1.7° for the GEV parameters).
- The shape parameter ξ is near zero across the region (Gumbel-like tails), with slight positive values along the coast.

## References

- Hrafnkelsson, B., Siegert, S., Huser, R., Bakka, H. & Jóhannesson, Á. V. (2021). Max-and-Smooth: a two-step approach for approximate Bayesian inference in latent Gaussian models. *Bayesian Analysis*, 16(2), 611–638. [doi:10.1214/20-BA1219](https://doi.org/10.1214/20-BA1219)

- Hazra, A., Huser, R. & Jóhannesson, Á. V. (2023). Bayesian spatial modelling of extreme precipitation return levels. In: Hrafnkelsson, B. (ed.) *Bayesian Latent Gaussian Models*. Chapman & Hall/CRC, Ch. 7. [doi:10.1201/9781003emi6985-7](https://doi.org/10.1201/9781003emi6985-7)

- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing priors that penalize the complexity of Gaussian random fields. *Journal of the American Statistical Association*, 114(525), 445–452. [doi:10.1080/01621459.2017.1415907](https://doi.org/10.1080/01621459.2017.1415907)

- Carpenter, B. et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1). [doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)
