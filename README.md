# Spatial Extreme Precipitation in Andalucía

Bayesian spatial extreme value analysis of daily precipitation across Andalucía, Spain, using the **Max-and-Smooth** two-stage framework with Matérn(5/2) Gaussian process spatial smoothing and penalised complexity (PC) priors.

## Method

The analysis follows the Max-and-Smooth approach of Hrafnkelsson et al. (2021) as presented in Hazra, Huser & Jóhannesson (2023, Ch. 7):

1. **Stage 1** — Fit a Poisson point process (PPP) GEV model independently at each station via maximum likelihood. The GEV parameters are reparametrised as (ψ, τ, φ) = (log μ, log σ/μ, h(ξ)), where h maps the shape parameter to an unbounded scale. Parametric bootstrap (1000 replicates per station) provides per-station covariance matrices for the surrogate likelihood used in Stage 2.

2. **Stage 2** — Smooth the Stage 1 estimates spatially using three independent Matérn(5/2) Gaussian processes (one per parameter), fitted jointly in Stan via Hamiltonian Monte Carlo. The model uses penalised complexity (PC) priors on the GP range and marginal standard deviation (Fuglstad et al., 2019), with a nugget term absorbing station-specific noise.

Predictions at unobserved locations are obtained via simple kriging from the posterior GP, using the posterior mean hyperparameters to compute kriging weights.

## Data

- **Source**: AEMET (Agencia Estatal de Meteorología) OpenData API
- **Coverage**: 127 stations across all 8 provinces of Andalucía
- **Period**: 1950–2024 (varies by station; minimum record length governed by Bayesian borrowing of strength)
- **QC**: Station-years require ≥90% daily completeness (≥330 days/year)

## Results

### Return level maps

Posterior mean return levels at *T* = 10, 20, 50, and 100 years, interpolated across Andalucía on a 0.025° grid. White contour lines mark AEMET alarm thresholds (80 and 120 mm).

![Return level maps](figures/return_level_maps.png)

### Exceedance probability

Posterior mean probability that the annual maximum daily rainfall exceeds a given threshold (100, 150, 200 mm) at least once within a planning horizon (20, 50, 100 years).

![Exceedance probability](figures/exceedance_prob.png)

### Station diagnostics

Return level curves at 6 selected stations. Points show observed annual maxima (Gringorten plotting positions). The dashed red line is the site-only MLE fit; the solid blue line is the spatially smoothed posterior mean with 90% credible band.

![Station return level curves](figures/station_return_level_curves.png)

## Pipeline

Pre-computed results are in `data/stage1_results.rds` and `data/stage2_matern_pc_results.rds` (gitignored; regenerate with steps 2–3 below).

```bash
Rscript R/01_acquire_data.R         # Download AEMET daily precipitation
Rscript R/02_stage1_mle.R           # Stage 1: per-station PPP GEV MLEs (~10 min)
Rscript R/03_stage2_smooth.R        # Stage 2: spatial GP smoothing in Stan (~18 min)
Rscript R/04_return_level_maps.R    # Return level maps (Figure 1)
Rscript R/05_exceedance_maps.R      # Exceedance probability maps (Figure 2)
Rscript R/06_station_diagnostics.R  # Station return level curves (Figure 3)
```

## References

- Hrafnkelsson, B., Siegert, S., Huser, R., Bakka, H. & Jóhannesson, Á. V. (2021). Max-and-Smooth: a two-step approach for approximate Bayesian inference in latent Gaussian models. *Bayesian Analysis*, 16(2), 611–638. [doi:10.1214/20-BA1219](https://doi.org/10.1214/20-BA1219)

- Hazra, A., Huser, R. & Jóhannesson, Á. V. (2023). Bayesian spatial modelling of extreme precipitation return levels. In: Hrafnkelsson, B. (ed.) *Bayesian Latent Gaussian Models*. Chapman & Hall/CRC, Ch. 7.

- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing priors that penalize the complexity of Gaussian random fields. *Journal of the American Statistical Association*, 114(525), 445–452. [doi:10.1080/01621459.2017.1415907](https://doi.org/10.1080/01621459.2017.1415907)

- Carpenter, B. et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1). [doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)
