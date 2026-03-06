# Spatial Extreme Precipitation in Andalucía

Bayesian spatial extreme value analysis of daily precipitation across Andalucía, Spain, using the **Max-and-Smooth** two-stage framework with Matérn(5/2) Gaussian process spatial smoothing and penalised complexity (PC) priors.

## Method

The analysis follows the Max-and-Smooth approach of Hrafnkelsson et al. (2021) as presented in Hazra, Huser & Jóhannesson (2023, Ch. 7):

1. **Stage 1** — Fit a Poisson point process (PPP) GEV model independently at each station via maximum likelihood. Parametric bootstrap provides per-station covariance matrices.

2. **Stage 2** — Smooth the Stage 1 estimates spatially using a Matérn(5/2) Gaussian process in Stan. The GEV parameters (ψ = log μ, τ = log σ/μ, φ = h(ξ)) are modelled as independent GPs with PC priors (Fuglstad et al., 2019) on the range and marginal variance.

Predictions at unobserved locations are obtained via kriging from the posterior GP.

## Data

- **Source**: AEMET (Agencia Estatal de Meteorología) OpenData API
- **Coverage**: 127 stations across all 8 provinces of Andalucía
- **Period**: 1950–2024 (varies by station)
- **QC**: Station-years require ≥90% daily completeness (≥330 days/year)

## Core pipeline

```bash
Rscript R/03_max_and_smooth.R      # Stage 1: per-station PPP GEV MLEs (~10 min)
Rscript R/03d_smooth_matern_pc.R   # Stage 2: spatial GP smoothing in Stan (~18 min)
```

Pre-computed results are in `data/stage1_results.rds` and `data/stage2_matern_pc_results.rds`.

## References

- Hrafnkelsson, B., Siegert, S., Huser, R., Bakka, H. & Jóhannesson, Á. V. (2021). Max-and-Smooth: a two-step approach for approximate Bayesian inference in latent Gaussian models. *Bayesian Analysis*, 16(2), 611–638. [doi:10.1214/20-BA1219](https://doi.org/10.1214/20-BA1219)

- Hazra, A., Huser, R. & Jóhannesson, Á. V. (2023). Bayesian spatial modelling of extreme precipitation return levels. In: Hrafnkelsson, B. (ed.) *Bayesian Latent Gaussian Models*. Chapman & Hall/CRC, Ch. 7.

- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing priors that penalize the complexity of Gaussian random fields. *Journal of the American Statistical Association*, 114(525), 445–452. [doi:10.1080/01621459.2017.1415907](https://doi.org/10.1080/01621459.2017.1415907)

- Carpenter, B. et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1). [doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)
