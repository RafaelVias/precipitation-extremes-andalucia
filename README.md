🇬🇧 English | [🇪🇸 Español](README.es.md)

# Spatial Extreme Precipitation in Andalucía

Bayesian spatial extreme value analysis of daily precipitation across Andalucía, Spain, using the **Max-and-Smooth** two-stage framework with Matérn(5/2) Gaussian process spatial smoothing and penalised complexity (PC) priors.

## Introduction

Extreme daily precipitation is a key input for flood risk assessment, infrastructure design, and water resource planning. Estimating how rare a given rainfall event is — and how that rarity varies across space — requires fitting statistical models to observed data. When individual station records are short or sparse, site-by-site estimates can be noisy and spatially inconsistent.

The goal of this project is to provide spatially continuous estimates of **return levels** (the rainfall amount expected to be exceeded on average once every *T* years) and **exceedance probabilities** (the chance of exceeding a critical threshold within a planning horizon) across Andalucía, to support flood risk assessment and hydrological planning in the region. Using daily precipitation records from 127 AEMET stations, the two-stage **Max-and-Smooth** framework (Hrafnkelsson et al., 2021) fits a GEV distribution at each station independently, then borrows strength across stations through a Matérn Gaussian process prior, producing spatially coherent maps with full uncertainty quantification.

## Data

- **Source**: AEMET (Agencia Estatal de Meteorología) OpenData API
- **Coverage**: 127 stations across all 8 provinces of Andalucía
- **Period**: 1950–2024 (varies by station; minimum record length governed by Bayesian borrowing of strength)
- **QC**: Station-years require ≥90% daily completeness (≥330 days/year)

![Station network](figures/station_map.png)

## Method

The analysis follows the **Max-and-Smooth** approach of Hrafnkelsson et al. (2021) as presented in Hazra, Huser & Jóhannesson (2023, Ch. 7). Annual maximum daily rainfall at each station is modelled by a generalised extreme value (GEV) distribution whose parameters vary smoothly across space through Gaussian process priors.

### GEV model and reparametrisation

The GEV distribution function is

$$F(z) = \exp\left[-\left(1 + \xi \left(\frac{z - \mu}{\sigma}\right)\right)^{-1/\xi}\right]$$

with location $\mu$, scale $\sigma > 0$, and shape $\xi$. The parameters are reparametrised as

$$(\psi, \tau, \phi) = (\log(\mu), \log(\sigma / \mu), h(\xi))$$

where $h$ is a smooth bijection mapping $\xi \in (-0.5, 0.5)$ to the real line:

$$h(\xi) = a + b \cdot \log(-\log(1 - u^c)),\quad u = \frac{\xi - \xi_{\min}}{\xi_{\max} - \xi_{\min}}$$

with constants $c = 0.8$, $a$ and $b$ chosen so that $h(0) \approx 0$. This ensures all three working parameters are unbounded, which is essential for the Gaussian surrogate likelihood and the GP prior in Stage 2.

### Stage 1 — Site-wise maximum likelihood

At each station, a Poisson point process (PPP) likelihood is maximised over daily observations exceeding a site-specific threshold (the 75th percentile of positive precipitation). Multi-start optimisation over five initial shape values avoids local optima, with interior solutions ($\xi$ away from the bounds) preferred over boundary solutions to ensure reliable Hessians.

A parametric bootstrap (1000 replicates per station) provides per-station 3 × 3 covariance matrices $\hat{\Sigma}_i$ for the surrogate Gaussian likelihood used in Stage 2. Exceedances are simulated from the fitted PPP model and re-estimated, yielding empirical covariances of the MLE across replicates.

### Stage 2 — Spatial smoothing

The Stage 1 estimates $\hat{\eta}_i = (\hat{\psi}_i, \hat{\tau}_i, \hat{\phi}_i)$ are treated as noisy observations of a latent spatial field. For each of the three GEV parameters, an independent Gaussian process prior is placed on the latent values:

$$\eta_k \sim \text{GP}(\mu_k, C_k)$$

with a Matérn covariance of smoothness 5/2:

$$C_k(d) = \sigma_k^2 \left(1 + \frac{\sqrt{5} d}{\rho_k} + \frac{5 d^2}{3 \rho_k^2}\right) \exp\left(-\frac{\sqrt{5} d}{\rho_k}\right)$$

where $d$ is inter-station distance (in degrees), $\sigma_k$ is the marginal standard deviation, and $\rho_k$ is the practical range. A nugget term $\nu_k^2$ is added to the diagonal to absorb station-specific noise not captured by the spatial GP.

#### Surrogate likelihood

The Stage 1 bootstrap covariances form a block-diagonal precision matrix $Q$ with 3 × 3 blocks $\hat{\Sigma}_i^{-1}$. The surrogate likelihood is

$$\hat{\eta} \mid \eta \sim \mathcal{N}(\eta, Q^{-1})$$

which decouples the computationally expensive per-station PPP fits from the spatial smoothing.

#### Non-centered parameterisation

To improve HMC sampling efficiency, the GP is parameterised as $\eta_k = \mu_k + L_k z_k$, where $L_k$ is the Cholesky factor of the covariance matrix and $z_k \sim \mathcal{N}(0, I)$.

#### Penalised complexity priors

Following Fuglstad et al. (2019), the GP hyperparameters receive PC priors calibrated as:

| Parameter | Prior | Calibration |
|-----------|-------|-------------|
| $\sigma_k$ (marginal SD) | Exponential | $P(\sigma > 1) = 0.05 \Rightarrow \lambda_\sigma = 3.0$ |
| $\rho_k$ (range) | PC prior on range | $P(\rho < 0.1°) = 0.05 \Rightarrow \lambda_\rho = 0.30$ |
| $\nu_k$ (nugget SD) | Exponential | $P(\nu > 0.5) = 0.05 \Rightarrow \lambda_\nu = 6.0$ |

The intercepts receive vague priors: $\mu_k \sim \mathcal{N}(0, 100^2)$. The model is fitted jointly in Stan (Carpenter et al., 2017) using NUTS with 4 chains × 2000 iterations (1000 warmup), `adapt_delta = 0.9`.

### Prediction

Return levels and exceedance probabilities at unobserved locations are obtained from the **conditional distribution of the Matérn GP**. The GP hyperparameters $(\hat{\sigma}_k, \hat{\rho}_k, \hat{\nu}_k)$ are fixed at their posterior means to compute the conditional moments at a prediction location $s^*$:

$$\eta_k(s^*) \mid \eta_k \sim \mathcal{N}(\mu_{\text{cond}}, \sigma^2_{\text{cond}})$$

where

$$\mu_{\text{cond}} = \mu_k + \gamma_k^\top \Sigma_k^{-1} (\eta_k - \mu_k), \qquad \sigma^2_{\text{cond}} = \sigma_k^2 - \gamma_k^\top \Sigma_k^{-1} \gamma_k$$

with $\gamma_k$ the cross-covariance vector between $s^*$ and the stations, and $\Sigma_k$ the station-station covariance matrix (including nugget). For each posterior draw of the latent field $\eta_k$, a prediction is **sampled** from this conditional distribution, so the interpolation uncertainty $\sigma^2_{\text{cond}}$ is fully propagated into the posterior predictive return levels. Predictions are computed on a 0.025° grid across Andalucía.

## Results

### Return level maps

Posterior mean return levels at *T* = 10, 20, 50, and 100 years, interpolated across Andalucía on a 0.025° grid. White contour lines mark AEMET alarm thresholds (80 and 120 mm).

![Return level maps](figures/return_level_maps.png)

Posterior standard deviation of the return level estimates, reflecting uncertainty from both the GEV parameter estimation and the spatial interpolation.

![Return level SD](figures/return_level_maps_sd.png)

### Exceedance probability

Posterior mean probability that the annual maximum daily rainfall exceeds a given threshold (100, 150, 200 mm) at least once within a planning horizon (20, 50, 100 years).

![Exceedance probability](figures/exceedance_prob.png)

Posterior standard deviation of the exceedance probabilities, capturing uncertainty in both the tail behaviour and the spatial prediction.

![Exceedance probability SD](figures/exceedance_prob_sd.png)

### Station diagnostics

Return level curves at 6 selected stations. Points show observed annual maxima (Gringorten plotting positions). The dashed red line is the site-only MLE fit; the solid blue line is the spatially smoothed posterior mean with 90% credible band.

![Station return level curves](figures/station_return_level_curves.png)

## Dependencies

- **R** (≥ 4.2)
- **Stan**: [CmdStan](https://mc-stan.org/cmdstanr/) (≥ 2.33)
- **R packages**: cmdstanr, bayesplot, ggplot2, patchwork, sf, rnaturalearth, rnaturalearthdata, dplyr, lubridate, Matrix, climaemet

## Pipeline

Pre-computed results are in `data/stage1_results.rds` and `data/stage2_matern_pc_results.rds` (gitignored; regenerate with steps 2–3 below).

```bash
Rscript R/00_station_map.R          # Station network map (Figure 0)
Rscript R/01_acquire_data.R         # Download AEMET daily precipitation
Rscript R/02_stage1_mle.R           # Stage 1: per-station PPP GEV MLEs (~10 min)
Rscript R/03_stage2_smooth.R        # Stage 2: spatial GP smoothing in Stan (~18 min)
Rscript R/04_return_level_maps.R    # Return level maps (Figure 1)
Rscript R/05_exceedance_maps.R      # Exceedance probability maps (Figure 2)
Rscript R/06_station_diagnostics.R  # Station return level curves (Figure 3)
Rscript R/07_convergence_diagnostics.R  # MCMC convergence diagnostics (diagnostics/)
Rscript R/08_spanish_figures.R      # Figuras en español (figures/es/)
```

## References

- Hrafnkelsson, B., Siegert, S., Huser, R., Bakka, H. & Jóhannesson, Á. V. (2021). Max-and-Smooth: a two-step approach for approximate Bayesian inference in latent Gaussian models. *Bayesian Analysis*, 16(2), 611–638. [doi:10.1214/20-BA1219](https://doi.org/10.1214/20-BA1219)

- Hazra, A., Huser, R. & Jóhannesson, Á. V. (2023). Bayesian spatial modelling of extreme precipitation return levels. In: Hrafnkelsson, B. (ed.) *Bayesian Latent Gaussian Models*. Chapman & Hall/CRC, Ch. 7. [doi:10.1007/978-3-031-39791-2_7](https://doi.org/10.1007/978-3-031-39791-2_7)

- Fuglstad, G.-A., Simpson, D., Lindgren, F. & Rue, H. (2019). Constructing priors that penalize the complexity of Gaussian random fields. *Journal of the American Statistical Association*, 114(525), 445–452. [doi:10.1080/01621459.2017.1415907](https://doi.org/10.1080/01621459.2017.1415907)

- Carpenter, B. et al. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1). [doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

## Acknowledgements

The Stage 2 Stan implementation — in particular the sparse Cholesky surrogate likelihood (`normal_prec_chol_lpdf`) and the non-centered spatial parameterisation — is adapted from Brynjólfur Gauti Guðmundsson's [maxandsmooth](https://github.com/bgautijonsson/maxandsmooth) R package. The Stage 1 fitting code (`vendor/max_and_smooth/`) is from the companion repository to Hazra, Huser & Jóhannesson (2023): [arnabstatswithR/max_and_smooth](https://github.com/arnabstatswithR/max_and_smooth).
