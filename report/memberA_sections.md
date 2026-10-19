# Member A — report sections (draft)

> Draft text for the sections Member A owns: **Problem Statement**, **Data Description**, and the
> **Poisson-process model** (theory, estimation, simulation, validation). All numbers come from
> `notebooks/01_data_cleaning_eda.ipynb` and `notebooks/02_poisson_model.ipynb`; figures are in `figures/`.
> Rewrite in your own words before submitting.

---

## Problem Statement

Delhi has some of the worst urban air quality in the world. On a **Severe** air-quality day
(24-hour mean PM2.5 above 250 µg/m³, the top category of India's National AQI) the air is harmful even to
healthy people, schools are often closed, and emergency measures under the Graded Response Action Plan (GRAP)
are triggered. Whether a given day becomes Severe depends on many factors that cannot be predicted exactly:
wind speed, temperature inversions, crop-residue burning in neighbouring states, festival fireworks,
traffic and dust. The number and timing of Severe days is therefore naturally modelled as a **random
(stochastic) process**.

We ask:

1. At what rate do Severe days occur, and how uncertain is that rate?
2. Does the rate depend on the season, and has it changed between 2016 and 2025?
3. Are Severe days independent random events (a Poisson process), or do they come in clusters?
4. How well does a model fitted on past years predict the number of Severe days in later years?

## Data Description

**Source.** Central Pollution Control Board (CPCB), Central Control Room for Air Quality Management,
*CAAQMS Data Repository* (https://airquality.cpcb.gov.in/ccr/). Station: **R K Puram, Delhi – DPCC**
(site 124), a residential area in south Delhi. The data are hourly concentrations of 24 pollutants and
weather variables, one file per year, 2009–2026 (downloaded October 2026). We use **PM2.5 (µg/m³)**.

**Cleaning.**

| Step | Rule | Effect |
|---|---|---|
| Implausible values | readings at the analyser cap (1000 µg/m³) set to missing | 17 hours, all in 2011 |
| Daily mean | day is *observed* if ≥ 18 of 24 hours are valid (75 % rule) | |
| Short gaps | gaps of ≤ 3 consecutive days linearly interpolated | 194 days in the study window (5.3 %) |
| Long gaps | left missing, excluded from the exposure time | 223 days (6.1 %) |

**Study window.** 2009–2010 have no data, 2014 has none, 2013 has only 17 % and 2015 starts in April
(Figure A01). From **2016 to 2025** every year has at least 72 % observed days (2018–2025 at least 90 %),
so the models are estimated on these **10 years (3653 days, 3430 valid)**. January–September 2026 is kept
as an additional hold-out period.

**Summary (2016–2025, valid days).** Mean PM2.5 = 111.9 µg/m³, median 83.9, maximum 729.3. In terms
of the CPCB AQI categories: Good 11.6 %, Satisfactory 23.4 %, Moderate 18.1 %, Poor 12.4 %,
Very Poor 25.8 %, Severe 8.7 % of days. Pollution has a very strong annual cycle (Figures A02–A04):
monthly medians are lowest in July–August and highest in November–January.

**Event definition.** We define an *event* as a day with PM2.5 > *c*. Exceedances of the 60 µg/m³ national
24-hour standard happen on 65 % of days, and exceedances of 120 µg/m³ on 35 %. Neither is a rare event.
We use **c = 250 µg/m³ ("Severe")**: 300 days in 10 years (8.7 % of valid days). These form
**112 episodes** (runs of consecutive Severe days) with a mean length of 2.7 days.

## Poisson-Process Model

### Formulation

Let time *t* be measured in days. The counting process *N(t)* is the number of Severe days in (0, t], with
state space {0, 1, 2, …}. Other random variables are the inter-arrival times $X_k = S_k - S_{k-1}$ between
successive Severe days and the episode lengths $L_i$. Missing days are excluded: only valid days add to the
exposure time *T*.

* **Homogeneous Poisson process (HPP)** with rate λ: N(0) = 0, independent and stationary increments.
  Then $N(t)\sim\text{Poisson}(\lambda t)$, $E[N(t)]=\operatorname{Var}[N(t)]=\lambda t$, and
  $X_k \sim \text{Exp}(\lambda)$ i.i.d. with mean 1/λ.
* **Non-homogeneous Poisson process (NHPP)** with a piecewise-constant monthly rate $\lambda(t)=\lambda_{m(t)}$:
  $N(t)\sim\text{Poisson}(\Lambda(t))$, where $\Lambda(t)=\int_0^t\lambda(s)\,ds$.
* **Compound (clustered) NHPP**: episodes start as an NHPP with rate $\mu_{m(t)}$; each lasts
  $L_i\sim\text{Geometric}(1-q)$ days. $N(t)=\sum_{i=1}^{K(t)}L_i$, with $E[N]=\Lambda_\mu E[L]$ and
  $\operatorname{Var}[N]=\Lambda_\mu E[L^2]$ (so the dispersion index is $E[L^2]/E[L]$).

Because data are daily, at most one event per day is observed; when the daily rate is small, the day-by-day
Bernoulli process is well approximated by a Poisson process (law of rare events).

### Parameter estimation

The log-likelihood $\ell(\lambda)=n\log\lambda-\lambda T$ gives the maximum-likelihood estimate (MLE)
$\hat\lambda=n/T$. An exact 95 % confidence interval follows from the chi-square relation (Garwood):
$[\chi^2_{0.025,2n}/2T,\ \chi^2_{0.975,2n+2}/2T]$. For the NHPP the same is applied month by month.

**HPP.** n = 300, T = 3430 → **λ̂ = 0.0875 per day (95 % CI 0.078–0.098) ≈ 31.9 Severe days per year
(28.4–35.8)**. Mean waiting time 1/λ̂ = 11.4 days; P(at least one Severe day in a week) = 0.46.

**Monthly NHPP (Figure A08).**

| Month | Jan | Feb | Mar–Apr | May | Jun–Sep | Oct | Nov | Dec |
|---|---|---|---|---|---|---|---|---|
| Severe days | 85 | 5 | 0 | 3 | 0 | 11 | 102 | 94 |
| λ̂ₘ per day | 0.296 | 0.018 | 0 | 0.010 | 0 | 0.038 | **0.345** | 0.312 |

94 % of Severe days fall in November–January, where about one day in three is Severe. Apart from three
isolated days in May, there were none from March to September in ten years. By season: Winter (Dec–Feb)
0.214/day, Post-monsoon (Oct–Nov) 0.194/day, Summer 0.004/day, Monsoon 0.

### Model checking

**HPP is rejected.**

* *Dispersion.* Monthly counts have mean 2.5 and variance 21.4; the Pearson dispersion index is **8.1**
  (should be ≈ 1; χ² = 965, df = 119, p < 10⁻¹⁰⁰) (Figure A06).
* *Inter-arrival times* (Figure A07): the mean gap is 5.6 days, not 11.4. **66 % of gaps are exactly one day**,
  against 9 % under the HPP, and a few gaps last 100–280 days (the summer). KS test vs Exp(λ̂): D = 0.58, p ≈ 10⁻⁹³.

**Seasonality is the main effect.** Likelihood-ratio tests: HPP vs seasonal NHPP LR = 490 (df 3),
HPP vs monthly NHPP LR = 680 (df 11), seasonal vs monthly LR = 190 (df 8); all p < 10⁻³⁰.
AIC: HPP 928 → seasonal 443 → **monthly 270**. The monthly NHPP follows the observed staircase of
cumulative counts (Figure A09), and the dispersion falls from 8.1 to **1.8** (still > 1, p < 0.001).

**No trend over the years.** A Poisson regression with month and year effects and a log-exposure offset
shows that years differ (LR = 24.7, df 9, p = 0.003). After allowing for overdispersion (quasi-Poisson,
φ = 1.54) this is not significant (p = 0.10). The linear trend is **−0.4 % per year (95 % CI −5.6 % to +5.0 %)**.
Observed/expected ratios (Figure A10) range from 0.51 (2022) to 1.48 (2016) and 1.41 (2023).

**Clustering.** 63 % of Severe days immediately follow another Severe day.

* *Time-rescaling theorem.* If the model is correct, $\tau_k=\Lambda(S_k)-\Lambda(S_{k-1})$ are i.i.d. Exp(1)
  (Figure A12). KS distance: HPP 0.51, NHPP for Severe days 0.23 (p ≈ 10⁻¹⁴), **NHPP for episode onsets 0.17**
  (p = 0.004).
* *Episode onsets* behave like a Poisson process: dispersion index 0.89 (p = 0.71).
* *Episode lengths* (Figure A11) are roughly geometric with q̂ = 0.63 (P(Severe tomorrow | Severe today)),
  but with a heavier tail (up to 17 days; χ² GOF p = 0.02).

The remaining misfit of onsets is expected, because a new episode cannot start inside a running one.
Geometric run lengths are exactly what a two-state Markov chain produces. This is the link to the
Markov-chain model.

### Simulation

Each fitted model was simulated 2000 times on the same 3653-day calendar (Figures A13–A14).
All three reproduce the 10-year total (≈ 300; 95 % range 266–335 for HPP and NHPP, 230–375 for compound),
but they differ in the dispersion of monthly counts (observed 1.77):

| Model | simulated dispersion (mean, 95 % range) | P(sim ≥ observed) |
|---|---|---|
| HPP | 10.7 (7.0–15.8) | 1.00 |
| Monthly NHPP | 1.1 (0.7–1.7) | 0.02 |
| Compound NHPP | 8.1 (3.0–18.2) | 1.00 |

The NHPP is slightly too regular and the compound model is far too variable. In reality clustering is
self-limiting: a month has only ~30 days, and episodes cannot overlap.

### Out-of-sample validation

All models were re-fitted on **2016–2022** and used to predict 2023, 2024, 2025 and Jan–Sep 2026
(Figure A15).

| Year | Observed | HPP 95 % PI | Monthly NHPP 95 % PI | Compound 95 % PI |
|---|---|---|---|---|
| 2023 | 44 | 18–39 ✗ | 19–40 ✗ | 9–51 ✓ |
| 2024 | 34 | 19–39 ✓ | 18–40 ✓ | 10–53 ✓ |
| 2025 | 35 | 19–39 ✓ | 18–39 ✓ | 9–54 ✓ |
| 2026 (Jan–Sep) | 7 | 13–31 ✗ | 5–17 ✓ | 0–25 ✓ |

| Monthly predictions (45 test months) | HPP | **Monthly NHPP** | Compound |
|---|---|---|---|
| 95 % interval coverage | 0.80 | **0.956** | 0.978 |
| Mean log predictive score | −3.36 | **−0.99** | −1.05 |
| Mean absolute error (days) | 3.76 | **1.09** | 1.32 |

The monthly NHPP is the best predictive model. Its coverage is almost exactly the nominal 95 %, and it has the
best log score and the lowest error. Its only miss is the unusually bad winter of 2023.

### Sensitivity to the threshold

With c = 120 µg/m³ there are 1186 events (126 per year), mean episode length 7.4 days, and dispersion 2.4
after the seasonal fit (vs 1.8 for c = 250). The Poisson description works best for the rarer Severe threshold.

### Interpretation

* Severe days at R K Puram occur at about **32 per year**, almost all in **November–January**, when roughly
  one day in three is Severe. That is a pattern of winter meteorology (low mixing height, calm winds)
  combined with post-harvest crop burning, not of random year-round events.
* There is **no statistically detectable improvement** between 2016 and 2025.
* A **seasonal Poisson process** is a good model for *how many* Severe days a month or winter will have.
  It is not a good model for *which* days, because Severe days persist (P(Severe tomorrow | Severe today) ≈ 0.63).
  The Markov-chain model captures that persistence.

### Limitations (Poisson part)

* One station only; other parts of Delhi may differ.
* Daily resolution: at most one event per day, so the Poisson model is an approximation.
* 6 % of days missing and 5 % interpolated; interpolation can create or remove a few threshold crossings.
* Rates assumed constant within each calendar month; no weather covariates.
* Year-to-year variation (e.g. 2023) is larger than a Poisson process allows.

**Possible extensions:** Poisson regression with weather covariates; Cox (doubly stochastic) process for
year-to-year variation; Hawkes (self-exciting) process for clustering; multi-station spatial model.

---

## Slide notes (≈ 4 minutes)

1. **Problem + data (1 min)**: why Severe days matter; CPCB R K Puram hourly data → daily, 2016–2025
   (A01 coverage, A02 time series). Event = PM2.5 > 250.
2. **Poisson model (1 min)**: N(t), HPP formulas, λ̂ = 0.088/day ≈ 32/yr with CI. Then show that the HPP
   fails: dispersion 8.1 and 66 % one-day gaps (A06/A07).
3. **Seasonal NHPP (1 min)**: monthly rates (A08), cumulative staircase (A09), LRT/AIC; no trend (A10).
4. **Clustering + validation (1 min)**: episodes and time-rescaling (A11/A12), simulations (A13/A14),
   out-of-sample table (A15). Takeaway: NHPP for counts, Markov chain for persistence → hand over to Member B.

**Individual contribution (Member A):** downloaded and cleaned the CPCB data and produced the shared
daily dataset; did the exploratory analysis and chose the threshold; built the HPP, NHPP and compound
Poisson models (estimation, likelihood-ratio tests, time-rescaling checks, simulation and out-of-sample validation).

## References (Poisson part)

* Central Pollution Control Board, *National Air Quality Index* (2014) — AQI breakpoints.
* CPCB CCR, *CAAQMS Data Repository*, https://airquality.cpcb.gov.in/ccr/ (accessed Oct 2026).
* Ross, S. M., *Introduction to Probability Models*, Academic Press — Poisson and non-homogeneous Poisson processes.
* Brown, E. N., Barbieri, R., Ventura, V., Kass, R. E., Frank, L. M. (2002), "The time-rescaling theorem and
  its application to neural spike train data analysis", *Neural Computation* 14(2), 325–346.
* Garwood, F. (1936), "Fiducial limits for the Poisson distribution", *Biometrika* 28, 437–442.
* Achcar, J. A., Fernández-Bremauntz, A. A., Rodrigues, E. R., Tzintzun, G. (2008), "Estimating the number
  of ozone peaks in Mexico City using a non-homogeneous Poisson model", *Environmetrics* 19(5), 469–485.
