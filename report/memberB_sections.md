# Member B — report sections (draft)

> Draft text for the section Member B owns: the **Markov-chain model** of daily AQI categories (formulation,
> estimation, multi-step and long-run behaviour, spell lengths, first-passage times, order test, validation).
> All numbers come from `notebooks/03_markov_chain.ipynb`; figures are in `figures/` (B01–B07).
> Rewrite in your own words before submitting.

---

## Markov-Chain Model

### Why a Markov chain

The Poisson part shows that Severe days are not independent: 63 % of them follow another Severe day, and they
come in episodes of 2.7 days on average. A Poisson process can say *how many* Severe days a month will have,
but not *which* days. A Markov chain models exactly this day-to-day persistence. It also uses the full range
of air quality (not only "Severe or not"), so it can describe how the air gets worse step by step.

### Formulation

| Item | Choice |
|---|---|
| Time parameter | discrete, $n = 0, 1, 2, \dots$ (days) |
| State space | $S$ = {Good/Satisfactory, Moderate, Poor, Very Poor, Severe}: the CPCB AQI classes of daily PM2.5, with Good (≤ 30) and Satisfactory (30–60) merged |
| Random variables | $X_n$ = state on day $n$; spell length $D_i$ in state $i$; first-passage time to Severe $T_S$ |
| Missing days | a transition is used only when both days are valid; the chain is broken at missing days |

$\{X_n\}$ is modelled as a time-homogeneous first-order Markov chain,
$P(X_{n+1}=j \mid X_n=i, X_{n-1},\dots) = p_{ij}$. Good and Satisfactory are merged because Good days are
very rare in winter and both mean acceptable air.

Properties used:
* $n$-step transition probabilities: $P^{(n)} = P^n$ (Chapman–Kolmogorov).
* Stationary distribution $\pi = \pi P$: the chain is irreducible and aperiodic, so $p^{(n)}_{ij} \to \pi_j$, the
  long-run share of days in state $j$; mean recurrence time of state $j$ is $1/\pi_j$.
* Spell length in state $i$: $D_i \sim$ Geometric$(1-p_{ii})$, $E[D_i]=1/(1-p_{ii})$.
* Mean first-passage time to Severe: $m = (I-Q)^{-1}\mathbf 1$, with $Q$ = $P$ without the Severe row and column.

### Estimation

With $n_{ij}$ observed $i\to j$ transitions, the likelihood $\prod p_{ij}^{n_{ij}}$ gives the MLE
$\hat p_{ij} = n_{ij}/n_i$ with standard error $\sqrt{\hat p_{ij}(1-\hat p_{ij})/n_i}$ (Anderson & Goodman 1957).
In 2016–2025 there are **3394 usable transitions** (Figure B01):

| from \ to | G/S | Moderate | Poor | Very Poor | Severe | $n_i$ |
|---|---|---|---|---|---|---|
| G/S | **0.863** | 0.126 | 0.009 | 0.001 | 0.001 | 1186 |
| Moderate | 0.228 | **0.517** | 0.199 | 0.054 | 0.002 | 613 |
| Poor | 0.038 | 0.254 | **0.445** | 0.261 | 0.002 | 422 |
| Very Poor | 0.007 | 0.046 | 0.106 | **0.720** | 0.121 | 875 |
| Severe | 0 | 0.003 | 0.023 | 0.342 | **0.631** | 298 |

All standard errors are below 0.03. The matrix is strongly diagonal, and moves are almost always to a
neighbouring class. Severe is entered almost only from Very Poor and left almost only to Very Poor.

**Link to the Poisson part.** $\hat p_{\text{Severe,Severe}} = 188/298 = 0.631$ (95 % Wilson CI 0.575–0.684).
Member A's episode analysis gave $\hat q = 0.627$; a binomial test finds no difference (p = 0.90).
The implied mean Severe spell, $1/(1-0.631) = 2.71$ days, equals A's mean episode length of 2.68 days.

### Seasonal transition matrices

Transitions were split by season (Figure B02). The seasons are clearly different: a likelihood-ratio test of
equal matrices gives G² = 347, df = 54, p ≈ 10⁻⁴⁴ (winter vs summer alone: G² = 136, df = 18, p ≈ 10⁻²⁰).

| | Winter (Dec–Feb) | Summer (Mar–May) |
|---|---|---|
| Very Poor → Very Poor | 0.74 | 0.55 |
| Very Poor → Severe | 0.14 | 0.01 |
| Very Poor → better | 0.12 | 0.44 |
| Poor → Very Poor | 0.42 | 0.18 |
| Severe → Severe | 0.63 (n = 184) | — (3 Severe days) |

In winter the chain drifts towards Very Poor and Severe; in summer it drifts down. Persistence of Severe is the
same in winter (0.625) and post-monsoon (0.658): seasons change how easily the chain **enters** Severe, not how
long it **stays**. Summer and monsoon Severe rows (3 and 0 Severe days) and the winter Good/Satisfactory row
(13 transitions) cannot be estimated reliably.

### Behaviour over several days

From the winter matrix, starting on a Severe day (Figure B03):

| days ahead | 1 | 2 | 3 | 7 | 14 |
|---|---|---|---|---|---|
| P(Severe) | 0.63 | 0.44 | 0.34 | 0.23 | 0.21 |

From a Very Poor winter day, P(Severe) is 0.14 tomorrow and 0.21 within three days. The influence of today's
state is gone after about two weeks; the second eigenvalue modulus (0.70 in winter) gives a memory
half-life of about 2 days. The all-year chain forgets more slowly ($|\lambda_2| = 0.90$, half-life 6 days)
because its clean state is very sticky.

### Long-run share of days

The stationary distribution matches the observed shares (Figure B04):

| | G/S | Moderate | Poor | Very Poor | Severe |
|---|---|---|---|---|---|
| all-year π | 0.350 | 0.181 | 0.124 | 0.258 | **0.087** |
| winter π (observed) | 0.022 (0.019) | 0.092 (0.084) | 0.125 (0.121) | 0.555 (0.563) | **0.206** (0.214) |

A Severe day recurs on average every $1/\pi_5$ = **11.5 days** over the whole year, the same as Member A's
$1/\hat\lambda$ = 11.4 days, and every 4.9 days in winter. The post-monsoon chain does not reach equilibrium
within its season (October is clean, November is Severe), so a stationary distribution is only meaningful for
roughly stable seasons such as winter.

### Spell lengths

| State | $1/(1-\hat p_{ii})$ | observed mean | longest |
|---|---|---|---|
| Good/Satisfactory | 7.3 | 6.9 | 134 |
| Moderate | 2.1 | 2.1 | 13 |
| Poor | 1.8 | 1.8 | 7 |
| Very Poor | 3.6 | 3.5 | 15 |
| Severe | 2.7 | 2.7 | 17 |

Mean spell lengths agree with the geometric prediction for every state (Figure B05). Moderate and Poor are
short "transit" states. The Severe distribution has a slightly heavier tail than geometric (sd 2.7 vs 2.2;
two 17-day episodes), and clean spells mix short winter breaks with months-long monsoon periods.

### Expected time until a Severe day

Mean first-passage times $m_i$ (days):

| From | all-year chain | winter chain | post-monsoon chain |
|---|---|---|---|
| Good/Satisfactory | 53.0 | 15.6 | 23.4 |
| Moderate | 46.5 | 14.4 | 18.6 |
| Poor | 39.6 | 13.2 | 13.5 |
| Very Poor | 27.5 | 10.2 | 9.3 |

The all-year values average over a calendar the homogeneous chain ignores. To get a realistic answer for a
calendar date, the transition matrix was allowed to change by month (rows with < 10 transitions taken from the
season), and Severe was made absorbing. Starting on **1 October**, the first Severe day of the winter is
expected after **29 days** from a Moderate start (21 from Very Poor, 32 from Good/Satisfactory); P(by 15 Nov)
≈ 0.87–0.94. In the 10 observed seasons it came after 19–37 days (mean **29.6**). The model gets the typical date
right but is too spread out (Figure B06): it allows Severe days in early October that never happened, because
the October matrix averages the clean start and the polluted end of the month.

### Is one day of memory enough?

Anderson–Goodman likelihood-ratio tests ($G^2 = 2\sum n_{ijk}\log(\hat p_{ijk}/\hat p_{jk})$):

| Data | order 0 vs 1 | order 1 vs 2 |
|---|---|---|
| pooled 2016–2025 | G² = 4638, df 16, p ≈ 0 | G² = 186, df 72, p ≈ 10⁻¹² |
| Winter | p ≈ 10⁻¹⁰⁰ | G² = 39.8, df 45, **p = 0.69** |
| Post-monsoon | p ≈ 10⁻¹⁴⁸ | G² = 54.5, df 42, **p = 0.09** |
| Summer | p ≈ 10⁻¹⁰⁴ | G² = 82.1, df 54, p = 0.008 |
| Monsoon | p ≈ 10⁻⁷⁹ | G² = 47.6, df 26, p = 0.006 |

| Order | parameters | log-likelihood | AIC | BIC |
|---|---|---|---|---|
| 0 | 4 | −5036.9 | 10081.7 | 10106.2 |
| 1 | 20 | −2740.5 | 5520.9 | **5643.3** |
| 2 | 100 | −2647.6 | **5495.1** | 6107.1 |

Yesterday clearly matters. Pooled over the year, the first-order chain is rejected, but most of this is caused
by mixing seasons: the day before yesterday carries information about the season. Within winter and
post-monsoon, which contain almost all Severe days, first order is adequate. A weak second-order effect
remains in summer and the monsoon. BIC picks order 1. For Severe days, $P(S\mid S,S) = 0.65$ vs
$P(S\mid \text{not }S,S) = 0.61$ (Fisher p = 0.53; in winter 0.62 vs 0.62), so two Severe days in a row do not make a
third more likely.

### Out-of-sample validation

Fitted on **2016–2022** and used to predict each day of 2023 – Sep 2026 from the previous day
(1333 transitions, 120 Severe days):

| Model | mean log score | accuracy | Brier P(Severe) |
|---|---|---|---|
| Order 0, monthly category frequencies | −1.042 | 0.551 | 0.0598 |
| Order 1, all-year | −0.903 | 0.647 | 0.0521 |
| **Order 1, seasonal** | **−0.853** | 0.646 | **0.0508** |
| Order 2, seasonal | −0.898 | 0.630 | 0.0509 |
| Member A monthly NHPP | — | — | 0.0598 |

The seasonal first-order chain is the best model. The second-order chain is worse out of sample (it
overfits), which agrees with the order tests. For Severe days the Markov chain beats the Poisson NHPP
(Brier 0.051 vs 0.060), because it forecasts ≈ 0.6 after a Severe day and almost 0 after a Moderate day, instead
of the same monthly rate every day. These forecasts are roughly calibrated (Figure B07).

### Interpretation

* Delhi's air quality moves **one category at a time** and is very persistent: a Severe day is followed by
  another with probability 0.63, and a Severe spell lasts 2.7 days on average.
* The winter regime is different from the rest of the year: one winter day in five is Severe, and from a Very
  Poor winter day the next Severe day is on average 10 days away.
* The first Severe day of the season typically comes about four weeks after 1 October (around the end of
  October and the start of November).
* **Today's category is enough** to forecast tomorrow within a season. Seasonal matrices matter much more than
  extra memory.
* Poisson and Markov models answer different questions: the NHPP gives good counts per month, the Markov
  chain gives good day-ahead forecasts. Both estimate the same persistence (0.63).

### Limitations (Markov part)

* One station; categories ignore how far PM2.5 is inside a class.
* Seasonal/monthly matrices with few transitions in some rows; matrices are constant within a month.
* Year-to-year differences (e.g. the bad winter of 2023) and weather are not modelled.
* 6 % of days missing (the chain is broken there) and 5 % interpolated.

**Possible extensions:** transition probabilities depending on weather and fire counts (multinomial logistic
regression); hidden Markov model with latent meteorological regimes; continuous-state autoregressive model
for log PM2.5.

---

## Slide notes (≈ 4 minutes)

1. **Why a Markov chain (0.5 min):** Severe days cluster (A's 0.63); states = 5 AQI classes; daily steps;
   chain broken at missing days.
2. **Transition matrix (1 min):** B01 heatmap; diagonal, one step at a time; p_SS = 0.631 matches A's 0.627.
   Winter vs summer (B02), homogeneity test.
3. **Consequences (1.5 min):** n-step decay after a Severe day (B03); stationary shares match data, recurrence
   11.5 days = A's 1/λ (B04); spell lengths geometric (B05); time to Severe and first Severe day after 1 Oct (B06).
4. **Order + validation (1 min):** yesterday enough within seasons, BIC picks order 1; out-of-sample: seasonal
   first-order chain best, beats NHPP for day-ahead Severe forecasts (B07). Hand over to Member C.

**Individual contribution (Member B):** formulated and estimated the Markov chain on AQI categories
(overall, seasonal and monthly matrices); derived n-step probabilities, stationary distribution, spell
lengths and first-passage times; tested time-homogeneity and the order of the chain; validated the chain out of
sample against the Poisson model.

## References (Markov part)

* Anderson, T. W., Goodman, L. A. (1957), "Statistical inference about Markov chains", *Annals of
  Mathematical Statistics* 28(1), 89–110.
* Ross, S. M., *Introduction to Probability Models*, Academic Press — Markov chains, stationary distributions,
  first-passage times.
* Katz, R. W. (1981), "On some criteria for estimating the order of a Markov chain", *Technometrics* 23(3), 243–249.
* Central Pollution Control Board, *National Air Quality Index* (2014) — AQI breakpoints.
* Brier, G. W. (1950), "Verification of forecasts expressed in terms of probability", *Monthly Weather Review* 78(1), 1–3.
