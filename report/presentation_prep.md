# Presentation prep: 15 minutes (11 talk + 4 Q&A)

Deck: "Severe Air-Pollution Days in Delhi", 16 slides (version 2). Full speaker notes are on each slide,
already cut to these times. Every chart is drawn from `data/clean/` and `results/*.json`, so all
numbers match the notebooks. Rehearse twice with a timer. If you run long, cut the policy slide first.

## Running order

| # | Slide | Speaker | Time | Cumulative | Brief section |
|---|---|---|---|---|---|
| 1 | Cover (winter calendar) | Harsh | 0:15 | 0:15 | Introduction |
| 2 | The problem | Harsh | 0:35 | 0:50 | Introduction |
| 3 | Data | Harsh | 0:40 | 1:30 | Data |
| 4 | Poisson process | Harsh | 0:50 | 2:20 | Model |
| 5 | Seasonal Poisson | Harsh | 0:50 | 3:10 | Model |
| 6 | Markov chain (state diagram) | Vivian | 0:50 | 4:00 | Model |
| 7 | Memory fades in a week | Vivian | 0:45 | 4:45 | Model / theory |
| 8 | Smog season starts ~1 Nov | Vivian | 0:45 | 5:30 | Results |
| 9 | Methodology and code | Rudraksh | 1:00 | 6:30 | Methodology / code |
| 10 | Scoreboard | Rudraksh | 0:50 | 7:20 | Results |
| 11 | Extreme spells | Rudraksh | 0:40 | 8:00 | Results |
| 12 | Forecast test | Rudraksh | 0:50 | 8:50 | Results |
| 13 | Policy what-if | Rudraksh | 0:35 | 9:25 | Interpretation |
| 14 | Conclusion | Vivian | 0:40 | 10:05 | Interpretation |
| 15 | Contributions | all, 15 s each | 0:45 | 10:50 | Individual contribution |
| 16 | Questions | all | 4:00 | | Q&A |

Talk time is about 10:50: Harsh ≈ 3:20, Vivian ≈ 3:15, Rudraksh ≈ 4:10, including each person's 15 s
contribution turn. Have notebooks 02, 03 and 04 open in Jupyter, already run, in case
you're asked to show code. Slide 9 shows the real simulator lines from notebook 04.

## Numbers to know cold

| | |
|---|---|
| Severe | daily mean PM2.5 > 250 µg/m³ |
| Data | CPCB R K Puram, 2016–2025 fit; 3,653 days, 3,430 valid; 300 Severe days in 112 episodes |
| HPP | λ̂ = 0.0875/day = 31.9/yr, exact CI 28.4–35.8 |
| HPP rejected | dispersion 8.1; 66% of gaps are 1 day (HPP expects 9%) |
| Monthly NHPP | Nov–Jan rates 0.30–0.34/day; LR 680 (df 11); AIC 928 → 270; no trend (−0.4%/yr, CI −5.6 to +5.0) |
| NHPP out of sample | 95.6% of 2023–26 months inside 95% PI; missed 2023 (44 vs 19–40) |
| Markov | 5 states, 3,394 transitions; p_SS = 0.631 (CI 0.57–0.68) vs Poisson episodes 0.627 |
| Markov results | mean Severe spell 2.7 d; winter π(Severe) = 0.21; Very Poor → Severe 10 d; 2nd-order test p = 0.69 |
| Simulation | spells: Markov 2.64 vs real 2.68 (KS p = 0.32); dispersion real 1.77, NHPP 1.1, Markov 6.8 |
| Markov diagram | only 3.5% of days jump two or more AQI bands; winter P(Severe in n days \| Severe) = 0.625, 0.44, 0.34 (n = 1–3) → π = 0.21 |
| Extreme spells | Dec–Feb run of 7+ Severe days: real 4 of 11 winters (36%), Markov 33%, monthly Poisson ≈ 1% (20,000 simulated winters each) |
| Forecast | 1 day: Brier 0.051 vs 0.060, AUC 0.93 vs 0.90; Poisson better on Brier from 2 days |
| Policy (assumed effects) | winters with a 5-day crisis: 66% → 59% (entry −20%), 44% (persistence −15%), 38% (both) |

## Likely questions

**Harsh (data and Poisson)**

* *Why threshold 250, not 60?* At 60, 65% of days exceed, which isn't a rare event. At 250, 8.7% of days exceed, and it is the official Severe category.
* *Why is λ̂ = n/T the MLE?* The likelihood is e^{−λT}(λT)^n/n!; the log-likelihood n log λ − λT is maximised at λ = n/T.
* *What is the dispersion index?* Pearson X²/df of the monthly counts. A Poisson process has variance = mean, so it should be ≈ 1. We got 8.1 for the HPP and 1.8 for the NHPP.
* *What is the time-rescaling theorem?* If the intensity is right, Λ(Sₖ) − Λ(Sₖ₋₁) are iid Exp(1). We checked that with a KS test; it showed seasonality is fixed by the NHPP, but clustering is not.
* *Daily data, yet a continuous-time Poisson process?* At most one event per day is observed. With small rates a Bernoulli day process ≈ Poisson (law of rare events). It's a stated approximation; at a winter rate of 0.3/day it is weaker, which is one reason we add the Markov chain.
* *How did you handle missing data?* A day counts only if ≥ 18 of 24 hours are valid. Gaps ≤ 3 days are interpolated; longer gaps are excluded from the exposure T, so they don't bias λ̂.
* *Why did the Poisson model miss 2023?* Year-to-year variation (meteorology, fires) is larger than a Poisson process allows. A Cox process would be the fix.

**Vivian (Markov chain)**

* *Why is p̂ᵢⱼ = nᵢⱼ/nᵢ the MLE?* The likelihood is ∏ pᵢⱼ^{nᵢⱼ}; maximise each row under Σⱼ pᵢⱼ = 1 with a Lagrange multiplier.
* *Why a first-order chain?* The likelihood-ratio test of order 1 vs 2 within winter gives p = 0.69, and BIC also picks order 1. Pooled over the year, order 2 looks significant only because it mixes seasons.
* *Is the chain time-homogeneous?* No across seasons (p ≈ 10⁻²⁰), so we fit one matrix per season. Within a season we assume it is.
* *Why is the mean spell 1/(1 − pᵢᵢ)?* Time in state i is geometric with "leave" probability 1 − pᵢᵢ.
* *How do you get the stationary distribution?* Solve πP = π, Σπ = 1. The chain is irreducible and aperiodic within a season, so π is unique.
* *Do the Poisson and Markov models independently confirm p_SS ≈ 0.63?* No, and don't claim it. Both count the same thing: of the 300 Severe days, 300 − 112 = 188 were followed by another Severe day. The Markov estimate is 188/298 = 0.631 (298 Severe days had an observed next day); the episode estimate is 188/300 = 0.627. It is a consistency check only.
* *How do you get hitting times?* Solve kᵢ = 1 + Σⱼ≠Severe pᵢⱼ kⱼ, with k_Severe = 0.
* *How is the first Severe day after 1 Oct computed?* Make Severe absorbing and iterate vₙ₊₁ = vₙ P_month(n) from that year's 1 Oct state; the absorbed mass is P(first Severe by day n). Real waits were 19–37 days (median 30); the chain's medians are 26–32 days, and it gives 87–93% by 15 Nov (all 10 real seasons were).
* *Where does the season-start model fall short?* Real start dates are more tightly bunched. One October matrix spreads the late-October Diwali and crop-burning spike over the whole month, so the chain gives a 20–50% chance of Severe in early October, which never happened. Weekly or weather-driven matrices would fix it.

**Rudraksh (simulation and validation)**

* *Walk us through the code on slide 9.* Each path's row of today's matrix gives tomorrow's probabilities; cumsum turns them into a CDF; one uniform draw u picks the first state whose CDF exceeds u (inverse-CDF sampling). This is done for all 2,000 paths at once.
* *Why 33% on the spells slide when the 7+ day number was ≈ 30% before?* 33% comes from 20,000 simulated Dec–Feb winters started from the winter stationary distribution; earlier runs with a different set-up gave 28–30%. Either way it is close to the real 36%, and the Poisson model gives ≈ 1%.

* *How do you know the simulator is correct?* Parameter recovery: simulate from known P, re-estimate it, and check the bias goes to 0, the SE shrinks like 1/√N, and the 95% CI covers ≈ 95%. Simulated spell lengths and π also match the formulas.
* *Why doesn't any model match the dispersion 1.77?* The NHPP has no memory (1.1). The Markov chain lets spells pile up within a month (6.8). Real spells are spread more evenly than either model expects. Say this honestly; it is a finding, not a failure.
* *What is the Brier score?* The mean of (forecast probability − outcome)²; lower is better. AUC measures how well the model ranks Severe days above non-Severe days.
* *Why does the Markov forecast get worse with the horizon?* Pʰ → π (Chapman–Kolmogorov), so today's state is forgotten within a few days. The monthly NHPP has finer seasonal detail.
* *Are the policy effects real?* No. The −20% and −15% changes are assumed. The scenario shows only the mechanism: cutting persistence shortens spells, which matters most for long crises.

## Corrections needed in `report/memberC_sections.md` before the final report

Several numbers and claims there don't match the notebook outputs. Fix them from `results/C01_*.json` and
`results/C02_*.json` before merging the report.

1. **Dispersion "paradox resolved"**: the outputs say the opposite. The Markov chain gives 6.8 (monthly) and 20.3
   (seasonal) against the real 1.77, with P(D_sim ≥ D_obs) = 1.00. No model matches; rewrite as on slide 9.
   Also, the variance-inflation formula should use ρ_k = (p_SS − p_{not S→S})^k, not p_SS^k.
2. **Compound Poisson episode test (KS p < 10⁻⁵⁰)**: this is an artifact. The compound simulator puts all of an
   episode's days on its onset day, and `extract_episodes` counts only values equal to 1. Drop this claim, or
   simulate episodes as consecutive days and count values ≥ 1.
3. **Multi-step forecast**: the Markov chain has the lower Brier score only at h = 1. From h = 2 on the NHPP is
   lower (0.0597 vs 0.0608 at h = 2). The Markov chain has the higher AUC up to h = 3. The bold markers in the table are wrong.
4. **Parameter-recovery table**: the values don't match the output (for example N = 200 gives mean 0.614 and
   coverage 93.2%, not 0.622 and 94.2%). Copy from `C01_simulation_results.json`.
5. **Other mismatches**: HPP D = 10.60 (not 10.55); KS p = 0.322 (not 0.312); slide notes say p < 10⁻¹⁰ but the table says
   10⁻⁵⁰; policy baseline 18.4 (not 18.5); 11 winters, not 10.
6. **"17 days in November 2017"**: the 17-day run started on 20 Dec 2018. Another 17-day run started on
   27 Oct 2016, outside the Dec–Feb window.
7. **"(0.32)⁷ ≈ 3.4 × 10⁻⁴, once every three centuries"**: wrong reasoning, because a run can start on any of ~90
   days. Use the simulated probability instead (≈ 0.4% per winter for the NHPP).
8. **"20-year return level" from 10–11 winters**: empirical return levels beyond the sample length can't be
   estimated this way. Either drop it or fit an extreme-value model.
9. **GRAP section**: say plainly that the effect sizes are assumed scenarios, not estimates; "health benefits" overclaims.
