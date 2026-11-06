# Member C (Rudraksh): Simulation, model comparison & final report

Parts A (Poisson process) and B (Markov chain) are finished, and each already has its own out-of-sample
validation. **Don't refit or re-validate those models.** Your part is to simulate both, compare them
against each other and the data, and put everything together into the final report and talk.

## Inputs

| File | What it is |
|---|---|
| `data/clean/delhi_rkpuram_daily.csv` | shared daily dataset (see README for the columns) |
| `results/poisson_params.json` | Part A: HPP rate, monthly NHPP rates, episode onset rates, q = P(stay Severe) |
| `results/markov_params.json` | Part B: 5-state transition matrices (all-year, seasonal, monthly, train-period) |
| `results/A02_poisson_results.json`, `results/B01_markov_results.json` | all numbers from A and B |
| `report/memberA_sections.md`, `report/memberB_sections.md` | A's and B's report drafts |

Event definition (fixed): a **Severe** day is one with daily PM2.5 > 250 µg/m³. Study period 2016–2025; train
2016–2022, test 2023 → Sep 2026.

## Tasks

### 1. Markov-chain simulator (theory vs simulation)
- Simulate at least 1000 daily category paths over the 2016–2025 calendar, using the seasonal (or monthly) matrices
  and switching matrix when the season changes.
- Check that the simulations reproduce Part B's theory:
  - probability of each category several days ahead (Pⁿ)
  - long-run share of days in each category (stationary distribution π)
  - mean spell length in each category, 1/(1 − pᵢᵢ)
  - expected time until a Severe day
- Show it as a theory vs simulation table plus plots. A close match verifies both the formulas and the code.

### 2. Head-to-head: which model reproduces the data?
Simulate, on the same calendar, Part A's monthly NHPP, Part A's compound (clustered) NHPP, and the seasonal
Markov chain. For each model, compare the simulated distribution with the observed value of:

| Feature | Observed value |
|---|---|
| Severe days per year | see `A02_poisson_results.json` → `yearly_OE` |
| Dispersion index of monthly counts | 1.77 |
| Longest Severe episode in each winter | maximum 17 days |
| Date of the first Severe day each season | compute from the data |

Report where the observed value falls in each simulated distribution (a percentile or p-value), and rank which
model gets which feature right.

### 3. Parameter uncertainty (bootstrap)
- Use a block bootstrap that resamples whole winters (Oct–Feb seasons) to get intervals for the Poisson rate λ and for
  P(Severe → Severe).
- Compare them with the formula-based intervals from A (Garwood) and B. Explain any difference: resampling
  whole winters captures the year-to-year variation, such as the bad winter of 2023.

### 4. One combined results table
- Put A's and B's out-of-sample results (fit 2016–2022, test 2023 → Sep 2026) into a single comparison table.
- Write the Discussion. The expected message: the Poisson NHPP answers *how many* Severe days a month or season
  will have; the Markov chain answers *which days* and *how long* a bad spell lasts. Both estimate the same
  persistence (≈ 0.63).

### 5. Final report and slides (you own these)
- Write: **Abstract, Introduction, Literature review** (4–6 papers on stochastic models of air-pollution
  exceedances), **Methodology overview, Discussion, Conclusion**.
- Merge `memberA_sections.md` and `memberB_sections.md` into one **10–12 page report**, following the
  structure in the project brief, with the code in an appendix.
- Build the **15-minute deck**: each member presents their own part and their own 1-minute individual
  contribution slide; leave 4 minutes for Q&A.

## Deliverables
- `notebooks/04_simulation_validation.ipynb`
- figures `figures/C01_…`, `C02_…`, …
- `results/C01_*.json`
- `report/memberC_sections.md`
- final report (PDF) and slides

## Conventions (same as notebooks 01–03)
- Fixed random seed; the notebook runs from top to bottom (`Restart & Run All`) without errors.
- A short markdown cell explaining each step and a "Reading" cell interpreting each result.
- Every plot has a title, axis labels with units, and a legend; save each one with `savefig` to `figures/`.
- Run 01 → 02 → 03 first if any `results/` file is missing.
