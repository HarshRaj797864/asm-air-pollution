# Stochastic modelling of Severe air-pollution days in Delhi

ASM project (2026M): data-driven stochastic modelling of air-pollution exceedance events, using
CPCB hourly PM2.5 data from **R K Puram, Delhi (DPCC, site 124)**, 2016–2026.

| Part | Owner | Notebook |
|---|---|---|
| Data cleaning, EDA, Poisson process | Member A | `01_data_cleaning_eda.ipynb`, `02_poisson_model.ipynb` |
| Markov chain on AQI categories | Member B | `03_markov_chain.ipynb` |
| Simulation engine + validation | Member C | `04_simulation_validation.ipynb` (to do) |

## Layout

```
data/raw/      CPCB hourly files, one per year, unchanged as downloaded
data/clean/    delhi_rkpuram_daily.csv: the shared daily dataset (use this)
notebooks/     analysis notebooks, run from inside this folder
figures/       all saved plots (A01-A15 = Member A, B01-B07 = Member B)
results/       JSON outputs: summaries, fitted parameters
report/        report section drafts
```

## Setup

```bash
pip install -r requirements.txt
cd notebooks
jupyter notebook
```

Run `01` first; `02` and then `03` (which reads `results/` from `02`). Each notebook runs top to bottom in under a
minute, and the random seed is fixed.

## Shared dataset: `data/clean/delhi_rkpuram_daily.csv`

One row per day from 2016-01-01 to 2026-09-30.

| Column | Meaning |
|---|---|
| `date` | calendar day |
| `pm25` | daily mean PM2.5 (µg/m³); empty when missing |
| `status` | `observed` (≥ 18 valid hours), `interpolated` (gap ≤ 3 days), or `missing` |
| `valid_hours` | number of valid hourly readings that day |
| `aqi_category` | CPCB AQI class from PM2.5: Good ≤ 30 < Satisfactory ≤ 60 < Moderate ≤ 90 < Poor ≤ 120 < Very Poor ≤ 250 < Severe |
| `year`, `month`, `season` | season: Winter (Dec–Feb), Summer (Mar–May), Monsoon (Jun–Sep), Post-monsoon (Oct–Nov) |
| `period` | `study` (2016–2025, used for fitting) or `holdout` (Jan–Sep 2026) |

Notes for the Markov chain:
* Drop or break the chain at `missing` days. Don't count a transition across a missing day.
* Use the same train/test split as the Poisson part (fit on 2016–2022, test on 2023 onward) so the
  models can be compared.
* `results/poisson_params.json` has P(Severe tomorrow | Severe today) ≈ 0.63 from the episode analysis. The
  Severe→Severe entry of your transition matrix should come out close to this.

## Key results so far (Member A)

* 300 Severe days (PM2.5 > 250) in 2016–2025: λ̂ ≈ 32 per year (95 % CI 28–36).
* Strong seasonality: about 0.3 Severe days per day in Nov–Jan, almost none from March to September. The monthly
  NHPP beats the homogeneous Poisson process (LR = 680, df = 11).
* No significant trend over 2016–2025.
* Severe days cluster in episodes (mean 2.7 days); episode starts behave like a Poisson process.
* Out of sample (fit on 2016–2022), the monthly NHPP gives 95.6 % interval coverage on 2023–2026 months.

Full write-up: `report/memberA_sections.md`.

## Key results (Member B, Markov chain)

* 5 states (Good+Satisfactory, Moderate, Poor, Very Poor, Severe), 3394 day-to-day transitions in 2016–2025.
* P(Severe → Severe) = 0.631 (95 % CI 0.58–0.68), matching the Poisson part's 0.627; mean Severe spell 2.7 days.
* Winter and summer matrices differ strongly (p ≈ 10⁻²⁰). One winter day in five is Severe in the long run.
* Within winter, yesterday alone is enough (first-order test p = 0.69); BIC and out-of-sample scores pick a
  seasonal first-order chain. Day-ahead Severe forecasts beat the monthly NHPP (Brier 0.051 vs 0.060).
* `results/markov_params.json` has the all-year, seasonal and monthly transition matrices for the simulation part.

Full write-up: `report/memberB_sections.md`.

## Data source

Central Pollution Control Board, Central Control Room for Air Quality Management, CAAQMS Data Repository,
https://airquality.cpcb.gov.in/ccr/ (downloaded October 2026).
