# Stochastic modelling of Severe air-pollution days in Delhi

ASM project (2026M): data-driven stochastic modelling of air-pollution exceedance events, using CPCB hourly
PM2.5 data from **R K Puram, Delhi (DPCC, site 124)**, 2016–2026. A Severe day is a day with mean PM2.5 above
250 µg/m³.

| Part | Member | Notebook |
|---|---|---|
| Data cleaning, exploratory analysis, Poisson process | Harsh Raj | `01_data_cleaning_eda.ipynb`, `02_poisson_model.ipynb` |
| Markov chain on AQI categories | Vivian Mathew | `03_markov_chain.ipynb` |
| Simulation engine and validation | Rudraksh Kapil | `04_simulation_validation.ipynb` |

## Folder layout

```
data/raw/        CPCB hourly files, one per year (2009–2026), unchanged as downloaded
data/clean/      delhi_rkpuram_daily.csv: the cleaned daily dataset used by all notebooks
notebooks/       analysis notebooks 01–04 (run from inside this folder, in order)
figures/         every plot saved by the notebooks (A = Poisson, B = Markov, C = simulation)
results/         JSON outputs: summaries, fitted parameters, simulation and validation results
report/          report source (report.md; report_full.md is the uncut version) and the report PDF
presentation/    presentation slides (PDF)
```

## How to run

```bash
pip install -r requirements.txt
cd notebooks
jupyter notebook
```

Run `01`, `02`, `03` and `04` in order; each reads the outputs of the previous ones from `results/`. Each
notebook runs top to bottom in under a minute, with a fixed random seed.

## Cleaned dataset: `data/clean/delhi_rkpuram_daily.csv`

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

## Key results

**Poisson process (notebook 02).** 300 Severe days in 3,430 valid days: λ̂ ≈ 32 per year (95% CI 28–36). The
homogeneous process is rejected (monthly dispersion 8.1; 66% of gaps are one day vs 9% expected). Monthly rates
fix seasonality: Nov–Jan hold 94% of Severe days (LR = 680, df = 11), with no significant trend. Fitted on
2016–2022, the monthly NHPP covers 95.6% of 2023–2026 months.

**Markov chain (notebook 03).** Five states, 3,394 transitions. P(Severe → Severe) = 0.631 (95% CI 0.58–0.68),
mean Severe spell 2.7 days; one winter day in five is Severe in the long run. A first-order seasonal chain is
adequate within winter (order test p = 0.69). The first Severe day comes about 29 days after 1 October (observed
mean 29.6). The chain gives the best next-day forecast (Brier 0.051 vs 0.060).

**Simulation and validation (notebook 04).** One engine simulates all models 2,000 times on the real calendar;
parameter recovery confirms unbiased estimates with standard error ∝ N^(−1/2). All models reproduce the
ten-year total, but none reproduces the observed monthly dispersion of 1.77 (NHPP 1.1, Markov 6.8). The Markov
chain reproduces episode lengths (KS p = 0.32) and week-long winter crises (33% of winters vs 4 of 11 observed;
NHPP ≈ 1%). Out of sample it wins one day ahead; the NHPP has the lower Brier score from two days on. Policy
scenarios (with assumed effect sizes) show that reducing persistence removes far more 5-day crises (−45%) than
reducing entry into Severe (−17%).

Known issue: the compound-NHPP simulator in notebook 04 places all of an episode's days on its onset day, so
that model's episode-length comparison is not valid; its totals and dispersion are unaffected.

Full write-up: `report/report.md` (PDF: `report/ASM_Project_Report.pdf`).

## Data source

Central Pollution Control Board, Central Control Room for Air Quality Management, CAAQMS Data Repository,
https://airquality.cpcb.gov.in/ccr/ (downloaded October 2026).
