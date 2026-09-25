# Analysis pipeline

The canonical, self-contained pipeline for the manuscript. Run all of it
from this directory:

```
Rscript 00_run_all.R
```

`00_run_all.R` runs each step below as an isolated process (the steps
`setwd("..")` internally, so they must not share an R session). A fresh
checkout with the raw inputs in `../../data` rebuilds every intermediate
CSV in `../data` and every figure in `../figures`.

Required R packages: `data.table`, `arrow`, `ggplot2`, `patchwork`,
`lubridate`, `fixest`, `mgcv`, `heatwaveR`, `rfishprod`, `strucchange`,
`sf`, `rnaturalearth` (with `rnaturalearthdata`). The ingest step needs
Python 3 with `openpyxl`.

| Step | Purpose | Key outputs |
|------|---------|-------------|
| `00a_download_oisst.py` | NOAA OISST v2.1 daily grids for the Gulf, native 0.25°, via NCEI OPeNDAP (cached, ~30 min) | `data/env/oisst_daily/` |
| `00b_build_reef_sst.R` | Basin series, heatwave catalogue, and per-reef SST products from the nearest ocean cell | `sst_gulf_monthly_1981_2026.csv`, `sst_reef_{longterm,year}.csv`, `sst_reef_daily.rds` |
| `01a_ingest_conapesca_raw.py` | Parse the raw government avisos de arribo (three file generations) into one tidy table; back-fill scientific names (cached, ~20 min) | `conapesca_avisos_pacifico_2000_2026.csv.gz`, `conapesca_species_lookup.csv`, ingest log |
| `01_data_preparation.R` | Artisanal (MENORES) filter; rocky-reef classification by genus against the LTEM transects; effort as reef trips (receipts landing ≥1 reef species) with all receipts kept as comparison | `artisanal_bcs_annual.csv`, `artisanal_bcs_yearly_totals.csv`, `artisanal_5state_{annual,effort}.csv` |
| `02_climate_analysis.R` | Whole-Gulf SST trend; warm-season (May–Oct) anomaly; per-latitude trends | `warm_season_anomaly_annual.csv`, `sst_trend_by_latitude.csv` |
| `03_ltem_analysis.R` | Reef-fixed-effect trends and warming betas for invertebrates, fish groups and top species; Figure-2 trajectories | `*_warming*.csv`, `ltem_*_annual.csv` |
| `03b_productivity.R` | Biomass production and turnover on the balanced 26-reef panel (rfishprod Kmax at each reef's climatology) | `productivity_*.csv` |
| `03c_pathways.R` | Five consumer pathways; the pelagic-subsidy share and its warm-year interaction | `pathway_*.csv` |
| `03d_buffer.R` | Buffer over time as % of baseline; regime rates (fishing only, climate only, both) | `buffer_timeseries.csv`, `buffer_regimes.csv` |
| `03e_buffer_climate.R` | Buffer strength Φ; the summer climate interaction with the full inference battery (two-way clustering, wild cluster bootstrap, leave-one-out) | `buffer_reef_year.csv`, `buffer_phi_models.csv` |
| `03f_buffer_nonlinear.R` | GAM alternative; the model-free stratified arbiter by thermal quartile | `buffer_phi_stratified.csv`, verdict |
| `03g_thermal_confounding.R` | Pre-registered checks that thermal exposure is heat, not geography | `thermal_confounding.csv` |
| `03h_buffer_winter.R` | The same buffer battery against winter (Dec–Mar) exposure, the canopy season | `buffer_phi_by_season.csv`, `buffer_winter_*.csv` |
| `03i_growth_scenarios.R` | Kmax temperature scenarios (conservative / tracking / extreme; 30–365 d windows): the primary estimates are a floor | `growth_scenario_*.csv` |
| `03j_buffer_era.R` | Φ before vs within the heatwave era (pre-registered 2014 split) | `buffer_phi_era.csv` |
| `04_fishery_analysis.R` | CPUE (reef landings per reef trip); single-lag scan 0–8 yr; two-mode climate model | `lag_scan_*.csv`, `*_two_mode.csv` |
| `05_variation_partitioning.R` | Borcard partitioning of survey biomass into climate vs fishery components | `varpart_two_mode.csv` |
| `06_economic_value.R` | Ex-vessel value at constant 2022–2025 prices, five states and La Paz/Loreto | `economic_*.csv` |
| `07_decoupling.R` | The hyperstability test (β against the independent survey index); Bai–Perron breaks dating when each fishery series turned | `decoupling_*.csv` |
| `08_gap_analysis.R` | Climate vs fishing vs combined; protection contrasts; Figure S14 | `gap_*.csv` |
| `09_figures.R` | All figures (1–3, S1–S13, S15–S18) and the traceable in-text numbers | `Figure*.pdf/.png`, `in_text_statistics.csv` |

Document builds (`build_word.sh`, `build_gcb.sh`, `docx_linenums_pagenums.py`)
also live here; the archived Science Advances builds are in
`../../archive/sciadv_submission/`.
