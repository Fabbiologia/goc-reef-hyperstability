# Hyperstable catch records hid two decades of reef decline in the Gulf of California

Analysis companion repository. It contains the complete, runnable analysis
behind the manuscript by Favoretto, Carmona-Ruiz, López-Sagástegui, Haddock,
Fraser, Tamayo, Sánchez-Ortíz, Guidetti, León-Solorzano and Aburto-Oropeza
(in review): the pipeline that turns the raw survey and landings data into
every intermediate table, every figure and every number quoted in the text.
The manuscript itself is not distributed here.

Two companion papers use the same monitoring programme with non-overlapping
analyses: Carmona-Ruiz et al. (in review at Marine Ecology Progress Series),
on the erosion of size structure, and a manuscript in preparation on where
reef fish production originates along the Gulf's latitudinal gradient.

## What the analysis shows

Twenty-seven years of underwater surveys on Gulf of California rocky reefs,
read against the complete official landing record of the small-scale fishery.
The reef community reorganised along temperature lines through four basin-scale
marine heatwaves and the average monitored reef lost about a quarter of its
fish, while reef landings more than doubled and catch per unit effort carried
no usable signal about the stock (the proportionality assumption behind
catch-based management is rejected at p = 1.3e-6). The pipeline also measures
the engine of that illusion: biomass turnover rising as stocks fall, carried
by a fish community whose production runs almost entirely on the reef's own
energy, with a pelagic subsidy that is small, declining fastest and failing
in warm years.

## Layout

```
data/                       inputs
  conapesca_species_lookup.csv   species code -> scientific name
  ltem_fish_traits.csv           growth-model trait table (186 taxa)
  ltem_name_lookup.csv           species-name harmonisation
  conapesca_ingest_log.csv       provenance log of the raw ingest
  env/                           SST series derived from NOAA OISST v2.1
  (ltem.parquet and CONAPESCA Raw/ must be obtained; see Data access)
analysis/
  scripts/                  the numbered pipeline (below)
  data/                     derived intermediates (committed, inspectable)
  figures/                  all main and supplementary figures, plus
                            in_text_statistics.csv with every number
                            quoted in the manuscript
convert_ltem.R              builds data/ltem.parquet from the LTEM database
```

## The pipeline

Run from `analysis/scripts/`:

```
Rscript 00_run_all.R
```

which executes, each as an isolated process:

| Step | What it does |
|---|---|
| `01a_ingest_conapesca_raw.py` | tidy Pacific wild-capture extract from the raw government landing receipts (skipped if cached) |
| `01_data_preparation.R` | artisanal landings, reef-species classification, SST/MHW copy |
| `02_climate_analysis.R` | SST trend and warm-season anomaly |
| `03_ltem_analysis.R` | reef-fixed-effect biomass trends, Figure 2 trajectories |
| `03b_productivity.R` | biomass production and turnover on the balanced 26-reef panel |
| `03c_pathways.R` | energy pathways and the subsidy interaction |
| `04_fishery_analysis.R` | CPUE models and lag scans |
| `05_variation_partitioning.R` | climate-versus-fishery variance partitioning |
| `06_economic_value.R` | ex-vessel value, five states and La Paz/Loreto |
| `07_decoupling.R` | survey versus landings: the hyperstability test |
| `08_gap_analysis.R` | climate, fishing and their combined effect |
| `09_figures.R` | all figures and `in_text_statistics.csv` |

A fresh checkout with the raw inputs in `data/` rebuilds every intermediate
CSV in `analysis/data/` and every figure in `analysis/figures/`. Because the
derived intermediates are committed, steps `09` (figures) and most statistical
summaries can be re-run and inspected without any raw data.

## Data access

- **NOAA OISST v2.1** (sea surface temperature): freely available from NCEI.
  The derived Gulf series used by the pipeline are committed in `data/env/`.
- **CONAPESCA landing receipts** ("avisos de arribo", ~5.4 GB, 2000 to 2026):
  public records available on request from CONAPESCA. Place the yearly files
  in `data/CONAPESCA Raw/` and `01a` rebuilds the tidy extract.
- **LTEM visual census**: available on request from dataMares
  (www.datamares.org). Build `data/ltem.parquet` with `convert_ltem.R`.

## Requirements

R (>= 4.3) with `data.table`, `arrow`, `ggplot2`, `patchwork`, `lubridate`,
and `rfishprod` (github.com/renatoamorais/rfishprod); Python 3 with
`openpyxl` for the raw ingest step.

## License

Code is released under the MIT License (see `LICENSE`). Data files derived
from third-party sources (NOAA, CONAPESCA, LTEM) remain subject to their
providers' terms.
