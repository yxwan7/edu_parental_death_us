Quantify the burden of children experiencing parental death by parental education in the United States, 2010-2024
=================================================================================================================

## Repository layout

```
01_data/        inputs and derived panels; see "Data" below
02_codes/       run in numeric order, from the repository root
├── 00_functions.r                 shared definitions and helpers
├── 01_download_nchs_files.r
├── 02_process_births.r
├── 03_process_deaths.r
├── 04_process_population.r
├── 05_mortality_rates.r
├── 06_mortality_rates_by_cause.r
├── 07_fertility_rates.r
├── 08_male_fertility_rates.r
├── 09_education_transitions.r
├── 10_kinship_model.r
├── 11_parental_death.r
├── 12_point_estimates.r
├── 13_tables_figures_point_estimates.r
├── 14_tables_figures_publication.r
├── uncertainty/                   replicate weights, draws, aggregation, Slurm scripts
└── sensitivity/                   Supplementary Material analyses
03_outputs/     the reported results
├── 00_point_results/              Table 1 and Figures 1-4, point estimates only
├── 01_tables/                     Table 1, as CSV, numeric CSV and LaTeX
├── 02_figures/                    Figures 1-4
├── 03_appendix/                   Tables S4-S7, Figures S1-S10, sensitivity CSVs
└── 04_figure_data/                the table behind each figure, as CSV
```

## Data

The inputs are public but not redistributed here. Place each at the path shown
and the scripts will find it. `01_data/02_processed/` keeps only the aggregated
panels and estimates read from `10_kinship_model.r` onward, including the
interval bundle of the reported batch (`draws/ui_20260827/final`).

```
01_data/01_raw/
├── nchs/natality/    one .zip per year, 1990-2024
├── nchs/mortality/   one .zip per year, 1990-2024
├── acs/              one IPUMS-USA extract: a .xml companion and its .dat file
└── external/nsch_validation_race.csv
```

NVSS natality and mortality come from
[NCHS](https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm#Mortality_Multiple),
and the ACS from [IPUMS-USA](https://usa.ipums.org/usa/).

## Running the analysis

Run everything from the repository root.

### Point estimates

```bash
Rscript 02_codes/01_download_nchs_files.r … 02_codes/13_tables_figures_point_estimates.r
```

`01` builds a file inventory but downloads nothing unless
`NCHS_RUN_DOWNLOAD=true`. `10_kinship_model.r` takes about fifteen minutes; the
rest are faster. `13` writes Table 1 and Figures 1-4 without intervals to
`03_outputs/00_point_results`.

### Uncertainty

Reported values are draw medians with 2.5th/97.5th percentiles over 1000 Monte
Carlo draws. Each draw re-runs stages `05`-`11` on perturbed inputs — about an
hour and 25 GB each, so a cluster is needed.

```bash
Rscript 02_codes/uncertainty/01_replicate_weights.r  # once

export EDU_PARENT_LOSS_SSH_HOST=<login host>  # Slurm cluster
02_codes/uncertainty/hpc/run_batch.sh smoke ci_YYYYMMDD 1000
02_codes/uncertainty/hpc/run_batch.sh draws ci_YYYYMMDD 1 1000 40 1000
02_codes/uncertainty/hpc/run_batch.sh status ci_YYYYMMDD  # wait for the array
02_codes/uncertainty/hpc/run_batch.sh fetch ci_YYYYMMDD 01_data/02_processed/draws/ci_YYYYMMDD

export EDU_PARENT_LOSS_SHARD_ROOT=01_data/02_processed/draws/ci_YYYYMMDD
export EDU_PARENT_LOSS_CI_ROOT=$EDU_PARENT_LOSS_SHARD_ROOT/final
Rscript 02_codes/uncertainty/03_aggregate_draws.r
Rscript 02_codes/14_tables_figures_publication.r
```

Draw and smoke jobs request one core, 32G and twelve hours; the prepare job that
precedes the array asks for 2G and twenty minutes. All of them expect a conda
environment named `edu_parental_death` (override with
`EDU_PARENT_LOSS_CONDA_ENV`). A single draw also runs locally through
`uncertainty/02_run_draw.r`. `03_outputs` holds the results of the final
1,000-draw batch.

### Sensitivity

Pass a leg name to both `10_kinship_model.r` and `11_parental_death.r`, then run
`sensitivity/01_male_fertility.r`, which writes the two annual RD/RR tables and
regenerates Figures S9 and S10:

```
s1_primary  s1_female_borrow  s1_zero_gradient
s2_complete_case  s2_mnar_all_ba  s2_mnar_all_noba
```

## Requirements

R 4.4 or later with `data.table`, `here`, `DemoKin`, `ipumsr`, `tidyselect`, and
`rvest` (download step only). Figures also need `ggplot2`, `scales`, `patchwork`
and `ragg`; cluster draws do not. The published estimates used R 4.6.1,
`data.table` 1.18.4, `here` 1.0.2 and `DemoKin` 1.0.3.
