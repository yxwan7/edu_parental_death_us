Quantify the burden of children experiencing parental death by parental education in the United States, 2010-2024
=================================================================================================================

Replication code for a kinship-model estimate of how many children in the United
States lose a parent, by parental education, race and ethnicity, and cause of
death.

## Repository layout

```
01_data/        inputs and derived panels; not distributed, see "Data" below
02_codes/       run in numeric order, from the repository root
├── 00_functions.r                 shared definitions and helpers
├── 01_download_nchs_files.r
├── 02_process_births.r
├── 03_process_deaths.r
├── 04_process_population.r
├── 05_mortality_rates.r
├── 06_mortality_rates_by_cause.r
├── 07_fertility_rates.r
├── 08_male_fertility_proxy.r
├── 09_education_transitions.r
├── 10_kinship_model.r
├── 11_parental_death.r
├── 12_point_estimates.r
├── 13_replicate_weights.r
├── 14_run_draw.r
├── 15_aggregate_draws.r
├── 16_tables_figures.r
├── hpc/                           Slurm submission scripts for the draw array
└── sensitivity/                   Supplementary Material analyses
03_outputs/     the reported results
├── 01_tables/                     Table 1, as CSV, numeric CSV and LaTeX
├── 02_figures/                    Figures 1-4
└── 03_appendix/                   Tables S4-S7, Figures S1-S10, sensitivity CSVs
```

`16_tables_figures.r` writes Table 1, Tables S4-S7, Figures 1-4 and Figures
S1-S7. Figure S8 comes from `sensitivity/02_model_validation.r`, and Figures S9
and S10 from `sensitivity/01_male_fertility.r`. Appendix tables are provided as
CSV and LaTeX; publication tables are not written as RDS files.

## Data

The inputs are public but not redistributed here. Place each at the path shown
and the scripts will find it.

```
01_data/01_raw/
├── nchs/natality/    one .zip per year, 1990-2024
├── nchs/mortality/   one .zip per year, 1990-2024
├── acs/              one IPUMS-USA extract: a .xml companion and its .dat file
├── external/hfc/m_ASFR_USA.txt
└── external/nsch_validation_race.csv
```

NVSS natality and mortality come from
[NCHS](https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm#Mortality_Multiple),
the ACS from [IPUMS-USA](https://usa.ipums.org/usa/), and male age-specific
fertility rates from the
[Human Fertility Collection](https://www.fertilitydata.org/).

Keep the NVSS archives zipped, with a four-digit year in each file name — that
is how the scripts identify them. One per year, no gaps: processing stops if a
year in 1990-2024 is missing or duplicated. The public-use mortality files
suppress geography but carry the decedent's education, which is what this
analysis uses.

The ACS extract covers the 2000-2024 samples and must include `YEAR`, `SEX`,
`AGE`, `RACE`, `HISPAN`, `EDUCD`, `PERWT`, and `REPWTP1`-`REPWTP80` (the
replicate weights, needed only for the uncertainty layer). Save the fixed-width
data file and its DDI companion, exactly one of each: the scripts refuse to
guess column positions.

`nsch_validation_race.csv` holds the published prevalence Figure S8 is checked
against — Schlüter et al. (2024) and the National Survey of Children's Health,
2016-2020. Transcribe it with the header
`source,race_cat,year,median,lower95,upper95`, where `source` is either
`Matrix Projection Model` or `National Survey of Children's Health` and the
three estimate columns are percentages. Only Figure S8 reads it.

## Running the analysis

Run everything from the repository root.

### Point estimates

```bash
Rscript 02_codes/01_download_nchs_files.r … 02_codes/12_point_estimates.r
```

`01` builds a file inventory but downloads nothing unless
`NCHS_RUN_DOWNLOAD=true`. `10_kinship_model.r` takes about fifteen minutes; the
rest are faster.

### Uncertainty

Reported values are draw medians with 2.5th/97.5th percentiles over 1000 Monte
Carlo draws. Each draw re-runs stages `05`-`11` on perturbed inputs — about an
hour and 25 GB each, so a cluster is needed.

```bash
Rscript 02_codes/13_replicate_weights.r                  # once

export EDU_PARENT_LOSS_SSH_HOST=<login host>             # Slurm cluster
02_codes/hpc/run_batch.sh smoke ci_YYYYMMDD 1000
02_codes/hpc/run_batch.sh draws ci_YYYYMMDD 1 1000 40 1000
02_codes/hpc/run_batch.sh status ci_YYYYMMDD             # wait for the array
02_codes/hpc/run_batch.sh fetch ci_YYYYMMDD 01_data/02_processed/draws/ci_YYYYMMDD

export EDU_PARENT_LOSS_SHARD_ROOT=01_data/02_processed/draws/ci_YYYYMMDD
export EDU_PARENT_LOSS_CI_ROOT=$EDU_PARENT_LOSS_SHARD_ROOT/final
Rscript 02_codes/15_aggregate_draws.r
Rscript 02_codes/16_tables_figures.r
```

Draw and smoke jobs request one core, 32G and six hours; the prepare job that
precedes the array asks for 2G and twenty minutes. All of them expect a conda
environment named `edu_parental_death` (override with
`EDU_PARENT_LOSS_CONDA_ENV`). A single draw also runs locally through
`14_run_draw.r`. `03_outputs` holds the results of the final 1,000-draw batch.

### Sensitivity

Pass a leg name to both `10_kinship_model.r` and `11_parental_death.r`, then run
`sensitivity/01_male_fertility.r`, which writes the two annual RD/RR tables and
regenerates Figures S9 and S10:

```
s1_primary  s1_female_borrow  s1_zero_gradient
s2_complete_case  s2_mnar_all_ba  s2_mnar_all_noba
```

`sensitivity/02_model_validation.r` takes no leg. It reads the pooled
parent-survival exposure and ACS population panels left by the mainline run and
writes Figure S8.

## Requirements

R 4.4 or later with `data.table`, `here`, `DemoKin`, `ipumsr`, `tidyselect`, and
`rvest` (download step only). Figures also need `ggplot2`, `scales`, `patchwork`
and `ragg`; cluster draws do not. Versions are not pinned — the published
estimates used R 4.4.1, `data.table` 1.16.2, `here` 1.0.1 and `DemoKin` 1.0.3.
