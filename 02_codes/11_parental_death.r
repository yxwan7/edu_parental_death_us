# ==============================================================================
# 11_parental_death.r
#
# Combine cause-specific parental mortality with living-parent exposure
# following Schlüter et al. (2024).
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/11_parental_death.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Parent-loss definitions ----

# Flag-based 3-way split (rate/02b) and the ICD-10 cause_15/cause_6 panel (rate/02).
CAUSE_LEVELS <- c("drugs", "firearms", "other")
CAUSE_15_LEVELS <- c(
  "drug_poisoning_overdose",
  "firearm_homicide",
  "firearm_suicide",
  "other_firearm_deaths",
  "covid19",
  "heart_disease",
  "cerebrovascular_disease",
  "diabetes_mellitus",
  "cancer",
  "chronic_lower_respiratory_disease",
  "alcohol_attributable",
  "nonfirearm_suicide",
  "nonfirearm_homicide",
  "other_unintentional_injuries",
  "all_remaining_causes"
)
TOL <- 1e-8

# Male-fertility sensitivity leg, matching the argument given to 10. With no
# argument this is the mainline; a leg reads 10's leg output and writes beside
# it under 01_data/02_processed/edu_parent_loss/sensitivity/<leg>/.
leg <- commandArgs(trailingOnly = TRUE)[1]
parent_loss_dir <- if (is.na(leg)) {
  here("01_data", "02_processed", "edu_parent_loss")
} else {
  here("01_data", "02_processed", "edu_parent_loss", "sensitivity", leg)
}
dir.create(parent_loss_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Parent-loss files ----

rate_dir <- here("01_data", "02_processed", "rate")
exposure_input_path <- file.path(
  parent_loss_dir, "02_education_parent_alive_exposure_pooled.rds"
)
# Flag-based drugs/firearms/other cause rate panel (rate/02b) feeds the 3-way
# family; the ICD-10 cause_15/cause_6 panel (rate/02) feeds the
# 6/15-cause family. Both are shipped rate products of 06.
mortality_by_cause_input_path <- file.path(
  rate_dir, "02b_education_mortality_rate_panel_by_cause_flag3.rds"
)
cause_6_15_mortality_input_path <- file.path(
  rate_dir, "02_education_mortality_rate_panel_by_cause.rds"
)
population_input_path <- here(
  "01_data", "02_processed", "population", "01_pop_main_panel.rds"
)

annual_output_rds_path <- file.path(parent_loss_dir, "03b_parent_loss_by_3_cause_annual.rds")
child_age_output_rds_path <- file.path(
  parent_loss_dir, "03c_parent_loss_by_3_cause_child_age.rds"
)
cause_6_15_annual_output_rds_path <- file.path(
  parent_loss_dir, "04b_parent_loss_by_6_15_cause_annual.rds"
)
cause_6_15_child_age_output_rds_path <- file.path(
  parent_loss_dir, "04c_parent_loss_by_6_15_cause_child_age.rds"
)

# ---- Parent-loss calculations ----

parent_type_from_sex <- function(sex) {
  fifelse(sex == "female", "mother", fifelse(sex == "male", "father", NA_character_))
}

prepare_child_population <- function(pop_dt) {
  pop_use <- copy(as.data.table(pop_dt))[
    race_eth %chin% RACES_R3 &
      sex %chin% SEX_LEVELS &
      age_years %in% AGES_MINOR &
      year %in% YEARS_MODEL,
    .(
      source = tolower(as.character(source)),
      year = as.integer(year),
      race_eth,
      sex,
      age_child = as.integer(age_years),
      n_persons_weighted = as.numeric(n_persons_weighted)
    )
  ]

  pop_use[, expected_source := "acs"]
  pop_use <- pop_use[source == expected_source]

  child_pop <- pop_use[
    ,
    .(n_children = sum(n_persons_weighted, na.rm = TRUE)),
    by = .(year, race_eth, age_child)
  ]

  expected_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    age_child = AGES_MINOR,
    unique = TRUE
  )
  child_pop <- child_pop[expected_grid, on = .(year, race_eth, age_child)]

  child_pop[order(year, race_eth, age_child)]
}

prepare_parent_exposure <- function(exposure_dt, child_population) {
  exposure_use <- copy(as.data.table(exposure_dt))[
    child_is_minor == TRUE &
      race_eth %chin% RACES_R3 &
      sex_parent %chin% SEX_LEVELS &
      parent_educ2 %chin% EDUC_LEVELS &
      age_child %in% AGES_MINOR &
      year %in% YEARS_MODEL,
    .(
      year = as.integer(year),
      race_eth,
      age_child = as.integer(age_child),
      age_parent = as.integer(age_parent),
      parent_type = parent_type_from_sex(sex_parent),
      parent_educ2,
      n_parents_alive = as.numeric(n_parents_alive)
    )
  ]

  exposure_use <- child_population[
    exposure_use,
    on = .(year, race_eth, age_child)
  ]

  exposure_use[
    ,
    child_parent_dyads := n_children * n_parents_alive
  ]

  dyads <- exposure_use[
    ,
    .(child_parent_dyads = sum(child_parent_dyads, na.rm = TRUE)),
    by = .(year, race_eth, age_child, age_parent, parent_educ2, parent_type)
  ]

  dyads[order(year, race_eth, parent_educ2, parent_type, age_parent, age_child)]
}

aggregate_parent_age_dyads <- function(parent_child_age_dyads) {
  parent_child_age_dyads[
    ,
    .(child_parent_dyads = sum(child_parent_dyads, na.rm = TRUE)),
    by = .(year, race_eth, age_parent, parent_educ2, parent_type)
  ][order(year, race_eth, parent_educ2, parent_type, age_parent)]
}

# Cause-specific parent deaths, for either cause family. The ICD-10 family also
# carries cause_15 / cause_6 through; the 3-way family has only `cause`.
prepare_cause_mortality <- function(mortality_dt, rate_definition_id, cause_col,
                                    cause_levels, carry_cols = character(0)) {
  mort <- copy(as.data.table(mortality_dt))[
    rate_definition == rate_definition_id &
      year %in% YEARS_MODEL &
      sex %chin% SEX_LEVELS &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS &
      get(cause_col) %chin% cause_levels
  ]

  mort <- mort[
    ,
    c(
      list(
        year = as.integer(year),
        race_eth = race_eth,
        age_parent = as.integer(age_years),
        parent_type = parent_type_from_sex(sex),
        parent_educ2 = educ2,
        cause = get(cause_col)
      ),
      lapply(setNames(carry_cols, carry_cols), function(nm) get(nm)),
      list(
        n_parent_deaths = as.numeric(n_deaths),
        n_parent_deaths_all_cause = as.numeric(n_deaths_all_cause),
        cause_share = as.numeric(cause_share),
        cause_share_source = cause_share_source,
        parent_population = as.numeric(n_pop_weighted)
      )
    )
  ]

  mort[order(year, race_eth, parent_educ2, parent_type, age_parent, cause)]
}

build_parent_age_base <- function(parent_dyads, cause_mortality) {
  parent_pop <- unique(
    cause_mortality[
      ,
      .(
        year, race_eth, age_parent, parent_educ2, parent_type,
        parent_population
      )
    ]
  )

  base <- merge(
    parent_pop,
    parent_dyads,
    by = c("year", "race_eth", "age_parent", "parent_educ2", "parent_type"),
    all.x = TRUE,
    sort = FALSE
  )
  base[is.na(child_parent_dyads), child_parent_dyads := 0]

  base[
    ,
    children_per_parent := fifelse(
      parent_population > 0,
      child_parent_dyads / parent_population,
      0
    )
  ]

  base[]
}

build_parent_child_age_base <- function(parent_child_age_dyads, cause_mortality) {
  parent_pop <- unique(
    cause_mortality[
      ,
      .(
        year, race_eth, age_parent, parent_educ2, parent_type,
        parent_population
      )
    ]
  )

  parent_pop_grid <- parent_pop[
    ,
    .(age_child = AGES_MINOR),
    by = .(year, race_eth, age_parent, parent_educ2, parent_type, parent_population)
  ]

  base <- merge(
    parent_pop_grid,
    parent_child_age_dyads,
    by = c("year", "race_eth", "age_parent", "parent_educ2", "parent_type", "age_child"),
    all.x = TRUE,
    sort = FALSE
  )
  base[is.na(child_parent_dyads), child_parent_dyads := 0]
  setnames(base, "child_parent_dyads", "child_parent_dyads_age_child")

  base[
    ,
    children_per_parent_at_child_age := fifelse(
      parent_population > 0,
      child_parent_dyads_age_child / parent_population,
      0
    )
  ]

  base[]
}

add_total_education_rows <- function(dt) {
  total_dt <- dt[
    parent_educ2 %chin% EDUC_LEVELS,
    .(
      n_parent_deaths = sum(n_parent_deaths, na.rm = TRUE),
      n_parent_deaths_all_cause = sum(n_parent_deaths_all_cause, na.rm = TRUE),
      parent_population = sum(parent_population, na.rm = TRUE),
      child_parent_dyads = sum(child_parent_dyads, na.rm = TRUE),
      n_children_losing_parent = sum(n_children_losing_parent, na.rm = TRUE),
      n_children_losing_parent_all_cause_expected = sum(
        n_children_losing_parent_all_cause_expected,
        na.rm = TRUE
      )
    ),
    by = .(year, race_eth, age_parent, cause, parent_type)
  ]
  total_dt[
    ,
    `:=`(
      parent_educ2 = "Total",
      children_per_parent = fifelse(
        parent_population > 0,
        child_parent_dyads / parent_population,
        0
      ),
      cause_share = fifelse(
        n_parent_deaths_all_cause > 0,
        n_parent_deaths / n_parent_deaths_all_cause,
        NA_real_
      ),
      cause_share_source = "education_total_sum"
    )
  ]

  total_dt[
    ,
    parental_loss_per_1000_age := fifelse(
      child_parent_dyads > 0,
      1000 * n_children_losing_parent / child_parent_dyads,
      NA_real_
    )
  ]

  rbindlist(list(dt, total_dt), use.names = TRUE, fill = TRUE)
}

add_total_education_child_age_rows <- function(dt) {
  total_dt <- dt[
    parent_educ2 %chin% EDUC_LEVELS,
    .(
      n_parent_deaths = sum(n_parent_deaths, na.rm = TRUE),
      n_parent_deaths_all_cause = sum(n_parent_deaths_all_cause, na.rm = TRUE),
      parent_population = sum(parent_population, na.rm = TRUE),
      child_parent_dyads_age_child = sum(child_parent_dyads_age_child, na.rm = TRUE),
      n_children_losing_parent_at_child_age = sum(
        n_children_losing_parent_at_child_age,
        na.rm = TRUE
      ),
      n_children_losing_parent_all_cause_expected_at_child_age = sum(
        n_children_losing_parent_all_cause_expected_at_child_age,
        na.rm = TRUE
      ),
      parent_years_lost_under18 = sum(parent_years_lost_under18, na.rm = TRUE),
      parent_years_lost_under18_all_cause_expected = sum(
        parent_years_lost_under18_all_cause_expected,
        na.rm = TRUE
      )
    ),
    by = .(year, race_eth, age_parent, age_child, cause, parent_type)
  ]
  total_dt[
    ,
    `:=`(
      parent_educ2 = "Total",
      remaining_childhood_years = 18 - age_child,
      children_per_parent_at_child_age = fifelse(
        parent_population > 0,
        child_parent_dyads_age_child / parent_population,
        0
      ),
      cause_share = fifelse(
        n_parent_deaths_all_cause > 0,
        n_parent_deaths / n_parent_deaths_all_cause,
        NA_real_
      ),
      cause_share_source = "education_total_sum"
    )
  ]

  total_dt[
    ,
    parental_loss_per_1000_child_age := fifelse(
      child_parent_dyads_age_child > 0,
      1000 * n_children_losing_parent_at_child_age / child_parent_dyads_age_child,
      NA_real_
    )
  ]

  rbindlist(list(dt, total_dt), use.names = TRUE, fill = TRUE)
}

aggregate_annual <- function(age_dt, child_age_dt) {
  annual <- age_dt[
    ,
    .(
      n_parent_deaths = sum(n_parent_deaths, na.rm = TRUE),
      n_parent_deaths_all_cause = sum(n_parent_deaths_all_cause, na.rm = TRUE),
      parent_population = sum(parent_population, na.rm = TRUE),
      n_children_losing_parent = sum(n_children_losing_parent, na.rm = TRUE),
      n_children_losing_parent_all_cause_expected = sum(
        n_children_losing_parent_all_cause_expected,
        na.rm = TRUE
      ),
      n_children_with_alive_parent_under18 = sum(child_parent_dyads, na.rm = TRUE)
    ),
    by = .(year, race_eth, parent_educ2, cause, parent_type)
  ]

  parent_years <- child_age_dt[
    ,
    .(
      total_parent_years_lost_under18 = sum(parent_years_lost_under18, na.rm = TRUE),
      total_parent_years_lost_under18_all_cause_expected = sum(
        parent_years_lost_under18_all_cause_expected,
        na.rm = TRUE
      )
    ),
    by = .(year, race_eth, parent_educ2, cause, parent_type)
  ]

  annual <- merge(
    annual,
    parent_years,
    by = c("year", "race_eth", "parent_educ2", "cause", "parent_type"),
    all.x = TRUE,
    sort = FALSE
  )

  annual[
    ,
    `:=`(
      parental_loss_per_1000 = fifelse(
        n_children_with_alive_parent_under18 > 0,
        1000 * n_children_losing_parent / n_children_with_alive_parent_under18,
        NA_real_
      ),
      all_cause_parental_loss_per_1000_expected = fifelse(
        n_children_with_alive_parent_under18 > 0,
        1000 * n_children_losing_parent_all_cause_expected /
          n_children_with_alive_parent_under18,
        NA_real_
      ),
      mean_parent_years_lost_under18 = fifelse(
        n_children_losing_parent > TOL,
        total_parent_years_lost_under18 / n_children_losing_parent,
        NA_real_
      ),
      parent_years_lost_per_1000_alive_parent_exposure = fifelse(
        n_children_with_alive_parent_under18 > 0,
        1000 * total_parent_years_lost_under18 /
          n_children_with_alive_parent_under18,
        NA_real_
      )
    )
  ]

  setcolorder(
    annual,
    c(
      "year", "race_eth", "parent_educ2", "cause", "parent_type",
      "n_children_losing_parent", "n_children_with_alive_parent_under18",
      "parental_loss_per_1000", "total_parent_years_lost_under18",
      "mean_parent_years_lost_under18",
      "parent_years_lost_per_1000_alive_parent_exposure",
      "n_parent_deaths", "n_parent_deaths_all_cause",
      "n_children_losing_parent_all_cause_expected",
      "all_cause_parental_loss_per_1000_expected",
      "total_parent_years_lost_under18_all_cause_expected",
      "parent_population"
    )
  )
  setorderv(
    annual,
    c("year", "race_eth", "parent_educ2", "cause", "parent_type")
  )

  annual[]
}

finalize_cause_6_15_output <- function(dt, age_cols = character()) {
  out <- copy(dt)
  out[
    ,
    `:=`(
      cause_15 = cause,
      cause_6 = classify_cause_6_from_cause_15(cause)
    )
  ]
  out[, cause := NULL]

  front_cols <- c("year", "race_eth", "parent_educ2", "cause_15", "cause_6", "parent_type", age_cols)
  setcolorder(out, c(front_cols, setdiff(names(out), front_cols)))
  setorderv(out, front_cols)
  out[]
}

# ---- Build the parent-loss inputs ----

exposure_pooled <- readRDS(exposure_input_path)
cause_mortality_panel <- readRDS(mortality_by_cause_input_path)
cause_6_15_mortality_panel <- readRDS(cause_6_15_mortality_input_path)
population_panel <- readRDS(population_input_path)

child_population <- prepare_child_population(population_panel)
rm(population_panel)

parent_child_age_dyads <- prepare_parent_exposure(exposure_pooled, child_population)
parent_dyads <- aggregate_parent_age_dyads(parent_child_age_dyads)
rm(exposure_pooled, child_population)
gc()

cause_mortality <- prepare_cause_mortality(
  cause_mortality_panel,
  rate_definition_id = "cause_specific",
  cause_col = "cause",
  cause_levels = CAUSE_LEVELS
)
rm(cause_mortality_panel)

parent_age_base <- build_parent_age_base(parent_dyads, cause_mortality)
parent_child_age_base <- build_parent_child_age_base(parent_child_age_dyads, cause_mortality)

age_dt <- cause_mortality[
  parent_age_base,
  on = .(year, race_eth, age_parent, parent_educ2, parent_type)
]

age_dt[
  ,
  `:=`(
    n_children_losing_parent = n_parent_deaths * children_per_parent,
    n_children_losing_parent_all_cause_expected =
      n_parent_deaths_all_cause * children_per_parent
  )
]
age_dt[
  ,
  parental_loss_per_1000_age := fifelse(
    child_parent_dyads > 0,
    1000 * n_children_losing_parent / child_parent_dyads,
    NA_real_
  )
]

age_dt <- add_total_education_rows(age_dt)

child_age_dt <- cause_mortality[
  parent_child_age_base,
  on = .(year, race_eth, age_parent, parent_educ2, parent_type),
  allow.cartesian = TRUE
]
rm(parent_child_age_base, cause_mortality)

if ("i.parent_population" %chin% names(child_age_dt)) {
  child_age_dt[, i.parent_population := NULL]
}
child_age_dt[
  ,
  `:=`(
    remaining_childhood_years = 18 - age_child,
    n_children_losing_parent_at_child_age =
      n_parent_deaths * children_per_parent_at_child_age,
    n_children_losing_parent_all_cause_expected_at_child_age =
      n_parent_deaths_all_cause * children_per_parent_at_child_age
  )
]
child_age_dt[
  ,
  `:=`(
    parent_years_lost_under18 =
      n_children_losing_parent_at_child_age * remaining_childhood_years,
    parent_years_lost_under18_all_cause_expected =
      n_children_losing_parent_all_cause_expected_at_child_age *
        remaining_childhood_years,
    parental_loss_per_1000_child_age = fifelse(
      child_parent_dyads_age_child > 0,
      1000 * n_children_losing_parent_at_child_age /
        child_parent_dyads_age_child,
      NA_real_
    )
  )
]

child_age_dt <- add_total_education_child_age_rows(child_age_dt)
setcolorder(
  child_age_dt,
  c(
    "year", "race_eth", "parent_educ2", "cause", "parent_type",
    "age_parent", "age_child", "remaining_childhood_years",
    "n_children_losing_parent_at_child_age",
    "child_parent_dyads_age_child", "children_per_parent_at_child_age",
    "parent_years_lost_under18", "n_parent_deaths",
    "n_parent_deaths_all_cause",
    "n_children_losing_parent_all_cause_expected_at_child_age",
    "parent_years_lost_under18_all_cause_expected", "parent_population",
    "cause_share", "cause_share_source", "parental_loss_per_1000_child_age"
  )
)
setorderv(
  child_age_dt,
  c("year", "race_eth", "parent_educ2", "cause", "parent_type", "age_parent", "age_child")
)

setcolorder(
  age_dt,
  c(
    "year", "race_eth", "parent_educ2", "cause", "parent_type",
    "age_parent", "n_children_losing_parent", "child_parent_dyads",
    "children_per_parent", "n_parent_deaths", "n_parent_deaths_all_cause",
    "n_children_losing_parent_all_cause_expected", "parent_population",
    "cause_share", "cause_share_source", "parental_loss_per_1000_age"
  )
)
setorderv(
  age_dt,
  c("year", "race_eth", "parent_educ2", "cause", "parent_type", "age_parent")
)

annual_dt <- aggregate_annual(age_dt, child_age_dt)

saveRDS(annual_dt, annual_output_rds_path)
saveRDS(child_age_dt, child_age_output_rds_path)

cause_6_15_mortality <- prepare_cause_mortality(
  cause_6_15_mortality_panel,
  rate_definition_id = "cause_6_15",
  cause_col = "cause_15",
  cause_levels = CAUSE_15_LEVELS,
  carry_cols = c("cause_15", "cause_6")
)
rm(cause_6_15_mortality_panel)

cause_6_15_parent_age_base <- build_parent_age_base(parent_dyads, cause_6_15_mortality)
cause_6_15_parent_child_age_base <- build_parent_child_age_base(
  parent_child_age_dyads,
  cause_6_15_mortality
)

cause_6_15_age_dt <- cause_6_15_mortality[
  cause_6_15_parent_age_base,
  on = .(year, race_eth, age_parent, parent_educ2, parent_type)
]

cause_6_15_age_dt[
  ,
  `:=`(
    n_children_losing_parent = n_parent_deaths * children_per_parent,
    n_children_losing_parent_all_cause_expected =
      n_parent_deaths_all_cause * children_per_parent
  )
]
cause_6_15_age_dt[
  ,
  parental_loss_per_1000_age := fifelse(
    child_parent_dyads > 0,
    1000 * n_children_losing_parent / child_parent_dyads,
    NA_real_
  )
]

cause_6_15_age_dt <- add_total_education_rows(cause_6_15_age_dt)

cause_6_15_child_age_dt <- cause_6_15_mortality[
  cause_6_15_parent_child_age_base,
  on = .(year, race_eth, age_parent, parent_educ2, parent_type),
  allow.cartesian = TRUE
]
rm(cause_6_15_parent_child_age_base, cause_6_15_mortality)

if ("i.parent_population" %chin% names(cause_6_15_child_age_dt)) {
  cause_6_15_child_age_dt[, i.parent_population := NULL]
}
cause_6_15_child_age_dt[
  ,
  `:=`(
    remaining_childhood_years = 18 - age_child,
    n_children_losing_parent_at_child_age =
      n_parent_deaths * children_per_parent_at_child_age,
    n_children_losing_parent_all_cause_expected_at_child_age =
      n_parent_deaths_all_cause * children_per_parent_at_child_age
  )
]
cause_6_15_child_age_dt[
  ,
  `:=`(
    parent_years_lost_under18 =
      n_children_losing_parent_at_child_age * remaining_childhood_years,
    parent_years_lost_under18_all_cause_expected =
      n_children_losing_parent_all_cause_expected_at_child_age *
        remaining_childhood_years,
    parental_loss_per_1000_child_age = fifelse(
      child_parent_dyads_age_child > 0,
      1000 * n_children_losing_parent_at_child_age /
        child_parent_dyads_age_child,
      NA_real_
    )
  )
]

cause_6_15_child_age_dt <- add_total_education_child_age_rows(cause_6_15_child_age_dt)

setcolorder(
  cause_6_15_child_age_dt,
  c(
    "year", "race_eth", "parent_educ2", "cause", "parent_type",
    "age_parent", "age_child", "remaining_childhood_years",
    "n_children_losing_parent_at_child_age",
    "child_parent_dyads_age_child", "children_per_parent_at_child_age",
    "parent_years_lost_under18", "n_parent_deaths",
    "n_parent_deaths_all_cause",
    "n_children_losing_parent_all_cause_expected_at_child_age",
    "parent_years_lost_under18_all_cause_expected", "parent_population",
    "cause_share", "cause_share_source", "parental_loss_per_1000_child_age"
  )
)
setorderv(
  cause_6_15_child_age_dt,
  c("year", "race_eth", "parent_educ2", "cause", "parent_type", "age_parent", "age_child")
)

setcolorder(
  cause_6_15_age_dt,
  c(
    "year", "race_eth", "parent_educ2", "cause", "parent_type",
    "age_parent", "n_children_losing_parent", "child_parent_dyads",
    "children_per_parent", "n_parent_deaths", "n_parent_deaths_all_cause",
    "n_children_losing_parent_all_cause_expected", "parent_population",
    "cause_share", "cause_share_source", "parental_loss_per_1000_age"
  )
)
setorderv(
  cause_6_15_age_dt,
  c("year", "race_eth", "parent_educ2", "cause", "parent_type", "age_parent")
)

cause_6_15_annual_dt <- aggregate_annual(cause_6_15_age_dt, cause_6_15_child_age_dt)

cause_6_15_annual_dt <- finalize_cause_6_15_output(cause_6_15_annual_dt)
cause_6_15_child_age_dt <- finalize_cause_6_15_output(
  cause_6_15_child_age_dt,
  age_cols = c("age_parent", "age_child")
)
rm(parent_dyads, parent_child_age_dyads, cause_6_15_parent_age_base)

saveRDS(cause_6_15_annual_dt, cause_6_15_annual_output_rds_path)
saveRDS(cause_6_15_child_age_dt, cause_6_15_child_age_output_rds_path)
