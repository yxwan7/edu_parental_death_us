# ==============================================================================
# 10_kinship_model.r
#
# Run the education-stratified DemoKin model and build living-parent exposure
# panels for the main analysis and named sensitivity scenarios.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(DemoKin)
})

here::i_am("02_codes/10_kinship_model.r")

source(here::here("02_codes", "00_functions.r"))

# ---- DemoKin definitions ----

MALE_FERTILITY_RATE_DEFINITION <- "male_proxy"
FOCAL_SEXES <- c("female", "male")
ALLOWED_PARENT_SEX <- c("female", "male")
FEMALE_REPRO_AGES <- 10:54
MALE_REPRO_AGES <- 15:59
BIRTH_FEMALE <- 1 / 2.04
OUTPUT_KIN <- "m"

MALE_PROXY_BY_LEG <- c(
  s1_primary       = "04_education_male_fertility_proxy_panel.rds",
  s1_female_borrow = "04a_education_male_fertility_proxy_panel_s1_female_borrow.rds",
  s1_zero_gradient = "04b_education_male_fertility_proxy_panel_s1_zero_gradient.rds",
  s2_complete_case = "04c_education_male_fertility_proxy_panel_s2_complete_case.rds",
  s2_mnar_all_noba = "04d_education_male_fertility_proxy_panel_s2_mnar_all_noba.rds",
  s2_mnar_all_ba   = "04e_education_male_fertility_proxy_panel_s2_mnar_all_ba.rds"
)

leg <- commandArgs(trailingOnly = TRUE)[1]

# ---- Kinship model files ----

rate_dir <- here("01_data", "02_processed", "rate")
population_dir <- here("01_data", "02_processed", "population")
parent_loss_dir <- if (is.na(leg)) {
  here("01_data", "02_processed", "edu_parent_loss")
} else {
  here("01_data", "02_processed", "edu_parent_loss", "sensitivity", leg)
}
dir.create(parent_loss_dir, recursive = TRUE, showWarnings = FALSE)

mortality_input_path <- file.path(rate_dir, "01_education_mortality_rate_panel.rds")
fertility_input_path <- file.path(rate_dir, "03_education_fertility_rate_panel.rds")
transition_input_path <- file.path(rate_dir, "05_education_transition_panel.rds")
male_fertility_proxy_input_path <- file.path(
  rate_dir,
  if (is.na(leg)) "04_education_male_fertility_proxy_panel.rds" else MALE_PROXY_BY_LEG[[leg]]
)
population_input_path <- file.path(population_dir, "01_pop_main_panel.rds")

exposure_pooled_output_path <- file.path(
  parent_loss_dir, "01_education_parent_alive_exposure_pooled.rds"
)

# ---- DemoKin input preparation ----

fill_age_series <- function(age, value) {
  non_missing <- which(!is.na(value))

  if (length(non_missing) == 1L) {
    return(rep(value[non_missing], length(value)))
  }

  approx(
    x = age[non_missing],
    y = value[non_missing],
    xout = age,
    method = "linear",
    rule = 2,
    ties = "ordered"
  )$y
}

mx_to_lifetable <- function(mortality_dt) {
  mortality_dt[
    order(age_years),
    {
      mx <- as.numeric(mortality_rate_per_100k) / 100000
      ax <- rep(0.5, .N)
      ax[1] <- min(1, 0.07 + 1.7 * mx[1])

      qx <- mx / (1 + (1 - ax) * mx)
      qx <- pmin(pmax(qx, 0), 1)
      qx[.N] <- 1
      px <- 1 - qx
      px[.N] <- 0

      .(
        age_years = age_years,
        mx = mx,
        qx = qx,
        px = px
      )
    },
    by = .(year, race_eth, sex, educ2)
  ]
}

complete_population_weights <- function(pop_dt) {
  pop_use <- copy(pop_dt)[
    race_eth %chin% RACES_R3 &
      sex %chin% FOCAL_SEXES &
      age_years %in% AGES_MODEL,
    .(
      source,
      year,
      race_eth,
      sex,
      age_years = as.integer(age_years),
      n_persons_weighted = as.numeric(n_persons_weighted)
    )
  ]

  pop_use[, source_expected := "acs"]
  pop_use <- pop_use[source == source_expected]

  pop_agg <- pop_use[
    ,
    .(n_children_pop = sum(n_persons_weighted, na.rm = TRUE)),
    by = .(year, race_eth, age_child = age_years, focal_sex = sex)
  ]

  expected_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    age_child = AGES_MODEL,
    focal_sex = FOCAL_SEXES,
    unique = TRUE
  )

  pop_agg <- pop_agg[expected_grid, on = .(year, race_eth, age_child, focal_sex)]
  setorder(pop_agg, year, race_eth, focal_sex, age_child)
  pop_agg[
    ,
    n_children_pop := {
      values <- as.numeric(n_children_pop)
      values[!is.finite(values) | values <= 0] <- NA_real_
      fill_age_series(age_child, values)
    },
    by = .(year, race_eth, focal_sex)
  ]

  pop_agg[
    ,
    child_sex_weight := n_children_pop / sum(n_children_pop),
    by = .(year, race_eth, age_child)
  ]

  pop_agg[order(year, race_eth, age_child, focal_sex)]
}

prepare_mortality_inputs <- function(mortality_dt) {
  mortality_use <- copy(mortality_dt)[
    rate_definition == "all_cause" &
      race_eth %chin% RACES_R3 &
      sex %chin% FOCAL_SEXES &
      educ2 %chin% EDUC_LEVELS &
      age_years %in% AGES_MODEL &
      year %in% YEARS_MODEL,
    .(
      year,
      race_eth,
      sex,
      age_years = as.integer(age_years),
      educ2,
      mortality_rate_per_100k = as.numeric(mortality_rate_per_100k)
    )
  ]

  expected_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    sex = FOCAL_SEXES,
    age_years = AGES_MODEL,
    educ2 = EDUC_LEVELS,
    unique = TRUE
  )

  mortality_use <- mortality_use[expected_grid, on = .(year, race_eth, sex, age_years, educ2)]

  noba_reference <- mortality_use[
    educ2 == "No_BA",
    .(year, race_eth, sex, age_years, noba_mortality_rate_per_100k = mortality_rate_per_100k)
  ]

  mortality_use <- noba_reference[mortality_use, on = .(year, race_eth, sex, age_years)]
  mortality_use[
    age_years <= 17L &
      educ2 != "No_BA" &
      is.na(mortality_rate_per_100k),
    mortality_rate_per_100k := noba_mortality_rate_per_100k
  ]
  mortality_use[, noba_mortality_rate_per_100k := NULL]

  mortality_use[
    ,
    mortality_rate_per_100k := fill_age_series(age_years, mortality_rate_per_100k),
    by = .(year, race_eth, sex, educ2)
  ]

  mx_to_lifetable(mortality_use)
}

prepare_fertility_inputs <- function(
    fertility_dt,
    rate_definition_id,
    repro_ages,
    allow_existing_fx = FALSE) {
  fertility_rate_input_col <- if ("fertility_rate_per_1000_model" %chin% names(fertility_dt)) {
    "fertility_rate_per_1000_model"
  } else {
    "fertility_rate_per_1000"
  }

  fertility_use <- copy(fertility_dt)[
    rate_definition == rate_definition_id &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS &
      age_years %in% repro_ages &
      year %in% YEARS_MODEL,
    .(
      year,
      race_eth,
      age_years = as.integer(age_years),
      educ2,
      fertility_rate_per_1000 = as.numeric(get(fertility_rate_input_col)),
      fx_input = if ("fx" %chin% names(fertility_dt) && allow_existing_fx) {
        as.numeric(fx)
      } else {
        as.numeric(NA_real_)
      }
    )
  ]

  repro_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    age_years = repro_ages,
    educ2 = EDUC_LEVELS,
    unique = TRUE
  )

  fertility_use <- fertility_use[repro_grid, on = .(year, race_eth, age_years, educ2)]
  fertility_use[
    ,
    fertility_rate_per_1000 := fill_age_series(age_years, fertility_rate_per_1000),
    by = .(year, race_eth, educ2)
  ]
  if (allow_existing_fx) {
    fertility_use[
      ,
      fx_input := fill_age_series(age_years, fx_input),
      by = .(year, race_eth, educ2)
    ]
  }

  full_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    age_years = AGES_MODEL,
    educ2 = EDUC_LEVELS,
    unique = TRUE
  )

  fertility_use <- fertility_use[full_grid, on = .(year, race_eth, age_years, educ2)]
  fertility_use[is.na(fertility_rate_per_1000), fertility_rate_per_1000 := 0]
  if (allow_existing_fx) {
    fertility_use[
      ,
      fx := fifelse(
        !is.na(fx_input),
        fx_input,
        fertility_rate_per_1000 / 1000
      )
    ]
  } else {
    fertility_use[, fx := fertility_rate_per_1000 / 1000]
  }
  fertility_use[, fx_input := NULL]

  fertility_use[order(year, race_eth, age_years, educ2)]
}

prepare_transition_inputs <- function(transition_dt) {
  transition_use <- copy(transition_dt)[
    race_eth %chin% RACES_R3 &
      sex %chin% FOCAL_SEXES &
      age_years %in% AGES_MODEL &
      year %in% YEARS_MODEL,
    .(
      year,
      race_eth,
      sex,
      age_years = as.integer(age_years),
      t_noba_to_noba = as.numeric(t_noba_to_noba),
      t_noba_to_ba = as.numeric(t_noba_to_ba),
      t_ba_to_noba = as.numeric(t_ba_to_noba),
      t_ba_to_ba = as.numeric(t_ba_to_ba)
    )
  ]

  expected_grid <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    sex = FOCAL_SEXES,
    age_years = AGES_MODEL,
    unique = TRUE
  )

  transition_use <- transition_use[expected_grid, on = .(year, race_eth, sex, age_years)]

  transition_use[order(year, race_eth, sex, age_years)]
}

build_year_state_matrix_list <- function(dt, value_col, years, ages, states) {
  out <- vector("list", length(years))
  names(out) <- as.character(years)

  for (i in seq_along(years)) {
    year_i <- years[i]
    wide <- dcast(
      dt[year == year_i],
      age_years ~ educ2,
      value.var = value_col
    )

    if (!identical(names(wide)[-1], states)) {
      setcolorder(wide, c("age_years", states))
    }

    mat <- as.matrix(wide[, ..states])
    storage.mode(mat) <- "double"
    rownames(mat) <- as.character(ages)
    colnames(mat) <- states
    out[[i]] <- mat
  }

  out
}

build_transition_matrix_list <- function(dt, years, ages) {
  out <- vector("list", length(years))
  names(out) <- as.character(years)

  for (i in seq_along(years)) {
    year_i <- years[i]
    by_age <- vector("list", length(ages))
    names(by_age) <- as.character(ages)
    dt_year <- dt[year == year_i][order(age_years)]

    for (j in seq_along(ages)) {
      row <- dt_year[j]
      by_age[[j]] <- matrix(
        c(
          row$t_noba_to_noba, row$t_ba_to_noba,
          row$t_noba_to_ba, row$t_ba_to_ba
        ),
        nrow = length(EDUC_LEVELS),
        ncol = length(EDUC_LEVELS),
        byrow = TRUE,
        dimnames = list(EDUC_LEVELS, EDUC_LEVELS)
      )
    }

    out[[i]] <- by_age
  }

  out
}

extract_parent_exposure <- function(kin_full_dt, race_eth, focal_sex) {
  dt <- copy(as.data.table(kin_full_dt))

  dt <- dt[group == OUTPUT_KIN]
  dt[, parent_educ2 := EDUC_LEVELS[as.integer(as.character(stage_kin))]]
  dt[
    ,
    `:=`(
      year = as.integer(year),
      race_eth = race_eth,
      focal_sex = focal_sex,
      age_child = as.integer(age_focal),
      age_parent = as.integer(age_kin),
      sex_parent = tolower(as.character(sex_kin)),
      child_is_minor = as.integer(age_focal) %in% AGES_MINOR,
      n_parents_alive = as.numeric(count)
    )
  ]

  dt[
    ,
    .(
      year,
      race_eth,
      focal_sex,
      age_child,
      child_is_minor,
      age_parent,
      sex_parent,
      parent_educ2,
      n_parents_alive
    )
  ][order(year, race_eth, focal_sex, age_child, sex_parent, parent_educ2, age_parent)]
}

# ---- Run the kinship model ----

mortality_panel <- as.data.table(readRDS(mortality_input_path))
fertility_panel <- as.data.table(readRDS(fertility_input_path))
transition_panel <- as.data.table(readRDS(transition_input_path))
population_panel <- as.data.table(readRDS(population_input_path))
male_fertility_proxy_panel <- as.data.table(readRDS(male_fertility_proxy_input_path))

mortality_lt <- prepare_mortality_inputs(mortality_panel)

fertility_input_females <- prepare_fertility_inputs(
  fertility_panel,
  rate_definition_id = "mother_based",
  repro_ages = FEMALE_REPRO_AGES
)
fertility_input_males <- prepare_fertility_inputs(
  male_fertility_proxy_panel,
  rate_definition_id = MALE_FERTILITY_RATE_DEFINITION,
  repro_ages = MALE_REPRO_AGES,
  allow_existing_fx = TRUE
)

transition_input <- prepare_transition_inputs(transition_panel)

child_sex_weights <- complete_population_weights(population_panel)

exposure_list <- vector("list", length(RACES_R3) * length(FOCAL_SEXES))
exposure_idx <- 0L

for (race_eth_i in RACES_R3) {

  mortality_race <- mortality_lt[race_eth == race_eth_i]
  fertility_race_females <- fertility_input_females[race_eth == race_eth_i]
  fertility_race_males <- fertility_input_males[race_eth == race_eth_i]
  transition_race <- transition_input[race_eth == race_eth_i]

  U_list_females <- build_year_state_matrix_list(
    mortality_race[sex == "female"],
    value_col = "px",
    years = YEARS_MODEL,
    ages = AGES_MODEL,
    states = EDUC_LEVELS
  )
  U_list_males <- build_year_state_matrix_list(
    mortality_race[sex == "male"],
    value_col = "px",
    years = YEARS_MODEL,
    ages = AGES_MODEL,
    states = EDUC_LEVELS
  )

  F_list_females <- build_year_state_matrix_list(
    fertility_race_females,
    value_col = "fx",
    years = YEARS_MODEL,
    ages = AGES_MODEL,
    states = EDUC_LEVELS
  )
  F_list_males <- build_year_state_matrix_list(
    fertility_race_males,
    value_col = "fx",
    years = YEARS_MODEL,
    ages = AGES_MODEL,
    states = EDUC_LEVELS
  )

  T_list_females <- build_transition_matrix_list(
    transition_race[sex == "female"],
    years = YEARS_MODEL,
    ages = AGES_MODEL
  )
  T_list_males <- build_transition_matrix_list(
    transition_race[sex == "male"],
    years = YEARS_MODEL,
    ages = AGES_MODEL
  )

  for (focal_sex_i in FOCAL_SEXES) {
    kin_full <- DemoKin::kin_multi_stage_time_variant_2sex(
      U_list_females = U_list_females,
      U_list_males = U_list_males,
      F_list_females = F_list_females,
      F_list_males = F_list_males,
      T_list_females = T_list_females,
      T_list_males = T_list_males,
      birth_female = BIRTH_FEMALE,
      parity = FALSE,
      output_kin = OUTPUT_KIN,
      summary_kin = FALSE,
      sex_Focal = fifelse(focal_sex_i == "female", "Female", "Male"),
      initial_stage_Focal = 1,
      output_years = YEARS_MODEL
    )

    exposure_idx <- exposure_idx + 1L
    exposure_list[[exposure_idx]] <- extract_parent_exposure(
      kin_full_dt = kin_full,
      race_eth = race_eth_i,
      focal_sex = focal_sex_i
    )
  }
}

exposure_by_child_sex <- rbindlist(exposure_list, use.names = TRUE)
setorder(
  exposure_by_child_sex,
  year, race_eth, focal_sex, age_child, sex_parent, parent_educ2, age_parent
)

expected_by_sex_grid <- CJ(
  year = YEARS_MODEL,
  race_eth = RACES_R3,
  focal_sex = FOCAL_SEXES,
  age_child = AGES_MODEL,
  age_parent = AGES_MODEL,
  sex_parent = ALLOWED_PARENT_SEX,
  parent_educ2 = EDUC_LEVELS,
  unique = TRUE
)
expected_by_sex_grid[, child_is_minor := age_child %in% AGES_MINOR]

exposure_by_child_sex <- exposure_by_child_sex[
  expected_by_sex_grid,
  on = .(year, race_eth, focal_sex, age_child, age_parent, sex_parent, parent_educ2, child_is_minor)
]

exposure_pooled <- child_sex_weights[
  exposure_by_child_sex,
  on = .(year, race_eth, age_child, focal_sex)
]

exposure_pooled <- exposure_pooled[
  ,
  .(
    child_is_minor = first(child_is_minor),
    n_parents_alive = sum(child_sex_weight * n_parents_alive)
  ),
  by = .(year, race_eth, age_child, age_parent, sex_parent, parent_educ2)
]

setorder(
  exposure_pooled,
  year, race_eth, age_child, sex_parent, parent_educ2, age_parent
)

# ---- Living-parent exposure panels ----

saveRDS(exposure_pooled, exposure_pooled_output_path, compress = "xz")
