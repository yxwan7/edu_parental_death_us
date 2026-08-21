# ==============================================================================
# 09_education_transitions.r
#
# Estimate two-state education transitions from adjacent ACS cross-sections and
# education-specific survival.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/09_education_transitions.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Transition definitions ----

# Pairing needs year t+1; the last model year is carried forward from t-1.
YEARS_PAIR <- 2000:2023

# Ages over which q12 is estimated; outside this window the matrix is identity,
# so education is fixed from age 30 on.
AGES_SIGNAL <- 18:29

TRANSITION_METHOD <- "mortality_corrected_synthetic_cohort_monotone_flow"

# Minimum unweighted counts for a share to be estimated at all.
MIN_CELL_UNWEIGHTED <- 100
MIN_ORIGIN_UNWEIGHTED <- 30
SHARE_EPS <- 1e-10

# Survival reaches past AGES_SIGNAL; the extra ages are dropped by the merge
# in build_raw_pairs().
AGES_SURVIVAL <- 18:64

educ_code_map <- c(
  No_BA = "noba",
  BA_plus = "ba"
)

# ---- Transition files ----

rate_dir <- here("01_data", "02_processed", "rate")
dir.create(rate_dir, recursive = TRUE, showWarnings = FALSE)

population_path <- here("01_data", "02_processed", "population", "01_pop_main_panel.rds")
# Survival constraint from the all-cause mortality rate panel built by 05.
mortality_path <- file.path(rate_dir, "01_education_mortality_rate_panel.rds")
transition_panel_rds_path <- file.path(rate_dir, "05_education_transition_panel.rds")

# ---- Transition estimator ----

clamp01 <- function(x) {
  pmin(pmax(x, 0), 1)
}

rate_to_survival <- function(rate_per_100k) {
  exp(-(as.numeric(rate_per_100k) / 100000))
}

# A cell without origin support receives survival 1 and is excluded by the
# origin-support check.
fill_supported_survival <- function(share, origin_n, survival) {
  out <- as.numeric(survival)
  weak_or_no_origin <- !is.finite(origin_n) |
    origin_n < MIN_ORIGIN_UNWEIGHTED |
    !is.finite(share) |
    share <= SHARE_EPS
  out[!is.finite(out) & weak_or_no_origin] <- 1
  out
}

build_population_shares <- function(population_dt) {
  population_use <- copy(population_dt)
  population_use[, source := normalize_missing_character(source)]
  population_use[, sex := normalize_missing_character(sex)]
  population_use[, race_eth := normalize_missing_character(race_eth)]
  # Unknown education carries no state information -> NA (dropped by the
  # level filter below).
  population_use[, educ2 := normalize_missing_character(educ2)]
  population_use[!educ2 %chin% EDUC_LEVELS, educ2 := NA_character_]
  population_use[, year := as.integer(year)]
  population_use[, age_years := as.integer(age_years)]

  population_use <- population_use[
    source == "acs" &
      year %in% YEARS_MODEL &
      age_years %in% AGES_EDUC_TRANSITION &
      sex %chin% SEX_LEVELS &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS
  ]
  population_use[, educ_code := educ_code_map[educ2]]

  population_agg <- population_use[
    ,
    .(
      n_persons_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
      n_persons_weighted = sum(n_persons_weighted, na.rm = TRUE)
    ),
    by = .(source, year, age_years, sex, race_eth, educ_code)
  ]

  population_totals <- population_agg[
    ,
    .(
      cell_total_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
      cell_total_weighted = sum(n_persons_weighted, na.rm = TRUE)
    ),
    by = .(source, year, age_years, sex, race_eth)
  ]

  population_long <- merge(
    population_agg,
    population_totals,
    by = c("source", "year", "age_years", "sex", "race_eth"),
    all.x = TRUE
  )
  population_long[
    ,
    share := fifelse(
      cell_total_weighted > 0,
      n_persons_weighted / cell_total_weighted,
      NA_real_
    )
  ]

  population_wide <- dcast(
    population_long,
    source + year + age_years + sex + race_eth ~ educ_code,
    value.var = c("share", "n_persons_unweighted", "n_persons_weighted"),
    fill = 0
  )

  expected_names <- c(
    "share_noba", "share_ba",
    "n_persons_unweighted_noba", "n_persons_unweighted_ba",
    "n_persons_weighted_noba", "n_persons_weighted_ba"
  )
  missing_expected <- setdiff(expected_names, names(population_wide))
  if (length(missing_expected)) {
    for (col in missing_expected) {
      population_wide[, (col) := 0]
    }
  }

  setnames(
    population_wide,
    old = c(
      "n_persons_unweighted_noba", "n_persons_unweighted_ba",
      "n_persons_weighted_noba", "n_persons_weighted_ba"
    ),
    new = c(
      "n_noba_unweighted", "n_ba_unweighted",
      "n_noba_weighted", "n_ba_weighted"
    ),
    skip_absent = TRUE
  )

  population_wide <- merge(
    population_wide,
    population_totals,
    by = c("source", "year", "age_years", "sex", "race_eth"),
    all.x = TRUE
  )

  setorder(population_wide, source, year, sex, race_eth, age_years)
  population_wide[]
}

build_survival_dt <- function(mortality_dt) {
  mortality_use <- copy(mortality_dt)
  mortality_use[, sex := normalize_missing_character(sex)]
  mortality_use[, race_eth := normalize_missing_character(race_eth)]
  mortality_use[, educ2 := normalize_missing_character(educ2)]
  mortality_use[, year := as.integer(year)]
  mortality_use[, age_years := as.integer(age_years)]

  mortality_use <- mortality_use[
    year %in% YEARS_PAIR &
      age_years %in% AGES_SURVIVAL &
      sex %chin% SEX_LEVELS &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS
  ]
  mortality_use[, educ_code := educ_code_map[educ2]]
  mortality_use[, survival := rate_to_survival(mortality_rate_per_100k)]

  survival_dt <- dcast(
    mortality_use,
    year + age_years + sex + race_eth ~ educ_code,
    value.var = "survival",
    fill = NA_real_
  )

  for (col in c("noba", "ba")) {
    if (!col %chin% names(survival_dt)) {
      survival_dt[, (col) := NA_real_]
    }
  }
  setnames(
    survival_dt,
    old = c("noba", "ba"),
    new = c("survival_noba", "survival_ba"),
    skip_absent = TRUE
  )
  setcolorder(
    survival_dt,
    c(
      "year", "age_years", "sex", "race_eth",
      "survival_noba", "survival_ba"
    )
  )
  setorder(survival_dt, year, sex, race_eth, age_years)
  survival_dt[]
}

build_raw_pairs <- function(population_wide, survival_dt) {
  current_dt <- copy(
    population_wide[
      source == "acs" &
        year %in% YEARS_PAIR &
        age_years %in% AGES_SIGNAL
    ]
  )

  next_dt <- copy(population_wide[source == "acs"])
  next_dt[, `:=`(year = year - 1L, age_years = age_years - 1L)]
  setnames(
    next_dt,
    old = c(
      "share_noba", "share_ba",
      "n_noba_unweighted", "n_ba_unweighted",
      "n_noba_weighted", "n_ba_weighted",
      "cell_total_unweighted", "cell_total_weighted"
    ),
    new = c(
      "share_noba_next", "share_ba_next",
      "n_noba_unweighted_next", "n_ba_unweighted_next",
      "n_noba_weighted_next", "n_ba_weighted_next",
      "cell_total_unweighted_next", "cell_total_weighted_next"
    )
  )

  raw_pairs <- merge(
    current_dt,
    next_dt[
      ,
      .(
        source,
        year,
        age_years,
        sex,
        race_eth,
        share_noba_next,
        share_ba_next,
        n_noba_unweighted_next,
        n_ba_unweighted_next,
        n_noba_weighted_next,
        n_ba_weighted_next,
        cell_total_unweighted_next,
        cell_total_weighted_next
      )
    ],
    by = c("source", "year", "age_years", "sex", "race_eth"),
    all = FALSE
  )
  raw_pairs[, source_component := "acs"]

  raw_pairs <- merge(
    raw_pairs,
    survival_dt,
    by = c("year", "age_years", "sex", "race_eth"),
    all.x = TRUE
  )

  raw_pairs[]
}

estimate_monotone_raw <- function(raw_pairs) {
  out <- copy(raw_pairs)

  out[
    ,
    `:=`(
      survival_noba_used = fill_supported_survival(
        share_noba,
        n_noba_unweighted,
        survival_noba
      ),
      survival_ba_used = fill_supported_survival(
        share_ba,
        n_ba_unweighted,
        survival_ba
      )
    )
  ]

  out[
    ,
    survival_blocking_missing :=
      (!is.finite(survival_noba_used) &
        n_noba_unweighted >= MIN_ORIGIN_UNWEIGHTED &
        share_noba > SHARE_EPS) |
      (!is.finite(survival_ba_used) &
        n_ba_unweighted >= MIN_ORIGIN_UNWEIGHTED &
        share_ba > SHARE_EPS)
  ]

  out[
    ,
    raw_inputs_ok :=
      is.finite(cell_total_unweighted) &
      cell_total_unweighted >= MIN_CELL_UNWEIGHTED &
      is.finite(cell_total_unweighted_next) &
      cell_total_unweighted_next >= MIN_CELL_UNWEIGHTED &
      is.finite(share_noba) &
      is.finite(share_ba) &
      is.finite(share_noba_next) &
      is.finite(share_ba_next) &
      !survival_blocking_missing
  ]

  # Roll the origin cross-section forward one year of education-specific
  # mortality, so the share change that remains is transition, not selection.
  out[
    ,
    `:=`(
      survivor_noba = share_noba * survival_noba_used,
      survivor_ba = share_ba * survival_ba_used
    )
  ]
  out[, survivor_total := survivor_noba + survivor_ba]
  out[
    raw_inputs_ok & is.finite(survivor_total) & survivor_total > 0,
    `:=`(
      r_noba = survivor_noba / survivor_total,
      r_ba = survivor_ba / survivor_total
    )
  ]

  out[
    ,
    q12_eligible :=
      raw_inputs_ok &
      n_noba_unweighted >= MIN_ORIGIN_UNWEIGHTED &
      is.finite(r_noba) &
      r_noba > SHARE_EPS
  ]

  # Two-state monotone flow: BA_plus is absorbing, so the rise in its
  # mortality-corrected share across the paired cross-sections identifies the
  # single upward flow exactly. The floor at zero is structural, not a smoothing
  # choice -- a degree cannot be un-earned, so a measured decline (sampling
  # noise or net out-migration) is recorded as no upward flow.
  out[
    q12_eligible == TRUE,
    q12_raw := clamp01(pmax(share_ba_next - r_ba, 0) / r_noba)
  ]
  out[q12_eligible != TRUE, q12_raw := NA_real_]

  out[
    ,
    raw_status := fcase(
      !raw_inputs_ok, "insufficient_inputs",
      !q12_eligible, "insufficient_origin_support",
      default = "ok"
    )
  ]

  setorder(out, year, sex, race_eth, age_years)
  out[]
}

build_transition_panel <- function(q_estimates) {
  full_grid <- CJ(
    year = YEARS_MODEL,
    sex = SEX_LEVELS,
    race_eth = RACES_R3,
    age_years = AGES_EDUC_TRANSITION,
    unique = TRUE
  )

  panel <- merge(
    full_grid,
    q_estimates,
    by = c("year", "sex", "race_eth", "age_years"),
    all.x = TRUE
  )

  panel[is.na(q12), q12 := 0]
  panel[age_years < min(AGES_SIGNAL) | age_years > max(AGES_SIGNAL), q12 := 0]

  panel[, `:=`(
    t_noba_to_noba = 1 - q12,
    t_noba_to_ba = q12,
    t_ba_to_noba = 0,
    t_ba_to_ba = 1,
    transition_method = TRANSITION_METHOD
  )]

  setcolorder(
    panel,
    c(
      "year", "sex", "race_eth", "age_years",
      "q12",
      "t_noba_to_noba", "t_noba_to_ba",
      "t_ba_to_noba", "t_ba_to_ba",
      "transition_method", "source_component"
    )
  )

  setorder(panel, year, sex, race_eth, age_years)
  panel[]
}

# ---- Estimate the transition panel ----

population_dt <- as.data.table(readRDS(population_path))
mortality_dt <- as.data.table(readRDS(mortality_path))

population_shares <- build_population_shares(population_dt)
survival_dt <- build_survival_dt(mortality_dt)

raw_pairs <- build_raw_pairs(population_shares, survival_dt)
raw_estimates <- estimate_monotone_raw(raw_pairs)

# Used as read. Two states and an absorbing target make q12 point-identified by
# the paired cross-sections. Sampling noise is carried by re-estimating the
# panel inside each uncertainty draw (see 14), not by damping it out of the
# point estimate.
q_estimates <- raw_estimates[
  ,
  .(
    year,
    sex,
    race_eth,
    age_years,
    q12 = clamp01(q12_raw),
    source_component = "acs_monotone_flow"
  )
]
q_estimates[!is.finite(q12), `:=`(q12 = 0, source_component = "no_origin_support")]

# The final model year has no t+1 cross-section to pair with.
q_last_year <- copy(q_estimates[year == max(YEARS_PAIR)])
q_last_year[, `:=`(
  year = max(YEARS_MODEL),
  source_component = paste0("copied_from_", max(YEARS_PAIR))
)]

q_estimates <- rbindlist(list(q_estimates, q_last_year), use.names = TRUE)

transition_panel <- build_transition_panel(q_estimates)

# ---- Transition panel ----

saveRDS(transition_panel, transition_panel_rds_path)
