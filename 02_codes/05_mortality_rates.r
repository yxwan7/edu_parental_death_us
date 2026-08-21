# ==============================================================================
# 05_mortality_rates.r
#
# Build education-specific all-cause mortality rates after allocating deaths
# with unknown education.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/05_mortality_rates.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Education allocation ----

# Children (0-17) carry no completed education, so their Unknown is bridged to
# No_BA (symmetrically in numerator and denominator); adult Unknown deaths
# (18-85) are donor-allocated.
child_unknown_bridge_max_age <- 17L
unknown_numerator_ages <- 18:85

# ---- Mortality rate files ----

rate_dir <- here("01_data", "02_processed", "rate")
dir.create(rate_dir, recursive = TRUE, showWarnings = FALSE)

death_input_path <- here("01_data", "02_processed", "deaths", "01_death_main_panel.rds")
denominator_input_path <- here(
  "01_data", "02_processed", "population", "02_pop_denominator_death.rds"
)
panel_rds_path <- file.path(rate_dir, "01_education_mortality_rate_panel.rds")

# ---- Mortality allocation functions ----

get_death_adult_age_band <- function(age_years) {
  out <- rep(NA_character_, length(age_years))
  out[age_years >= 18L & age_years <= 24L] <- "18-24"
  out[age_years >= 25L & age_years <= 34L] <- "25-34"
  out[age_years >= 35L & age_years <= 44L] <- "35-44"
  out[age_years >= 45L & age_years <= 54L] <- "45-54"
  out[age_years >= 55L & age_years <= 64L] <- "55-64"
  out[age_years >= 65L & age_years <= 74L] <- "65-74"
  out[age_years >= 75L & age_years <= 85L] <- "75-85"
  out
}

# Denominator: bridge child Unknown, scope to the analysis window.
prepare_denominator <- function(path) {
  dt <- as.data.table(readRDS(path))
  dt[, year := as.integer(year)]
  dt[, age_years := as.integer(age_years)]
  dt[, sex := normalize_missing_character(sex)]
  dt[, race_eth := normalize_missing_character(race_eth)]
  dt[, educ2 := normalize_missing_character(educ2)]
  dt <- bridge_child_unknown_educ2(dt, max_age = child_unknown_bridge_max_age)
  dt[
    year %in% YEARS_MODEL &
      sex %chin% SEX_LEVELS &
      !is.na(age_years) &
      race_eth %chin% RACES_R5 &
      educ2 %chin% EDUC_LEVELS,
    .(
      n_pop_unweighted = sum(n_pop_unweighted, na.rm = TRUE),
      n_pop_weighted = sum(n_pop_weighted, na.rm = TRUE)
    ),
    by = .(year, sex, age_years, race_eth, educ2)
  ][order(year, sex, age_years, race_eth, educ2)]
}

build_all_cause_numerator <- function(death_main_panel) {
  dt <- as.data.table(death_main_panel)
  dt[, year := as.integer(year)]
  dt[, age_years := as.integer(age_years)]
  dt[, sex := normalize_missing_character(sex)]
  dt[, race_eth := normalize_missing_character(race_eth)]
  dt[, educ2 := normalize_missing_character(educ2)]
  dt <- bridge_child_unknown_educ2(dt, max_age = child_unknown_bridge_max_age)

  dt[
    cause_group == "all_cause" &
      year %in% YEARS_MODEL &
      sex %chin% SEX_LEVELS &
      !is.na(age_years) &
      race_eth %chin% RACES_R5,
    .(n_deaths = sum(n_deaths, na.rm = TRUE)),
    by = .(year, sex, age_years, race_eth, educ2)
  ][order(year, sex, age_years, race_eth, educ2)]
}

target_exposure_counts <- function(target_row, denominator_dt) {
  exposure_dt <- denominator_dt[
    year == as.integer(target_row$year) &
      sex == as.character(target_row$sex) &
      age_years == as.integer(target_row$age_years) &
      race_eth == as.character(target_row$race_eth)
  ]
  counts <- sum_col_by_educ(exposure_dt, "n_pop_weighted")
  counts
}

bayes_candidate <- function(rule, event_dt, exposure_dt, threshold) {
  list(
    rule = rule,
    event_dt = event_dt,
    exposure_dt = exposure_dt,
    threshold = threshold
  )
}

score_candidate <- function(candidate, target_exposure, alpha = 0.5) {
  events <- sum_col_by_educ(candidate$event_dt, "n_deaths")
  exposures <- sum_col_by_educ(candidate$exposure_dt, "n_pop_weighted")
  donor_total <- sum(events)
  exposure_total <- sum(exposures)

  if (!is.finite(donor_total) || donor_total < candidate$threshold ||
      !is.finite(exposure_total) || exposure_total <= 0) {
    return(NULL)
  }

  rates <- rep(0, length(EDUC_LEVELS))
  positive_exposure <- exposures > 0
  rates[positive_exposure] <- (events[positive_exposure] + alpha) / exposures[positive_exposure]
  scores <- rates * target_exposure

  if (!is.finite(sum(scores)) || sum(scores) <= 0) {
    return(NULL)
  }

  shares <- scores / sum(scores)
  list(
    counts = events,
    exposures = exposures,
    rates = rates,
    scores = scores,
    shares = shares,
    donor_total = donor_total,
    donor_exposure_total = exposure_total,
    rule = candidate$rule
  )
}

choose_bayes_donor <- function(candidates, target_exposure) {
  for (candidate in candidates) {
    scored <- score_candidate(candidate, target_exposure, alpha = 0.5)
    if (!is.null(scored)) {
      return(scored)
    }
  }

  # Graceful fallback: relax the broadest candidate's count threshold so a rare
  # race x age cell with a small but non-empty donor pool still resolves instead
  # of aborting the whole run.
  fallback <- candidates[[length(candidates)]]
  fallback$threshold <- 0
  scored <- score_candidate(fallback, target_exposure, alpha = 0.5)
  if (!is.null(scored)) {
    return(scored)
  }
}

get_mortality_donors <- function(target_row, known_dt, denominator_dt, target_exposure) {
  year_i <- as.integer(target_row$year)
  sex_i <- as.character(target_row$sex)
  race_i <- as.character(target_row$race_eth)
  age_i <- as.integer(target_row$age_years)
  adult_band <- get_death_adult_age_band(age_i)

  candidates <- list(
    bayes_candidate(
      rule = "same_year_same_race_same_sex_single_age",
      threshold = 30,
      event_dt = known_dt[
        year == year_i &
          sex == sex_i &
          race_eth == race_i &
          age_years == age_i
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          sex == sex_i &
          race_eth == race_i &
          age_years == age_i
      ]
    ),
    bayes_candidate(
      rule = "same_year_same_race_same_sex_age_band",
      threshold = 30,
      event_dt = known_dt[
        year == year_i &
          sex == sex_i &
          race_eth == race_i &
          get_death_adult_age_band(age_years) == adult_band
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          sex == sex_i &
          race_eth == race_i &
          get_death_adult_age_band(age_years) == adult_band
      ]
    ),
    bayes_candidate(
      rule = "centered_3y_same_race_same_sex_age_band",
      threshold = 30,
      event_dt = known_dt[
        year %between% c(year_i - 1L, year_i + 1L) &
          sex == sex_i &
          race_eth == race_i &
          get_death_adult_age_band(age_years) == adult_band
      ],
      exposure_dt = denominator_dt[
        year %between% c(year_i - 1L, year_i + 1L) &
          sex == sex_i &
          race_eth == race_i &
          get_death_adult_age_band(age_years) == adult_band
      ]
    ),
    bayes_candidate(
      rule = "same_year_all_race_same_sex_age_band",
      threshold = 30,
      event_dt = known_dt[
        year == year_i &
          sex == sex_i &
          get_death_adult_age_band(age_years) == adult_band
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          sex == sex_i &
          get_death_adult_age_band(age_years) == adult_band
      ]
    ),
    bayes_candidate(
      rule = "centered_5y_all_race_same_sex_age_band",
      threshold = 1,
      event_dt = known_dt[
        year %between% c(year_i - 2L, year_i + 2L) &
          sex == sex_i &
          get_death_adult_age_band(age_years) == adult_band
      ],
      exposure_dt = denominator_dt[
        year %between% c(year_i - 2L, year_i + 2L) &
          sex == sex_i &
          get_death_adult_age_band(age_years) == adult_band
      ]
    )
  )

  choose_bayes_donor(candidates, target_exposure)
}

allocate_mortality_unknown <- function(unknown_dt, known_dt, denominator_dt) {
  if (!nrow(unknown_dt)) {
    return(data.table())
  }

  rbindlist(
    lapply(seq_len(nrow(unknown_dt)), function(i) {
      target_row <- unknown_dt[i]
      target_exposure <- target_exposure_counts(target_row, denominator_dt)
      donor_info <- get_mortality_donors(
        target_row = target_row,
        known_dt = known_dt,
        denominator_dt = denominator_dt,
        target_exposure = target_exposure
      )
      unknown_total <- as.numeric(target_row$n_deaths)

      out <- target_row[rep(1L, length(EDUC_LEVELS))]
      out[
        ,
        `:=`(
          educ2 = EDUC_LEVELS,
          n_deaths = unknown_total * donor_info$shares,
          donor_rule = donor_info$rule,
          rate_source = donor_info$rule,
          denominator_source = "target_mortality_denominator",
          donor_total = donor_info$donor_total,
          donor_exposure_total = donor_info$donor_exposure_total,
          target_exposure = as.numeric(target_exposure),
          donor_rate = as.numeric(donor_info$rates),
          bayes_score = as.numeric(donor_info$scores),
          donor_share_noba = donor_info$shares[1],
          donor_share_ba = donor_info$shares[2],
          unknown_total = unknown_total
        )
      ]
      out[]
    }),
    use.names = TRUE,
    fill = TRUE
  )
}

# ---- Prepare denominator (no imputation) ----

denominator_dt <- prepare_denominator(denominator_input_path)

# ---- Redistribute numerator Unknown education ----

death_main_panel <- readRDS(death_input_path)
all_cause_numerator <- build_all_cause_numerator(death_main_panel)

known_dt <- all_cause_numerator[
  educ2 %chin% EDUC_LEVELS
]
unknown_dt <- all_cause_numerator[
  educ2 == "Unknown" &
    age_years %in% unknown_numerator_ages
]

# Drop death cells with no matching population exposure: a known cell whose
# (year, sex, age, race, educ2) is absent from the denominator, or an Unknown
# cell whose (year, sex, age, race) target has no denominator exposure to anchor
# the redistribution. Both are retained in the fit record.
known_support_keys <- unique(
  denominator_dt[
    ,
    .(year, sex, age_years, race_eth, educ2)
  ]
)
target_support_keys <- unique(
  denominator_dt[
    ,
    .(year, sex, age_years, race_eth)
  ]
)
known_dt <- known_dt[
  known_support_keys,
  on = .(year, sex, age_years, race_eth, educ2),
  nomatch = 0
]
unknown_dt <- unknown_dt[
  target_support_keys,
  on = .(year, sex, age_years, race_eth),
  nomatch = 0
]

allocations_dt <- allocate_mortality_unknown(
  unknown_dt = unknown_dt,
  known_dt = known_dt,
  denominator_dt = denominator_dt
)

completed_numerator <- rbindlist(
  list(
    known_dt[
      ,
      .(year, sex, age_years, race_eth, educ2, n_deaths)
    ],
    allocations_dt[
      ,
      .(year, sex, age_years, race_eth, educ2, n_deaths)
    ]
  ),
  use.names = TRUE
)[
  ,
  .(n_deaths = sum(n_deaths, na.rm = TRUE)),
  by = .(year, sex, age_years, race_eth, educ2)
][order(year, sex, age_years, race_eth, educ2)]

# ---- Build mortality rate panel ----

mortality_panel <- merge(
  denominator_dt,
  completed_numerator,
  by = c("year", "sex", "age_years", "race_eth", "educ2"),
  all.x = TRUE,
  sort = TRUE
)
mortality_panel[is.na(n_deaths), n_deaths := 0]
mortality_panel[
  ,
  `:=`(
    rate_definition = "all_cause",
    mortality_rate_per_100k = 1e5 * (n_deaths / n_pop_weighted)
  )
]

setcolorder(
  mortality_panel,
  c(
    "rate_definition",
    "year", "sex", "age_years", "race_eth", "educ2",
    "n_deaths", "n_pop_unweighted", "n_pop_weighted",
    "mortality_rate_per_100k"
  )
)

# ---- Mortality rate products ----

saveRDS(mortality_panel, panel_rds_path)
