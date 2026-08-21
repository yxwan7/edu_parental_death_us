# ==============================================================================
# 07_fertility_rates.r
#
# Build mother-based fertility rates by education, including unknown-education
# allocation and small-cell stabilization.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/07_fertility_rates.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Fertility allocation settings ----

# Children (10-17) carry no completed education, so their Unknown is bridged to
# No_BA (symmetrically in numerator and denominator); adult Unknown births
# (18-54) are donor-allocated.
child_unknown_bridge_max_age <- 17L
unknown_numerator_ages <- 18:54

# Mother education on the birth certificate is missing far above baseline during
# the staggered state rollout of the 2003 certificate revision (observed Unknown
# share: 32.9% in 2009, falling to 2.4-3.0% by 2014-2015, back to the ~1.3%
# baseline from 2016). Missingness is state-driven, i.e. not random, so within
# this window the observed known-education composition is polluted and Unknown
# births are instead allocated from education-specific rates log-interpolated
# between the stable anchor years on either side.
transition_high_missing_years <- 2009:2015
transition_anchor_years <- c(2000:2008, 2016:2024)

# Small-cell reliability thresholds for the assembled panel: cells with an
# extreme direct rate on negligible denominator support get a donor age-band
# rate in fertility_rate_per_1000_model; low-support cells are review-flagged.
FERTILITY_REVIEW_UNWEIGHTED_MIN <- 20
FERTILITY_REVIEW_WEIGHTED_MIN <- 1000
FERTILITY_STABILIZE_RATE_PER_1000_MIN <- 500
FERTILITY_STABILIZE_UNWEIGHTED_MIN <- 5
FERTILITY_STABILIZE_WEIGHTED_MIN <- 500
FERTILITY_DONOR_BIRTHS_MIN <- 50
FERTILITY_DONOR_UNWEIGHTED_MIN <- 50

# ---- Fertility rate files ----

rate_dir <- here("01_data", "02_processed", "rate")
dir.create(rate_dir, recursive = TRUE, showWarnings = FALSE)

birth_input_path <- here("01_data", "02_processed", "births", "01_birth_main_panel.rds")
denominator_input_path <- here(
  "01_data", "02_processed", "population", "03_pop_denominator_birth.rds"
)
panel_rds_path <- file.path(rate_dir, "03_education_fertility_rate_panel.rds")

# ---- Fertility allocation functions ----

get_fertility_age_band <- function(age_years) {
  out <- rep(NA_character_, length(age_years))
  out[age_years >= 10L & age_years <= 14L] <- "10-14"
  out[age_years >= 15L & age_years <= 19L] <- "15-19"
  out[age_years >= 20L & age_years <= 24L] <- "20-24"
  out[age_years >= 25L & age_years <= 29L] <- "25-29"
  out[age_years >= 30L & age_years <= 34L] <- "30-34"
  out[age_years >= 35L & age_years <= 39L] <- "35-39"
  out[age_years >= 40L & age_years <= 44L] <- "40-44"
  out[age_years >= 45L & age_years <= 49L] <- "45-49"
  out[age_years >= 50L & age_years <= 54L] <- "50-54"
  out
}

# Denominator: bridge child Unknown, scope to the analysis window. The ACS
# female denominator is treated as fully observed exposure (Sasson, 2016, does
# not impute the at-risk denominator); adult Unknown education is dropped and
# its weighted size is retained in the fit record.
prepare_denominator <- function(path) {
  dt <- as.data.table(readRDS(path))
  dt[, year := as.integer(year)]
  dt[, age_years := as.integer(age_years)]
  dt[, race_eth := normalize_missing_character(race_eth)]
  dt[, educ2 := normalize_missing_character(educ2)]
  dt <- bridge_child_unknown_educ2(dt, max_age = child_unknown_bridge_max_age)
  dt[
    year %in% YEARS_MODEL &
      age_years %in% AGES_FEMALE_REPRO &
      race_eth %chin% RACES_R5 &
      educ2 %chin% EDUC_LEVELS,
    .(
      n_women_unweighted = sum(n_women_unweighted, na.rm = TRUE),
      n_women_weighted = sum(n_women_weighted, na.rm = TRUE)
    ),
    by = .(year, age_years, race_eth, educ2)
  ][order(year, age_years, race_eth, educ2)]
}

# Collapse the birth main panel across father characteristics into one
# mother-based numerator per year x age x race x educ2 cell.
build_mother_based_numerator <- function(birth_main_panel) {
  dt <- as.data.table(birth_main_panel)
  dt[, year := as.integer(year)]
  dt[, resident_birth_flag := as.logical(resident_birth_flag)]
  dt[, mother_age_years := as.integer(mother_age_years)]
  dt[, mother_race_eth := normalize_missing_character(mother_race_eth)]
  dt[, mother_educ2 := normalize_missing_character(mother_educ2)]
  dt <- bridge_child_unknown_educ2(
    dt,
    age_col = "mother_age_years",
    educ_col = "mother_educ2",
    max_age = child_unknown_bridge_max_age
  )

  dt[
    resident_birth_flag %in% TRUE &
      year %in% YEARS_MODEL &
      mother_age_years %in% AGES_FEMALE_REPRO &
      mother_race_eth %chin% RACES_R5,
    .(n_births = sum(n_births, na.rm = TRUE)),
    by = .(
      year,
      age_years = mother_age_years,
      race_eth = mother_race_eth,
      educ2 = mother_educ2
    )
  ][order(year, age_years, race_eth, educ2)]
}

target_exposure_counts <- function(target_row, denominator_dt) {
  exposure_dt <- denominator_dt[
    year == as.integer(target_row$year) &
      age_years == as.integer(target_row$age_years) &
      race_eth == as.character(target_row$race_eth)
  ]
  counts <- sum_col_by_educ(exposure_dt, "n_women_weighted")
  counts
}

bayes_candidate <- function(rule, event_dt, exposure_dt, threshold) {
  list(
    type = "event_rate",
    rule = rule,
    event_dt = event_dt,
    exposure_dt = exposure_dt,
    threshold = threshold
  )
}

# Transition-window donor: education-specific fertility rates estimated in each
# anchor year, log-linearly interpolated to the target year (rates are strictly
# positive, so the interpolation runs on the log scale; for 2009-2015 targets
# this bridges the 2008 and 2016 anchors).
anchor_lograte_candidate <- function(
    rule,
    event_dt,
    exposure_dt,
    target_year,
    threshold) {
  events_by <- event_dt[
    year %in% transition_anchor_years,
    .(events = sum(n_births, na.rm = TRUE)),
    by = .(year, educ2)
  ]
  exposure_by <- exposure_dt[
    year %in% transition_anchor_years,
    .(exposure = sum(n_women_weighted, na.rm = TRUE)),
    by = .(year, educ2)
  ]
  grid <- CJ(year = transition_anchor_years, educ2 = EDUC_LEVELS)
  grid <- merge(grid, events_by, by = c("year", "educ2"), all.x = TRUE, sort = FALSE)
  grid <- merge(grid, exposure_by, by = c("year", "educ2"), all.x = TRUE, sort = FALSE)
  grid[is.na(events), events := 0]
  grid[is.na(exposure), exposure := 0]

  donor_total <- grid[, sum(events, na.rm = TRUE)]
  donor_exposure_total <- grid[, sum(exposure, na.rm = TRUE)]
  if (!is.finite(donor_total) || donor_total < threshold ||
      !is.finite(donor_exposure_total) || donor_exposure_total <= 0) {
    return(NULL)
  }

  rates <- setNames(rep(NA_real_, length(EDUC_LEVELS)), EDUC_LEVELS)
  for (educ in EDUC_LEVELS) {
    sub <- grid[educ2 == educ & exposure > 0]
    if (!nrow(sub)) {
      next
    }
    sub[, rate := (events + 0.5) / exposure]
    sub <- sub[is.finite(rate) & rate > 0]
    if (!nrow(sub)) {
      next
    }
    if (nrow(sub) == 1L) {
      rates[[educ]] <- sub$rate[[1]]
    } else {
      rates[[educ]] <- exp(
        stats::approx(
          x = sub$year,
          y = log(sub$rate),
          xout = target_year,
          method = "linear",
          rule = 2,
          ties = "ordered"
        )$y
      )
    }
  }

  rates[is.na(rates) | !is.finite(rates) | rates < 0] <- 0
  if (!is.finite(sum(rates)) || sum(rates) <= 0) {
    return(NULL)
  }

  list(
    type = "anchor_lograte",
    rule = rule,
    rates = as.numeric(rates),
    threshold = threshold,
    donor_total = donor_total,
    donor_exposure_total = donor_exposure_total
  )
}

score_candidate <- function(candidate, target_exposure, alpha = 0.5) {
  if (identical(candidate$type, "anchor_lograte")) {
    rates <- candidate$rates
    donor_total <- candidate$donor_total
    donor_exposure_total <- candidate$donor_exposure_total
  } else {
    events <- sum_col_by_educ(candidate$event_dt, "n_births")
    exposures <- sum_col_by_educ(candidate$exposure_dt, "n_women_weighted")
    donor_total <- sum(events)
    donor_exposure_total <- sum(exposures)

    if (!is.finite(donor_total) || donor_total < candidate$threshold ||
        !is.finite(donor_exposure_total) || donor_exposure_total <= 0) {
      return(NULL)
    }

    rates <- rep(0, length(EDUC_LEVELS))
    positive_exposure <- exposures > 0
    rates[positive_exposure] <- (events[positive_exposure] + alpha) / exposures[positive_exposure]
  }

  scores <- rates * target_exposure
  if (!is.finite(sum(scores)) || sum(scores) <= 0) {
    return(NULL)
  }

  shares <- scores / sum(scores)
  list(
    rates = rates,
    scores = scores,
    shares = shares,
    donor_total = donor_total,
    donor_exposure_total = donor_exposure_total,
    rule = candidate$rule
  )
}

choose_bayes_donor <- function(candidates, target_exposure) {
  candidates <- candidates[!vapply(candidates, is.null, logical(1))]
  for (candidate in candidates) {
    scored <- score_candidate(candidate, target_exposure, alpha = 0.5)
    if (!is.null(scored)) {
      return(scored)
    }
  }

  # Graceful fallback: relax the broadest candidate's count threshold so a rare
  # race x age cell with a small but non-empty donor pool still resolves instead
  # of aborting the whole run.
  if (length(candidates) > 0L) {
    fallback <- candidates[[length(candidates)]]
    if (identical(fallback$type, "event_rate")) {
      fallback$threshold <- 0
      scored <- score_candidate(fallback, target_exposure, alpha = 0.5)
      if (!is.null(scored)) {
        return(scored)
      }
    }
  }

}

get_regular_fertility_donors <- function(target_row, known_dt, denominator_dt, target_exposure) {
  year_i <- as.integer(target_row$year)
  race_i <- as.character(target_row$race_eth)
  age_i <- as.integer(target_row$age_years)
  age_band_i <- get_fertility_age_band(age_i)

  candidates <- list(
    bayes_candidate(
      rule = "same_year_same_race_single_age",
      threshold = 50,
      event_dt = known_dt[
        year == year_i &
          race_eth == race_i &
          age_years == age_i
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          race_eth == race_i &
          age_years == age_i
      ]
    ),
    bayes_candidate(
      rule = "same_year_same_race_age_band",
      threshold = 50,
      event_dt = known_dt[
        year == year_i &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ]
    ),
    bayes_candidate(
      rule = "centered_3y_same_race_age_band",
      threshold = 50,
      event_dt = known_dt[
        year %between% c(year_i - 1L, year_i + 1L) &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ],
      exposure_dt = denominator_dt[
        year %between% c(year_i - 1L, year_i + 1L) &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ]
    ),
    bayes_candidate(
      rule = "same_year_all_race_age_band",
      threshold = 1,
      event_dt = known_dt[
        year == year_i &
          get_fertility_age_band(age_years) == age_band_i
      ],
      exposure_dt = denominator_dt[
        year == year_i &
          get_fertility_age_band(age_years) == age_band_i
      ]
    )
  )

  choose_bayes_donor(candidates, target_exposure)
}

get_transition_fertility_donors <- function(target_row, known_dt, denominator_dt, target_exposure) {
  year_i <- as.integer(target_row$year)
  race_i <- as.character(target_row$race_eth)
  age_i <- as.integer(target_row$age_years)
  age_band_i <- get_fertility_age_band(age_i)

  candidates <- list(
    anchor_lograte_candidate(
      rule = "anchor_lograte_same_race_single_age",
      threshold = 50,
      target_year = year_i,
      event_dt = known_dt[
        year %in% transition_anchor_years &
          race_eth == race_i &
          age_years == age_i
      ],
      exposure_dt = denominator_dt[
        year %in% transition_anchor_years &
          race_eth == race_i &
          age_years == age_i
      ]
    ),
    anchor_lograte_candidate(
      rule = "anchor_lograte_same_race_age_band",
      threshold = 50,
      target_year = year_i,
      event_dt = known_dt[
        year %in% transition_anchor_years &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ],
      exposure_dt = denominator_dt[
        year %in% transition_anchor_years &
          race_eth == race_i &
          get_fertility_age_band(age_years) == age_band_i
      ]
    ),
    anchor_lograte_candidate(
      rule = "anchor_lograte_all_race_age_band",
      threshold = 1,
      target_year = year_i,
      event_dt = known_dt[
        year %in% transition_anchor_years &
          get_fertility_age_band(age_years) == age_band_i
      ],
      exposure_dt = denominator_dt[
        year %in% transition_anchor_years &
          get_fertility_age_band(age_years) == age_band_i
      ]
    )
  )

  choose_bayes_donor(candidates, target_exposure)
}

get_fertility_donors <- function(target_row, known_dt, denominator_dt, target_exposure) {
  if (as.integer(target_row$year) %in% transition_high_missing_years) {
    return(get_transition_fertility_donors(
      target_row,
      known_dt,
      denominator_dt,
      target_exposure
    ))
  }

  get_regular_fertility_donors(
    target_row,
    known_dt,
    denominator_dt,
    target_exposure
  )
}

allocate_fertility_unknown <- function(unknown_dt, known_dt, denominator_dt) {
  if (!nrow(unknown_dt)) {
    return(data.table())
  }

  rbindlist(
    lapply(seq_len(nrow(unknown_dt)), function(i) {
      target_row <- unknown_dt[i]
      target_exposure <- target_exposure_counts(target_row, denominator_dt)
      donor_info <- get_fertility_donors(
        target_row = target_row,
        known_dt = known_dt,
        denominator_dt = denominator_dt,
        target_exposure = target_exposure
      )
      unknown_total <- as.numeric(target_row$n_births)

      out <- target_row[rep(1L, length(EDUC_LEVELS))]
      out[
        ,
        `:=`(
          educ2 = EDUC_LEVELS,
          n_births = unknown_total * donor_info$shares,
          donor_rule = donor_info$rule,
          rate_source = donor_info$rule,
          denominator_source = "target_fertility_denominator",
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

summarize_fertility_donor <- function(donor_dt, rule, strict = TRUE) {
  donor_births <- donor_dt[, sum(n_births, na.rm = TRUE)]
  donor_women_weighted <- donor_dt[, sum(n_women_weighted, na.rm = TRUE)]
  donor_women_unweighted <- donor_dt[, sum(n_women_unweighted, na.rm = TRUE)]

  if (!is.finite(donor_women_weighted) || donor_women_weighted <= 0) {
    return(NULL)
  }
  if (strict) {
    if (!is.finite(donor_births) || donor_births < FERTILITY_DONOR_BIRTHS_MIN ||
        !is.finite(donor_women_unweighted) ||
        donor_women_unweighted < FERTILITY_DONOR_UNWEIGHTED_MIN) {
      return(NULL)
    }
  }

  list(
    rule = rule,
    donor_births = donor_births,
    donor_women_weighted = donor_women_weighted,
    donor_women_unweighted = donor_women_unweighted,
    donor_rate_per_1000 = 1000 * donor_births / donor_women_weighted
  )
}

choose_fertility_small_cell_donor <- function(target_row, donor_pool) {
  year_i <- as.integer(target_row$year)
  race_i <- as.character(target_row$race_eth)
  educ_i <- as.character(target_row$educ2)
  age_band_i <- as.character(target_row$fertility_age_band)

  candidates <- list(
    list(
      rule = "same_year_same_race_same_educ_age_band_excluding_flagged",
      strict = TRUE,
      dt = donor_pool[
        year == year_i &
          race_eth == race_i &
          educ2 == educ_i &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "centered_3y_same_race_same_educ_age_band_excluding_flagged",
      strict = TRUE,
      dt = donor_pool[
        year %between% c(year_i - 1L, year_i + 1L) &
          race_eth == race_i &
          educ2 == educ_i &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "same_year_all_race_same_educ_age_band_excluding_flagged",
      strict = TRUE,
      dt = donor_pool[
        year == year_i &
          educ2 == educ_i &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "same_year_all_race_all_educ_age_band_excluding_flagged",
      strict = TRUE,
      dt = donor_pool[
        year == year_i &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "same_year_all_race_all_educ_age_band_relaxed_positive_denominator",
      strict = FALSE,
      dt = donor_pool[
        year == year_i &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "centered_3y_all_race_all_educ_age_band_relaxed_positive_denominator",
      strict = FALSE,
      dt = donor_pool[
        year %between% c(year_i - 1L, year_i + 1L) &
          fertility_age_band == age_band_i
      ]
    ),
    list(
      rule = "all_year_all_race_all_educ_age_band_relaxed_positive_denominator",
      strict = FALSE,
      dt = donor_pool[
        fertility_age_band == age_band_i
      ]
    )
  )

  for (candidate in candidates) {
    donor_info <- summarize_fertility_donor(
      candidate$dt,
      rule = candidate$rule,
      strict = candidate$strict
    )
    if (!is.null(donor_info)) {
      return(donor_info)
    }
  }

}

add_fertility_small_cell_reliability <- function(fertility_panel) {
  dt <- data.table::copy(data.table::as.data.table(fertility_panel))
  dt[, row_id__ := .I]
  dt[, fertility_age_band := get_fertility_age_band(as.integer(age_years))]
  dt[
    ,
    fertility_rate_per_1000_direct := as.numeric(fertility_rate_per_1000)
  ]
  dt[
    ,
    small_cell_review_flag :=
      as.numeric(n_women_unweighted) < FERTILITY_REVIEW_UNWEIGHTED_MIN |
      as.numeric(n_women_weighted) < FERTILITY_REVIEW_WEIGHTED_MIN
  ]
  dt[
    ,
    small_cell_stabilized_flag :=
      fertility_rate_per_1000_direct > FERTILITY_STABILIZE_RATE_PER_1000_MIN &
      (
        as.numeric(n_women_unweighted) < FERTILITY_STABILIZE_UNWEIGHTED_MIN |
          as.numeric(n_women_weighted) < FERTILITY_STABILIZE_WEIGHTED_MIN
      )
  ]
  dt[
    ,
    `:=`(
      fertility_rate_per_1000_model = fertility_rate_per_1000_direct,
      fertility_rate_model_status = "direct",
      small_cell_reason = data.table::fifelse(
        small_cell_stabilized_flag,
        "extreme_direct_rate_and_low_denominator_support",
        data.table::fifelse(
          small_cell_review_flag,
          "low_denominator_support_review",
          NA_character_
        )
      ),
      small_cell_rule = NA_character_,
      small_cell_donor_rate_per_1000 = NA_real_,
      small_cell_donor_births = NA_real_,
      small_cell_donor_women_weighted = NA_real_,
      small_cell_donor_women_unweighted = NA_real_
    )
  ]

  stabilized_ids <- dt[small_cell_stabilized_flag == TRUE, row_id__]
  if (length(stabilized_ids) > 0L) {
    donor_pool <- dt[
      small_cell_stabilized_flag == FALSE &
        !is.na(fertility_age_band) &
        is.finite(as.numeric(n_women_weighted)) &
        as.numeric(n_women_weighted) > 0
    ]

    for (row_id_i in stabilized_ids) {
      donor_info <- choose_fertility_small_cell_donor(
        dt[row_id__ == row_id_i],
        donor_pool
      )
      dt[
        row_id__ == row_id_i,
        `:=`(
          fertility_rate_per_1000_model = donor_info$donor_rate_per_1000,
          fertility_rate_model_status = "stabilized_small_cell",
          small_cell_rule = donor_info$rule,
          small_cell_donor_rate_per_1000 = donor_info$donor_rate_per_1000,
          small_cell_donor_births = donor_info$donor_births,
          small_cell_donor_women_weighted = donor_info$donor_women_weighted,
          small_cell_donor_women_unweighted = donor_info$donor_women_unweighted
        )
      ]
    }
  }

  dt[, c("row_id__", "fertility_age_band") := NULL]
  dt[]
}

# ---- Prepare denominator (no imputation) ----

denominator_dt <- prepare_denominator(denominator_input_path)

# ---- Redistribute numerator Unknown education ----

birth_main_panel <- readRDS(birth_input_path)
mother_based_numerator <- build_mother_based_numerator(birth_main_panel)

known_dt <- mother_based_numerator[
  educ2 %chin% EDUC_LEVELS
]
unknown_dt <- mother_based_numerator[
  educ2 == "Unknown" &
    age_years %in% unknown_numerator_ages
]

# Drop birth cells with no matching population exposure: a known cell whose
# (year, age, race, educ2) is absent from the denominator, or an Unknown cell
# whose (year, age, race) target has no denominator exposure to anchor the
# redistribution. Both are retained in the fit record.
known_support_keys <- unique(
  denominator_dt[
    ,
    .(year, age_years, race_eth, educ2)
  ]
)
target_support_keys <- unique(
  denominator_dt[
    ,
    .(year, age_years, race_eth)
  ]
)
known_dt <- known_dt[
  known_support_keys,
  on = .(year, age_years, race_eth, educ2),
  nomatch = 0
]
unknown_dt <- unknown_dt[
  target_support_keys,
  on = .(year, age_years, race_eth),
  nomatch = 0
]

allocations_dt <- allocate_fertility_unknown(
  unknown_dt = unknown_dt,
  known_dt = known_dt,
  denominator_dt = denominator_dt
)

completed_numerator <- rbindlist(
  list(
    known_dt[
      ,
      .(year, age_years, race_eth, educ2, n_births)
    ],
    allocations_dt[
      ,
      .(year, age_years, race_eth, educ2, n_births)
    ]
  ),
  use.names = TRUE
)[
  ,
  .(n_births = sum(n_births, na.rm = TRUE)),
  by = .(year, age_years, race_eth, educ2)
][order(year, age_years, race_eth, educ2)]

# ---- Build fertility rate panel with small-cell stabilization ----

fertility_panel <- merge(
  denominator_dt,
  completed_numerator,
  by = c("year", "age_years", "race_eth", "educ2"),
  all.x = TRUE,
  sort = TRUE
)
fertility_panel[is.na(n_births), n_births := 0]
fertility_panel[
  ,
  `:=`(
    rate_definition = "mother_based",
    fertility_rate_per_1000 = 1e3 * (n_births / n_women_weighted)
  )
]

fertility_panel <- add_fertility_small_cell_reliability(fertility_panel)
setcolorder(
  fertility_panel,
  c(
    "rate_definition",
    "year", "age_years", "race_eth", "educ2",
    "n_births", "n_women_unweighted", "n_women_weighted",
    "fertility_rate_per_1000",
    "fertility_rate_per_1000_direct",
    "fertility_rate_per_1000_model",
    "fertility_rate_model_status",
    "small_cell_review_flag",
    "small_cell_stabilized_flag",
    "small_cell_reason",
    "small_cell_rule",
    "small_cell_donor_rate_per_1000",
    "small_cell_donor_births",
    "small_cell_donor_women_weighted",
    "small_cell_donor_women_unweighted"
  )
)

# ---- Fertility rate products ----

saveRDS(fertility_panel, panel_rds_path)
