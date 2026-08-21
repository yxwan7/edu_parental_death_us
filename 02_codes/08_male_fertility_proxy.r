# ==============================================================================
# 08_male_fertility_proxy.r
#
# Build the HFC-anchored male fertility proxy and its S1/S2 sensitivity panels
# from father-age-standardized education TFRs.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/08_male_fertility_proxy.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Shared age definitions ----

# The reproductive age spans and the father-age band width are shared with the
# support gate in 13, so they are defined once in 00_functions.r.

# ---- Fertility proxy files ----

rate_dir <- here("01_data", "02_processed", "rate")
dir.create(rate_dir, recursive = TRUE, showWarnings = FALSE)

hfc_input_path <- here("01_data", "01_raw", "external", "hfc", "m_ASFR_USA.txt")
female_panel_path <- file.path(rate_dir, "03_education_fertility_rate_panel.rds")
birth_input_path <- here("01_data", "02_processed", "births", "01_birth_main_panel.rds")
pop_input_path <- here("01_data", "02_processed", "population", "01_pop_main_panel.rds")
proxy_rds_path <- file.path(rate_dir, "04_education_male_fertility_proxy_panel.rds")
s1_female_borrow_path <- file.path(
  rate_dir, "04a_education_male_fertility_proxy_panel_s1_female_borrow.rds"
)
s1_zero_gradient_path <- file.path(
  rate_dir, "04b_education_male_fertility_proxy_panel_s1_zero_gradient.rds"
)
s2_complete_case_path <- file.path(
  rate_dir, "04c_education_male_fertility_proxy_panel_s2_complete_case.rds"
)
s2_all_noba_path <- file.path(
  rate_dir, "04d_education_male_fertility_proxy_panel_s2_mnar_all_noba.rds"
)
s2_all_ba_path <- file.path(
  rate_dir, "04e_education_male_fertility_proxy_panel_s2_mnar_all_ba.rds"
)
father_education_allocation_path <- file.path(
  rate_dir, "04_father_education_missing_sensitivity_allocation.rds"
)

# ---- Male fertility components ----

# Keep only the repo's binary BA levels plus Unknown; anything else -> NA.
normalize_educ2 <- function(educ2) {
  educ2 <- normalize_missing_character(educ2)
  data.table::fifelse(
    educ2 %chin% c("No_BA", "BA_plus", "Unknown"), educ2, NA_character_
  )
}

build_female_components <- function(fertility_panel) {
  fertility_rate_input_col <- if ("fertility_rate_per_1000_model" %chin% names(fertility_panel)) {
    "fertility_rate_per_1000_model"
  } else {
    "fertility_rate_per_1000"
  }

  dt <- as.data.table(fertility_panel)[
    rate_definition == "mother_based" &
      year %in% YEARS_MODEL &
      age_years %in% AGES_FEMALE_REPRO &
      race_eth %chin% RACES_R5 &
      educ2 %chin% EDUC_LEVELS,
    .(
      year = as.integer(year),
      age_years = as.integer(age_years),
      race_eth,
      educ2,
      n_births = as.numeric(n_births),
      n_women_weighted = as.numeric(n_women_weighted),
      fertility_rate_per_1000_input = as.numeric(get(fertility_rate_input_col)),
      fertility_rate_per_1000_direct = if ("fertility_rate_per_1000_direct" %chin% names(fertility_panel)) {
        as.numeric(fertility_rate_per_1000_direct)
      } else {
        as.numeric(fertility_rate_per_1000)
      },
      small_cell_stabilized_flag = if ("small_cell_stabilized_flag" %chin% names(fertility_panel)) {
        as.logical(small_cell_stabilized_flag)
      } else {
        FALSE
      }
    )
  ]
  dt[, model_births := fertility_rate_per_1000_input * n_women_weighted / 1000]

  female_educ_age <- dt[
    ,
    .(
      births = sum(n_births, na.rm = TRUE),
      model_births = sum(model_births, na.rm = TRUE),
      women = sum(n_women_weighted, na.rm = TRUE)
    ),
    by = .(year, race_eth, educ2, age_years)
  ]
  female_educ_age[, asfr := fifelse(women > 0, model_births / women, 0)]
  female_educ_age[
    ,
    female_educ_tfr := sum(asfr),
    by = .(year, race_eth, educ2)
  ]
  female_educ_age[
    ,
    female_educ_std := fifelse(female_educ_tfr > 0, asfr / female_educ_tfr, 0)
  ]

  female_total_age <- dt[
    ,
    .(
      births = sum(n_births, na.rm = TRUE),
      model_births = sum(model_births, na.rm = TRUE),
      women = sum(n_women_weighted, na.rm = TRUE)
    ),
    by = .(year, race_eth, age_years)
  ]
  female_total_age[, female_total_asfr := fifelse(women > 0, model_births / women, 0)]
  female_total_age[
    ,
    female_total_tfr := sum(female_total_asfr),
    by = .(year, race_eth)
  ]
  female_total_age[
    ,
    female_total_std := fifelse(female_total_tfr > 0, female_total_asfr / female_total_tfr, 0)
  ]

  female_total_tfr_by_race <- unique(
    female_total_age[, .(year, race_eth, female_total_tfr)]
  )

  female_educ_tfr_share <- unique(
    female_educ_age[, .(year, race_eth, educ2, female_educ_tfr)]
  )[
    female_total_tfr_by_race,
    on = .(year, race_eth)
  ]
  female_educ_tfr_share[
    ,
    female_educ_tfr_share := fifelse(
      female_total_tfr > 0,
      female_educ_tfr / female_total_tfr,
      0
    )
  ]

  female_ratio_age <- merge(
    female_educ_age[
      ,
      .(year, race_eth, educ2, age_years, female_educ_std)
    ],
    female_total_age[
      ,
      .(year, race_eth, age_years, female_total_std)
    ],
    by = c("year", "race_eth", "age_years"),
    all.x = TRUE,
    sort = FALSE
  )
  female_ratio_age[
    ,
    female_ratio_raw := fifelse(
      female_total_std > 0,
      female_educ_std / female_total_std,
      1
    )
  ]

  female_national_age <- dt[
    ,
    .(
      births = sum(n_births, na.rm = TRUE),
      model_births = sum(model_births, na.rm = TRUE),
      women = sum(n_women_weighted, na.rm = TRUE)
    ),
    by = .(year, age_years)
  ]
  female_national_age[, female_national_asfr := fifelse(women > 0, model_births / women, 0)]
  female_national_tfr <- female_national_age[
    ,
    .(female_national_tfr = sum(female_national_asfr)),
    by = .(year)
  ]

  list(
    female_ratio_age = female_ratio_age,
    female_ratio_modeled = build_smoothed_female_ratio(female_ratio_age),
    female_total_tfr_by_race = female_total_tfr_by_race,
    female_educ_tfr_share = female_educ_tfr_share,
    female_national_tfr = female_national_tfr,
    fertility_rate_input_col = fertility_rate_input_col,
    fertility_input_stabilized_rows = dt[small_cell_stabilized_flag == TRUE, .N],
    fertility_input_model_birth_delta = dt[
      ,
      sum(
        (fertility_rate_per_1000_input - fertility_rate_per_1000_direct) *
          n_women_weighted / 1000,
        na.rm = TRUE
      )
    ]
  )
}

build_smoothed_female_ratio <- function(female_ratio_age) {
  model_input <- copy(female_ratio_age)[
    is.finite(female_ratio_raw) &
      age_years %in% AGES_FEMALE_REPRO,
    .(
      year,
      race_eth,
      educ2,
      ratio_age_source = age_years,
      ratio = as.numeric(female_ratio_raw)
    )
  ]

  modeled <- model_input[
    ,
    {
      dt <- .SD[order(year)]
      if (.N >= 2L && length(unique(dt$ratio)) > 1L) {
        pred <- as.numeric(
          predict(
            stats::lm(ratio ~ year, data = dt),
            newdata = data.frame(year = YEARS_MODEL)
          )
        )
      } else if (.N >= 1L) {
        pred <- rep(dt$ratio[.N], length(YEARS_MODEL))
      } else {
        pred <- rep(1, length(YEARS_MODEL))
      }

      data.table(year = YEARS_MODEL, female_ratio_modeled = pmax(pred, 0))
    },
    by = .(race_eth, educ2, ratio_age_source)
  ]

  out <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R5,
    educ2 = EDUC_LEVELS,
    age_years = AGES_MALE_REPRO,
    unique = TRUE
  )
  out[
    ,
    ratio_age_source := pmin(age_years, max(AGES_FEMALE_REPRO))
  ]
  out <- modeled[
    out,
    on = .(year, race_eth, educ2, ratio_age_source)
  ]
  out[
    is.na(female_ratio_modeled) | !is.finite(female_ratio_modeled),
    female_ratio_modeled := 1
  ]
  out[order(year, race_eth, educ2, age_years)]
}

build_hfc_male_total <- function(hfc_path, female_national_tfr) {
  hfc_dt <- fread(hfc_path)
  setnames(hfc_dt, tolower(names(hfc_dt)))
  male_direct <- hfc_dt[
    year %in% 2000:2015 & age %in% AGES_MALE_REPRO,
    .(
      year = as.integer(year),
      age_years = as.integer(age),
      male_total_asfr = as.numeric(asfr)
    )
  ][order(year, age_years)]

  male_direct_tfr <- male_direct[
    ,
    .(male_total_tfr = sum(male_total_asfr)),
    by = .(year)
  ]

  female_national_2015 <- female_national_tfr[year %in% 2000:2015]
  ratio_obs <- merge(
    male_direct_tfr,
    female_national_2015,
    by = "year",
    all.x = TRUE,
    sort = TRUE
  )
  ratio_obs[
    ,
    national_male_female_tfr_ratio := fifelse(
      female_national_tfr > 0,
      male_total_tfr / female_national_tfr,
      NA_real_
    )
  ]

  ratio_fit <- stats::lm(
    national_male_female_tfr_ratio ~ year,
    data = ratio_obs[year %in% 2010:2015]
  )
  ratio_forecast <- data.table(
    year = 2016:2024,
    national_male_female_tfr_ratio = as.numeric(
      predict(ratio_fit, newdata = data.frame(year = 2016:2024))
    )
  )
  ratio_forecast[
    ,
    national_male_female_tfr_ratio := pmax(national_male_female_tfr_ratio, 0.01)
  ]

  male_forecast_raw <- male_direct[
    year %in% 2010:2015 & male_total_asfr > 0,
    {
      dt <- .SD[order(year)]
      if (.N >= 2L) {
        pred <- exp(
          predict(
            stats::lm(log(male_total_asfr) ~ year, data = dt),
            newdata = data.frame(year = 2016:2024)
          )
        )
      } else if (.N == 1L) {
        pred <- rep(dt$male_total_asfr[1L], length(2016:2024))
      } else {
        pred <- rep(0, length(2016:2024))
      }

      data.table(year = 2016:2024, male_total_asfr_raw = as.numeric(pred))
    },
    by = .(age_years)
  ]

  male_forecast_grid <- CJ(
    year = 2016:2024,
    age_years = AGES_MALE_REPRO,
    unique = TRUE
  )
  male_forecast_raw <- male_forecast_raw[
    male_forecast_grid,
    on = .(year, age_years)
  ]
  male_forecast_raw[
    is.na(male_total_asfr_raw) | !is.finite(male_total_asfr_raw),
    male_total_asfr_raw := 0
  ]
  male_forecast_raw[
    ,
    male_total_asfr_std := if (sum(male_total_asfr_raw) > 0) {
      male_total_asfr_raw / sum(male_total_asfr_raw)
    } else {
      rep(0, .N)
    },
    by = .(year)
  ]

  extended_total <- merge(
    female_national_tfr[year %in% 2016:2024],
    ratio_forecast,
    by = "year",
    all.x = TRUE,
    sort = TRUE
  )[
    ,
    .(
      year,
      male_total_tfr = female_national_tfr * national_male_female_tfr_ratio
    )
  ]

  male_extended <- merge(
    male_forecast_raw,
    extended_total,
    by = "year",
    all.x = TRUE,
    sort = FALSE
  )[
    ,
    .(
      year,
      age_years,
      male_total_asfr = male_total_tfr * male_total_asfr_std
    )
  ]

  ratio_all <- rbindlist(
    list(
      ratio_obs[, .(year, national_male_female_tfr_ratio)],
      ratio_forecast
    ),
    use.names = TRUE
  )[order(year)]

  male_total_all <- rbindlist(
    list(
      male_direct[, .(year, age_years, male_total_asfr)],
      male_extended
    ),
    use.names = TRUE
  )[order(year, age_years)]

  list(
    male_total_all = male_total_all,
    ratio_all = ratio_all
  )
}

allocate_missing_father_education <- function(
    reported,
    missing_by_mother,
    missing_ba_odds_ratio) {
  cell_keys <- c("year", "race_eth", "mother_educ")
  reported_by_cell <- reported[
    ,
    .(n_reported = sum(n_births)),
    by = c(cell_keys, "father_educ")
  ]
  reported_by_cell[
    ,
    reported_share := n_reported / sum(n_reported),
    by = cell_keys
  ]
  reported_totals <- reported_by_cell[
    ,
    .(
      n_reported_total = sum(n_reported),
      p_ba_reported = sum(n_reported[father_educ == "BA_plus"]) /
        sum(n_reported)
    ),
    by = cell_keys
  ]

  allocation_base <- merge(
    reported_totals,
    missing_by_mother,
    by = cell_keys,
    all.x = TRUE,
    sort = FALSE
  )
  allocation_base[is.na(n_missing), n_missing := 0]
  allocation_base[
    ,
    p_ba_missing := if (missing_ba_odds_ratio == 0) {
      rep(0, .N)
    } else if (is.infinite(missing_ba_odds_ratio)) {
      rep(1, .N)
    } else {
      plogis(qlogis(p_ba_reported) + log(missing_ba_odds_ratio))
    }
  ]

  allocation <- rbindlist(
    list(
      allocation_base[
        ,
        .(
          year,
          race_eth,
          mother_educ,
          father_educ = "BA_plus",
          n_reported_total,
          n_missing,
          p_ba_reported,
          p_ba_missing
        )
      ],
      allocation_base[
        ,
        .(
          year,
          race_eth,
          mother_educ,
          father_educ = "No_BA",
          n_reported_total,
          n_missing,
          p_ba_reported,
          p_ba_missing
        )
      ]
    )
  )
  allocation <- merge(
    allocation,
    reported_by_cell,
    by = c(cell_keys, "father_educ"),
    all.x = TRUE,
    sort = FALSE
  )
  allocation[is.na(n_reported), n_reported := 0]
  allocation[is.na(reported_share), reported_share := 0]
  if (missing_ba_odds_ratio == 1) {
    allocation[, missing_share := reported_share]
  } else {
    allocation[
      ,
      missing_share := fifelse(
        father_educ == "BA_plus",
        p_ba_missing,
        1 - p_ba_missing
      )
    ]
  }
  allocation[
    ,
    `:=`(
      n_imputed = n_missing * missing_share,
      n_completed = n_reported + n_missing * missing_share
    )
  ]

  imputed <- allocation[
    ,
    .(n_imputed = sum(n_imputed)),
    by = .(year, race_eth, father_educ)
  ]
  completed_births <- merge(
    reported[
      ,
      .(n_reported = sum(n_births)),
      by = .(year, race_eth, father_educ)
    ],
    imputed,
    by = c("year", "race_eth", "father_educ"),
    all = TRUE
  )
  completed_births[is.na(n_reported), n_reported := 0]
  completed_births[is.na(n_imputed), n_imputed := 0]
  completed_births[, n_births := n_reported + n_imputed]

  list(
    completed_births = completed_births,
    allocation = allocation[
      order(year, race_eth, mother_educ, father_educ)
    ]
  )
}

build_male_educ_share <- function(
    birth_main_panel,
    female_share,
    pop_path,
    father_observed_years = 2010:2024,
    runup_years = 2000:2009,
    anchor_years = 2010:2012,
    impute_father_educ = TRUE,
    missing_ba_odds_ratio = 1) {
  # ACS male exposure is grouped in five-year age bands before calculating
  # age-standardized education-specific fertility levels.
  population <- as.data.table(readRDS(pop_path))
  population[, educ2 := normalize_educ2(educ2)]
  male_exposure <- population[
    sex == "male" &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS &
      age_years %in% AGES_MALE_REPRO &
      year %in% c(runup_years, father_observed_years),
    .(n_men = sum(as.numeric(n_persons_weighted))),
    by = .(
      year,
      race_eth,
      educ2,
      ageg = FATHER_AGE_BAND_WIDTH * (age_years %/% FATHER_AGE_BAND_WIDTH)
    )
  ]

  # Use births from years with sufficiently complete father education.
  father_births <- as.data.table(birth_main_panel)[
    resident_birth_flag %in% TRUE &
      year %in% father_observed_years &
      father_race_eth %chin% RACES_R3,
    .(
      year = as.integer(year),
      race_eth = father_race_eth,
      father_educ = normalize_educ2(father_educ2),
      mother_educ = normalize_educ2(mother_educ2),
      father_age_years = as.integer(father_age_years),
      n_births = as.numeric(n_births)
    )
  ]
  father_births[!father_educ %chin% EDUC_LEVELS, father_educ := NA_character_]
  father_births[!mother_educ %chin% EDUC_LEVELS, mother_educ := "Unknown"]

  reported <- father_births[!is.na(father_educ)]
  if (impute_father_educ) {
    missing_by_mother <- father_births[
      is.na(father_educ),
      .(n_missing = sum(n_births)), by = .(year, race_eth, mother_educ)
    ]
    father_education_result <- allocate_missing_father_education(
      reported = reported,
      missing_by_mother = missing_by_mother,
      missing_ba_odds_ratio = missing_ba_odds_ratio
    )
    # Preserve the frozen MAR summation order.
    if (missing_ba_odds_ratio == 1) {
      father_educ_share_by_mother <- reported[
        ,
        .(n_births = sum(n_births)),
        by = .(year, race_eth, mother_educ, father_educ)
      ][
        ,
        share := n_births / sum(n_births),
        by = .(year, race_eth, mother_educ)
      ]
      imputed <- missing_by_mother[
        father_educ_share_by_mother,
        on = .(year, race_eth, mother_educ),
        allow.cartesian = TRUE,
        nomatch = NULL
      ][
        ,
        .(n_imputed = sum(n_missing * share)),
        by = .(year, race_eth, father_educ)
      ]
      completed_births <- merge(
        reported[
          ,
          .(n_reported = sum(n_births)),
          by = .(year, race_eth, father_educ)
        ],
        imputed,
        by = c("year", "race_eth", "father_educ"),
        all = TRUE
      )
      completed_births[is.na(n_reported), n_reported := 0]
      completed_births[is.na(n_imputed), n_imputed := 0]
      completed_births[, n_births := n_reported + n_imputed]
    } else {
      completed_births <- father_education_result$completed_births
    }
    setnames(completed_births, "father_educ", "educ2")
  } else {
    completed_births <- reported[
      , .(n_births = sum(n_births)), by = .(year, race_eth, educ2 = father_educ)
    ]
    father_education_result <- NULL
  }

  # Age structure of each education group's fertility, in 5-year father-age bands,
  # taken from births with the father's age reported (15-59). Father age is missing
  # on the same records as father education, so this assumes the reported-age
  # structure represents the group (MAR on father age within year x race x
  # education); it is weakest for Hispanic, where father reporting is least
  # complete. This is the single simplifying assumption of the age standardization.
  father_age_share_year <- reported[
    !is.na(father_age_years) & father_age_years %in% AGES_MALE_REPRO,
    .(n_births = sum(n_births)),
    by = .(
      year,
      race_eth,
      educ2 = father_educ,
      ageg = FATHER_AGE_BAND_WIDTH *
        (father_age_years %/% FATHER_AGE_BAND_WIDTH)
    )
  ][
    , age_share := n_births / sum(n_births), by = .(year, race_eth, educ2)
  ]
  # Father age is 100% missing in 2012-2013 (a natality-file artifact), and some
  # otherwise observed cells omit an isolated age band. Fall back at the band
  # level to the (race, education) structure pooled over all father-observed years,
  # then renormalize each year-race-education cell to sum to 1.
  father_age_share_pooled <- reported[
    !is.na(father_age_years) & father_age_years %in% AGES_MALE_REPRO,
    .(n_births = sum(n_births)),
    by = .(
      race_eth,
      educ2 = father_educ,
      ageg = FATHER_AGE_BAND_WIDTH *
        (father_age_years %/% FATHER_AGE_BAND_WIDTH)
    )
  ][
    , pooled_share := n_births / sum(n_births), by = .(race_eth, educ2)
  ]
  father_age_share <- father_age_share_pooled[
    , .(race_eth, educ2, ageg, pooled_share)
  ][
    unique(completed_births[, .(year, race_eth, educ2)]),
    on = .(race_eth, educ2), allow.cartesian = TRUE
  ]
  father_age_share[
    father_age_share_year,
    on = .(year, race_eth, educ2, ageg),
    year_share := i.age_share
  ]
  father_age_share[, age_share := fifelse(is.na(year_share), pooled_share, year_share)]
  father_age_share[, age_share := age_share / sum(age_share), by = .(year, race_eth, educ2)]

  # Each education group's fertility as its TFR relative to the race pooled TFR:
  # spread births over father-age bands, divide by age-matched male exposure for
  # band ASFRs, and sum. Both TFRs are age-composition invariant.
  observed_age <- merge(
    completed_births[, .(year, race_eth, educ2, n_births_educ = n_births)],
    father_age_share[, .(year, race_eth, educ2, ageg, age_share)],
    by = c("year", "race_eth", "educ2"),
    allow.cartesian = TRUE
  )
  observed_age[, n_births := n_births_educ * age_share]
  observed_age <- merge(
    observed_age[, .(year, race_eth, educ2, ageg, n_births)],
    male_exposure,
    by = c("year", "race_eth", "educ2", "ageg")
  )
  observed_age[, asfr := fifelse(n_men > 0, n_births / n_men, 0)]

  # The 5-year band width cancels in the education/pooled ratio;
  # kept explicit so educ_tfr and pooled_tfr are genuine TFRs (for naming and any
  # absolute checks).
  educ_tfr <- observed_age[
    , .(educ_tfr = FATHER_AGE_BAND_WIDTH * sum(asfr)),
    by = .(year, race_eth, educ2)
  ]
  pooled_tfr <- observed_age[
    , .(n_births = sum(n_births), n_men = sum(n_men)), by = .(year, race_eth, ageg)
  ][
    ,
    .(
      pooled_tfr = FATHER_AGE_BAND_WIDTH * sum(
        fifelse(n_men > 0, n_births / n_men, 0)
      )
    ),
    by = .(year, race_eth)
  ]
  observed_level <- merge(educ_tfr, pooled_tfr, by = c("year", "race_eth"))
  observed_level[
    ,
    male_educ_rel := fifelse(pooled_tfr > 0, educ_tfr / pooled_tfr, 1)
  ]

  female_level <- copy(female_share)[
    race_eth %chin% RACES_R3,
    .(year, race_eth, educ2, female_educ_rel = female_educ_tfr_share)
  ]

  # The run-up (pre-2010) has almost no reported father education, so the male
  # level cannot be observed there. Carry it back from the earliest observed
  # anchor years, scaled by how each run-up year's observed female relative level
  # departs from its own anchor. This preserves the year-specific female trend and
  # assumes only that the male-to-female level ratio is stable across the bridge;
  # it is continuous at the father-observed boundary. Relatives are centered on 1,
  # so no renormalization is needed.
  male_anchor <- observed_level[
    year %in% anchor_years,
    .(male_anchor = mean(male_educ_rel)), by = .(race_eth, educ2)
  ]
  female_anchor <- female_level[
    year %in% anchor_years,
    .(female_anchor = mean(female_educ_rel)), by = .(race_eth, educ2)
  ]
  runup_level <- merge(
    female_level[year %in% runup_years], male_anchor,
    by = c("race_eth", "educ2")
  )
  runup_level <- merge(runup_level, female_anchor, by = c("race_eth", "educ2"))
  runup_level[
    ,
    male_educ_rel := male_anchor * fifelse(
      female_anchor > 0, female_educ_rel / female_anchor, 1
    )
  ]

  male_educ_level <- rbind(
    observed_level[, .(year, race_eth, educ2, male_educ_tfr_share = male_educ_rel)],
    runup_level[, .(year, race_eth, educ2, male_educ_tfr_share = male_educ_rel)]
  )

  list(
    male_educ_level = male_educ_level[order(year, race_eth, educ2)],
    father_education_allocation = if (is.null(father_education_result)) {
      NULL
    } else {
      father_education_result$allocation
    }
  )
}

# Assemble a full male proxy panel from a given education-level share. Everything
# except educ_share is shared across legs, so swapping educ_share isolates the
# education-gradient assumption.
assemble_male_proxy <- function(female_components, hfc_components, educ_share,
                                proxy_definition_tag) {
  male_proxy_age <- merge(
    female_components$female_ratio_modeled[
      race_eth %chin% RACES_R3,
      .(
        year,
        race_eth,
        educ2,
        age_years,
        ratio_age_source,
        female_ratio_modeled
      )
    ],
    hfc_components$male_total_all,
    by = c("year", "age_years"),
    all.x = TRUE,
    sort = FALSE
  )
  male_proxy_age[
    ,
    male_proxy_raw := female_ratio_modeled * male_total_asfr
  ]
  male_proxy_age[
    ,
    male_proxy_std := if (sum(male_proxy_raw) > 0) {
      male_proxy_raw / sum(male_proxy_raw)
    } else {
      rep(0, .N)
    },
    by = .(year, race_eth, educ2)
  ]

  male_proxy_age <- merge(
    male_proxy_age,
    female_components$female_total_tfr_by_race[
      race_eth %chin% RACES_R3
    ],
    by = c("year", "race_eth"),
    all.x = TRUE,
    sort = FALSE
  )
  male_proxy_age <- merge(
    male_proxy_age,
    female_components$female_educ_tfr_share[
      race_eth %chin% RACES_R3,
      .(year, race_eth, educ2, female_educ_tfr_share)
    ],
    by = c("year", "race_eth", "educ2"),
    all.x = TRUE,
    sort = FALSE
  )
  male_proxy_age <- merge(
    male_proxy_age,
    educ_share,
    by = c("year", "race_eth", "educ2"),
    all.x = TRUE,
    sort = FALSE
  )
  male_proxy_age <- merge(
    male_proxy_age,
    hfc_components$ratio_all,
    by = "year",
    all.x = TRUE,
    sort = FALSE
  )
  male_proxy_age[
    ,
    `:=`(
      male_total_tfr_race = female_total_tfr * national_male_female_tfr_ratio,
      male_educ_tfr = female_total_tfr * national_male_female_tfr_ratio * male_educ_tfr_share
    )
  ]
  male_proxy_age[
    ,
    asfr_male_proxy := male_proxy_std * male_educ_tfr
  ]

  proxy_full <- CJ(
    year = YEARS_MODEL,
    age_years = AGES_MODEL,
    race_eth = RACES_R3,
    educ2 = EDUC_LEVELS,
    unique = TRUE
  )
  proxy_full <- merge(
    proxy_full,
    male_proxy_age[
      ,
      .(
        year,
        age_years,
        race_eth,
        educ2,
        asfr_male_proxy,
        female_ratio_modeled,
        ratio_age_source,
        female_total_tfr,
        female_educ_tfr_share,
        male_educ_tfr_share,
        national_male_female_tfr_ratio,
        male_total_tfr_race,
        male_educ_tfr
      )
    ],
    by = c("year", "age_years", "race_eth", "educ2"),
    all.x = TRUE,
    sort = TRUE
  )
  proxy_full[
    ,
    asfr_male_proxy := fifelse(is.na(asfr_male_proxy), 0, asfr_male_proxy)
  ]
  proxy_full[
    ,
    `:=`(
      rate_definition = "male_hfc_proxy",
      proxy_definition = proxy_definition_tag,
      fertility_rate_per_1000 = asfr_male_proxy * 1000,
      fx = asfr_male_proxy
    )
  ]

  setcolorder(
    proxy_full,
    c(
      "rate_definition",
      "proxy_definition",
      "year",
      "age_years",
      "race_eth",
      "educ2",
      "fertility_rate_per_1000",
      "fx",
      "female_ratio_modeled",
      "ratio_age_source",
      "female_total_tfr",
      "female_educ_tfr_share",
      "male_educ_tfr_share",
      "national_male_female_tfr_ratio",
      "male_total_tfr_race",
      "male_educ_tfr"
    )
  )
  proxy_full[]
}

# ---- Build HFC-anchored male proxy ----

# The proxy anchors the age shape on the HFC national male ASFR and sets the
# education level from father-observed births; the female fertility panel (07)
# still supplies the education age offset, the race TFR level, and the run-up
# bridge trend.

female_panel <- readRDS(female_panel_path)
birth_main_panel <- readRDS(birth_input_path)

female_components <- build_female_components(female_panel)
hfc_components <- build_hfc_male_total(
  hfc_path = hfc_input_path,
  female_national_tfr = female_components$female_national_tfr
)

primary_result <- build_male_educ_share(
  birth_main_panel = birth_main_panel,
  female_share = female_components$female_educ_tfr_share,
  pop_path = pop_input_path
)
male_educ_level <- primary_result$male_educ_level

run_sensitivity_outputs <-
  !nzchar(Sys.getenv("EDU_PARENT_LOSS_DRAW_ROOT", unset = ""))

if (run_sensitivity_outputs) {
  # S1 alternatives change only the education-level share.
  female_borrow_share <- female_components$female_educ_tfr_share[
    race_eth %chin% RACES_R3,
    .(year, race_eth, educ2, male_educ_tfr_share = female_educ_tfr_share)
  ]

  zero_gradient_share <- CJ(
    year = YEARS_MODEL,
    race_eth = RACES_R3,
    educ2 = EDUC_LEVELS,
    unique = TRUE
  )[, male_educ_tfr_share := 1]

  complete_case_result <- build_male_educ_share(
    birth_main_panel = birth_main_panel,
    female_share = female_components$female_educ_tfr_share,
    pop_path = pop_input_path,
    impute_father_educ = FALSE
  )
  complete_case_share <- complete_case_result$male_educ_level

  mnar_odds_ratio <- setNames(
    S2_MNAR_SCENARIOS$missing_ba_odds_ratio,
    S2_MNAR_SCENARIOS$leg
  )
  all_noba_result <- build_male_educ_share(
    birth_main_panel = birth_main_panel,
    female_share = female_components$female_educ_tfr_share,
    pop_path = pop_input_path,
    missing_ba_odds_ratio = mnar_odds_ratio[["s2_mnar_all_noba"]]
  )
  all_ba_result <- build_male_educ_share(
    birth_main_panel = birth_main_panel,
    female_share = female_components$female_educ_tfr_share,
    pop_path = pop_input_path,
    missing_ba_odds_ratio = mnar_odds_ratio[["s2_mnar_all_ba"]]
  )
}

proxy_definition_primary <-
  "hfc_male_age_shape_x_father_observed_imputed_educ_share_2010_2024_anchor_trend_bridge_prior"

proxy_full <- assemble_male_proxy(
  female_components,
  hfc_components,
  male_educ_level,
  proxy_definition_primary
)

published_proxy <- proxy_full[
  year %in% 2000:2024 &
    race_eth %chin% RACES_R3
]

# ---- Primary proxy panel ----

saveRDS(published_proxy, proxy_rds_path)

if (run_sensitivity_outputs) {
  # ---- S1 and S2 proxy panels ----

  sensitivity_variants <- list(
    list(
      share = female_borrow_share,
      tag = "hfc_male_age_shape_x_female_borrowed_educ_share_s1_upper",
      path = s1_female_borrow_path
    ),
    list(
      share = zero_gradient_share,
      tag = "hfc_male_age_shape_x_education_neutral_educ_share_s1_zero",
      path = s1_zero_gradient_path
    ),
    list(
      share = complete_case_share,
      tag = "hfc_male_age_shape_x_father_observed_reported_educ_share_complete_case",
      path = s2_complete_case_path
    ),
    list(
      share = all_noba_result$male_educ_level,
      tag = "hfc_male_age_shape_x_father_educ_share_s2_mnar_all_noba",
      path = s2_all_noba_path
    ),
    list(
      share = all_ba_result$male_educ_level,
      tag = "hfc_male_age_shape_x_father_educ_share_s2_mnar_all_ba",
      path = s2_all_ba_path
    )
  )

  for (variant in sensitivity_variants) {
    panel <- assemble_male_proxy(
      female_components,
      hfc_components,
      variant$share,
      variant$tag
    )[
      year %in% 2000:2024 & race_eth %chin% RACES_R3
    ]
    saveRDS(panel, variant$path)
  }

  father_education_allocation <- rbindlist(
    list(
      primary = primary_result$father_education_allocation,
      s2_mnar_all_noba = all_noba_result$father_education_allocation,
      s2_mnar_all_ba = all_ba_result$father_education_allocation
    ),
    idcol = "leg"
  )
  father_education_allocation[
    ,
    missing_ba_odds_ratio := fcase(
      leg == "primary", 1,
      leg == "s2_mnar_all_noba", 0,
      leg == "s2_mnar_all_ba", Inf
    )
  ]
  setcolorder(
    father_education_allocation,
    c(
      "leg",
      "missing_ba_odds_ratio",
      setdiff(
        names(father_education_allocation),
        c("leg", "missing_ba_odds_ratio")
      )
    )
  )
  saveRDS(
    father_education_allocation,
    father_education_allocation_path,
    compress = "xz"
  )
}
