# ==============================================================================
# 06_mortality_rates_by_cause.r
#
# Split all-cause mortality into the 24 atom causes, carrying the 6 causes
# derived from them.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/06_mortality_rates_by_cause.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Cause-rate files ----

rate_dir <- here("01_data", "02_processed", "rate")
dir.create(rate_dir, recursive = TRUE, showWarnings = FALSE)

death_input_path <- here("01_data", "02_processed", "deaths", "01_death_main_panel.rds")
all_cause_input_path <- file.path(rate_dir, "01_education_mortality_rate_panel.rds")

# The cause rate panel: the cause_atom split, with cause_6 derived from it.
cause_atom_output_rds_path <- file.path(
  rate_dir, "02_education_mortality_rate_panel_by_cause.rds"
)

# ---- Cause allocation functions ----

# Cause shares within each donor cell, one column per cause level. Cells with no
# deaths at all carry no share and are left for a coarser donor.
cause_share_table <- function(cause_counts, by_cols, cause_col, cause_levels,
                              share_cols, donor_source) {
  count_dt <- cause_counts[
    ,
    .(n_deaths_raw = sum(n_deaths_raw, na.rm = TRUE)),
    by = c(by_cols, cause_col)
  ]

  share_dt <- dcast(
    count_dt,
    as.formula(paste(paste(by_cols, collapse = " + "), "~", cause_col)),
    value.var = "n_deaths_raw",
    fill = 0
  )

  for (cause_level in cause_levels) {
    if (!cause_level %chin% names(share_dt)) {
      share_dt[, (cause_level) := 0]
    }
  }

  share_dt[, raw_total := rowSums(.SD), .SDcols = cause_levels]
  share_dt <- share_dt[raw_total > 0]
  share_dt[
    ,
    (share_cols) := lapply(.SD, function(x) x / raw_total),
    .SDcols = cause_levels
  ]
  share_dt[, cause_share_source := donor_source]

  share_dt[, c(by_cols, share_cols, "cause_share_source"), with = FALSE]
}

# Fill the cells still without a share from a coarser donor table.
fill_missing_shares <- function(target_dt, share_dt, by_cols, share_cols) {
  missing_idx <- which(is.na(target_dt[[share_cols[[1]]]]))
  if (!length(missing_idx)) {
    return(target_dt)
  }

  joined <- share_dt[target_dt[missing_idx], on = by_cols]
  matched <- !is.na(joined[[share_cols[[1]]]])
  if (!any(matched)) {
    return(target_dt)
  }

  row_idx <- missing_idx[matched]
  for (share_col in share_cols) {
    target_dt[row_idx, (share_col) := joined[[share_col]][matched]]
  }
  target_dt[row_idx, cause_share_source := joined$cause_share_source[matched]]

  target_dt[]
}

# Donor ladder: try the finest cell first, then progressively coarser ones,
# until every completed cell carries a set of cause shares.
CAUSE_SHARE_DONORS <- list(
  list(by_cols = c("year", "sex", "age_years", "race_eth", "educ2"),
       source = "direct_year_sex_age_race_educ"),
  list(by_cols = c("year", "sex", "age_years", "race_eth"),
       source = "donor_year_sex_age_race"),
  list(by_cols = c("year", "sex", "age_years", "educ2"),
       source = "donor_year_sex_age_educ"),
  list(by_cols = c("year", "sex", "age_years"),
       source = "donor_year_sex_age"),
  list(by_cols = c("year", "sex", "race_eth", "educ2"),
       source = "donor_year_sex_race_educ_all_age"),
  list(by_cols = c("year", "sex", "race_eth"),
       source = "donor_year_sex_race_all_age_educ"),
  list(by_cols = c("year", "sex"),
       source = "donor_year_sex_all_age_race_educ"),
  list(by_cols = c("sex", "age_years"),
       source = "donor_all_years_sex_age"),
  list(by_cols = c("sex"),
       source = "donor_all_years_sex")
)

assign_cause_shares <- function(completed_dt, cause_counts, cause_col, cause_levels,
                                share_cols) {
  target_dt <- copy(completed_dt)

  for (donor in CAUSE_SHARE_DONORS) {
    share_dt <- cause_share_table(
      cause_counts = cause_counts,
      by_cols = donor$by_cols,
      cause_col = cause_col,
      cause_levels = cause_levels,
      share_cols = share_cols,
      donor_source = donor$source
    )
    target_dt <- fill_missing_shares(target_dt, share_dt, donor$by_cols, share_cols)
  }

  target_dt[]
}

build_direct_cause_atom_counts <- function(death_main_panel) {
  dt <- as.data.table(death_main_panel)[
    year %in% YEARS_MODEL &
      sex %chin% SEX_LEVELS &
      !is.na(age_years),
    .(
      year = as.integer(year),
      sex = normalize_missing_character(sex),
      age_years = as.integer(age_years),
      race_eth = normalize_missing_character(race_eth),
      educ2 = normalize_missing_character(educ2),
      cause_atom = normalize_missing_character(cause_atom),
      n_deaths = as.numeric(n_deaths)
    )
  ]

  dt <- bridge_child_unknown_educ2(dt)
  dt <- dt[
    race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS
  ]

  dt[
    ,
    .(n_deaths_raw = sum(n_deaths, na.rm = TRUE)),
    by = .(year, sex, age_years, race_eth, educ2, cause_atom)
  ]
}

# The all-cause panel with empty share columns, one per cause level of whichever
# family is being split, plus a row id used to restore the original order.
prepare_all_cause_panel <- function(all_cause_panel, share_cols) {
  dt <- copy(as.data.table(all_cause_panel))[
    rate_definition == "all_cause" &
      year %in% YEARS_MODEL &
      sex %chin% SEX_LEVELS &
      race_eth %chin% RACES_R3 &
      educ2 %chin% EDUC_LEVELS
  ]

  dt[
    ,
    `:=`(
      year = as.integer(year),
      age_years = as.integer(age_years),
      n_deaths = as.numeric(n_deaths),
      n_pop_unweighted = as.numeric(n_pop_unweighted),
      n_pop_weighted = as.numeric(n_pop_weighted)
    )
  ]

  dt[, (share_cols) := NA_real_]
  dt[, cause_share_source := NA_character_]
  dt[, completed_cell_id := .I]

  dt[]
}

# The cause_atom split, with cause_6 derived from it. Every level is taken
# directly from its share.
split_completed_deaths_cause_atom <- function(completed_dt, cause_atom_counts) {
  share_cols <- paste0("share__", CAUSE_ATOM_LEVELS)
  target_dt <- assign_cause_shares(
    completed_dt = completed_dt,
    cause_counts = cause_atom_counts,
    cause_col = "cause_atom",
    cause_levels = CAUSE_ATOM_LEVELS,
    share_cols = share_cols
  )

  cause_atom_panel <- rbindlist(
    lapply(
      CAUSE_ATOM_LEVELS,
      function(cause_level) {
        share_col <- paste0("share__", cause_level)
        cause_dt <- copy(target_dt)
        cause_dt[
          ,
          `:=`(
            cause_atom = cause_level,
            cause_6 = classify_cause_6_from_cause_atom(cause_level),
            cause_share = get(share_col),
            n_deaths_cause = n_deaths * get(share_col)
          )
        ]
        cause_dt
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  cause_atom_panel[
    ,
    `:=`(
      n_deaths_all_cause = n_deaths,
      n_deaths = n_deaths_cause,
      rate_definition = "cause_atom",
      mortality_rate_per_100k = 1e5 * n_deaths_cause / n_pop_weighted
    )
  ]

  cause_atom_panel[
    ,
    c("completed_cell_id", share_cols, "n_deaths_cause") := NULL
  ]

  setcolorder(
    cause_atom_panel,
    c(
      "rate_definition", "cause_atom", "cause_6",
      "year", "sex", "age_years", "race_eth", "educ2",
      "n_deaths", "n_deaths_all_cause", "cause_share",
      "cause_share_source", "n_pop_unweighted", "n_pop_weighted",
      "mortality_rate_per_100k"
    )
  )
  setorderv(
    cause_atom_panel,
    c("year", "sex", "age_years", "race_eth", "educ2", "cause_atom")
  )

  cause_atom_panel[]
}

# ---- Split the all-cause panel by cause ----

all_cause_panel <- readRDS(all_cause_input_path)
completed_mortality_cause_atom <- prepare_all_cause_panel(
  all_cause_panel, paste0("share__", CAUSE_ATOM_LEVELS)
)
rm(all_cause_panel)

death_main_panel <- readRDS(death_input_path)
cause_atom_counts <- build_direct_cause_atom_counts(death_main_panel)
rm(death_main_panel)
gc()

cause_atom_panel <- split_completed_deaths_cause_atom(completed_mortality_cause_atom, cause_atom_counts)

# ---- Write the cause rate panel ----

saveRDS(cause_atom_panel, cause_atom_output_rds_path)
