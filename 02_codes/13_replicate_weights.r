# ==============================================================================
# 13_replicate_weights.r
#
# Build the analysis-cell population panel with the ACS point weight and 80
# replicate weights.
# ==============================================================================

# ---- Packages ----

# tidyselect is named explicitly because the variable selection below calls
# all_of(): relying on ipumsr to re-export it makes this script depend on another
# package's export list.
suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ipumsr)
})

here::i_am("02_codes/13_replicate_weights.r")

# Use the same ACS recodes as the point-weight panel.
source(here::here("02_codes", "00_functions.r"))

# ---- Replicate-panel definitions ----

# The extract is a 54 GB fixed-width file read in chunks; only the chunk size is
# tunable from the shell. The year window is fixed because a replicate panel
# covering less than the full window cannot serve the draws.
chunk_size <- {
  raw_value <- Sys.getenv("IPUMS_POP_CHUNK_SIZE", unset = "100000")
  parsed <- suppressWarnings(as.integer(raw_value))
  if (is.na(parsed) || parsed <= 0L) 100000L else parsed
}

# Mirrors the ACS window in 04; a divergence would surface immediately as cells
# present in one panel and absent from the other in the equivalence gate.
ACS_MIN_YEAR <- 2000L
ACS_MAX_YEAR <- 2024L
REQUESTED_YEARS <- ACS_MIN_YEAR:ACS_MAX_YEAR

EDUC_VAR <- "educd"
WEIGHT_VAR <- "perwt"

# REPWTP indicates availability; REPWTP1-80 contain the replicate weights.
REPLICATE_VARS <- sprintf("repwtp%d", seq_len(ACS_REPLICATE_COUNT))
REPLICATE_COLUMNS <- acs_replicate_columns()
WEIGHT_COLUMNS <- c("n_persons_weighted", REPLICATE_COLUMNS)

# The panel is reduced onto this fixed key, so the accumulator stays bounded at
# the number of analysis cells (about 90k rows) no matter how long the file is.
ANALYSIS_KEY <- c("source", "year", "sex", "age_years", "race_eth", "educ2")

# Chunk aggregates buffered before each collapse into the accumulator. Larger
# values re-aggregate less often; the buffer is what bounds peak memory, so it
# stays small enough that the pending list cannot outgrow the accumulator.
CHUNKS_PER_FLUSH <- 10L

# Analysis age spans are defined in 00_functions.r.

# Point-weight equivalence is expected to be exact; the tolerance only keeps a
# last-bit difference in summation order from failing the gate, and the observed
# maximum difference is reported either way.
EQUIVALENCE_TOLERANCE <- 1e-6

# ---- Replicate-panel files ----

raw_dir <- here("01_data", "01_raw", "acs")
proc_dir <- here("01_data", "02_processed", "population")

main_panel_rds <- file.path(proc_dir, "01_pop_main_panel.rds")
replicate_panel_rds <- file.path(proc_dir, "04_pop_replicate_panel.rds")

dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Replicate-panel functions ----

select_vars_to_read <- function(ddi) {
  available_vars <- tolower(ipumsr::ipums_var_info(ddi)$var_name)

  vars_to_read <- c("year", "sex", "age", "race", "hispan", EDUC_VAR, WEIGHT_VAR, REPLICATE_VARS)
  missing_vars <- setdiff(vars_to_read, available_vars)
  if (length(missing_vars)) {
    stop(
      "The ACS extract is missing required IPUMS variables: ",
      paste(toupper(missing_vars), collapse = ", "),
      call. = FALSE
    )
  }
  vars_to_read
}

# Harmonize one raw chunk onto the reference key, one row per person.
prepare_replicate_chunk <- function(chunk_df) {
  in_window <- safe_integer(chunk_df$year) %in% REQUESTED_YEARS
  if (!any(in_window)) {
    return(NULL)
  }
  chunk_df <- chunk_df[in_window, , drop = FALSE]

  chunk_dt <- data.table(
    source = "acs",
    year = safe_integer(chunk_df$year),
    sex = harmonize_acs_sex(chunk_df$sex),
    age_years = harmonize_acs_age(chunk_df$age),
    race_eth = harmonize_acs_race_eth(chunk_df$race, chunk_df$hispan),
    educ2 = recode_educ2_from_educd(chunk_df[[EDUC_VAR]]),
    n_persons_unweighted = 1L,
    n_persons_weighted = safe_numeric(chunk_df[[WEIGHT_VAR]])
  )
  for (replicate_index in seq_len(ACS_REPLICATE_COUNT)) {
    set(
      chunk_dt,
      j = REPLICATE_COLUMNS[replicate_index],
      value = safe_numeric(chunk_df[[REPLICATE_VARS[replicate_index]]])
    )
  }

  # A missing weight contributes nothing to its cell, as in 04.
  for (weight_col in WEIGHT_COLUMNS) {
    missing_weight <- which(is.na(chunk_dt[[weight_col]]))
    if (length(missing_weight)) set(chunk_dt, i = missing_weight, j = weight_col, value = 0)
  }

  chunk_dt
}

aggregate_to_analysis_cells <- function(dt) {
  dt[
    ,
    lapply(.SD, sum),
    by = ANALYSIS_KEY,
    .SDcols = c("n_persons_unweighted", WEIGHT_COLUMNS)
  ]
}

# The bounded reducer: chunk aggregates land in a short pending list and are
# collapsed into a single accumulator, so nothing proportional to file length is
# ever held in memory.
accumulator <- NULL
pending_chunks <- list()

flush_pending_chunks <- function() {
  if (!length(pending_chunks)) {
    return(invisible(NULL))
  }
  combined <- rbindlist(c(list(accumulator), pending_chunks), use.names = TRUE)
  accumulator <<- aggregate_to_analysis_cells(combined)
  pending_chunks <<- list()
  invisible(NULL)
}

# Equivalence gate: aggregating 04's panel onto the reference key must reproduce
# this panel's point-weight columns cell for cell. A difference means the two
# extraction passes disagree about who is in the population.
compare_point_weights <- function(replicate_panel, main_panel) {
  reference <- as.data.table(main_panel)[
    ,
    .(
      n_persons_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
      n_persons_weighted = sum(n_persons_weighted, na.rm = TRUE)
    ),
    by = ANALYSIS_KEY
  ]
  candidate <- replicate_panel[
    , c(ANALYSIS_KEY, "n_persons_unweighted", "n_persons_weighted"), with = FALSE
  ]

  merged <- merge(
    reference, candidate,
    by = ANALYSIS_KEY, all = TRUE, suffixes = c("_reference", "_candidate")
  )
  merged[is.na(n_persons_unweighted_reference), n_persons_unweighted_reference := 0L]
  merged[is.na(n_persons_unweighted_candidate), n_persons_unweighted_candidate := 0L]
  merged[is.na(n_persons_weighted_reference), n_persons_weighted_reference := 0]
  merged[is.na(n_persons_weighted_candidate), n_persons_weighted_candidate := 0]

  merged[, `:=`(
    diff_unweighted = n_persons_unweighted_candidate - n_persons_unweighted_reference,
    diff_weighted = n_persons_weighted_candidate - n_persons_weighted_reference
  )]

  list(
    n_cells_reference = nrow(reference),
    n_cells_candidate = nrow(candidate),
    n_cells_disagreeing = merged[
      abs(diff_unweighted) > 0 | abs(diff_weighted) > EQUIVALENCE_TOLERANCE, .N
    ],
    max_abs_diff_weighted = merged[, max(abs(diff_weighted))]
  )
}

# The two derived bases the rate scripts consume: 05 works off the age-capped
# denominator, 07 off women only, and both bridge child Unknown education into
# No_BA before dividing.
build_death_basis <- function(panel) {
  capped <- copy(panel)
  capped[!is.na(age_years) & age_years > max(AGES_MODEL), age_years := max(AGES_MODEL)]
  bridge_child_unknown_educ2(capped)
}

build_birth_basis <- function(panel) {
  bridge_child_unknown_educ2(panel[sex %chin% "female"])
}

# Every denominator used by the rate layers, stated as an aggregate of panel
# cells. Each target uses the age span of its analysis: 09 estimates shares
# over ages 0-97 and stops on a missing q12, so its target spans 0-97
# even though 10 reads only AGES_MODEL. The two summed-over-education exposures
# carry all five race groups because 05 and 07 hard-stop on them while
# allocating Unknown-education events, not because the endpoints use them.
SUPPORT_TARGETS <- list(
  list(
    id = "death_exposure_cell",
    basis = build_death_basis,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & sex %chin% SEX_LEVELS & age_years %in% AGES_MODEL &
          race_eth %chin% RACES_R5 & educ2 %chin% EDUC_LEVELS
      ]
    },
    key = c("year", "sex", "age_years", "race_eth")
  ),
  list(
    id = "death_denominator_educ",
    basis = build_death_basis,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & sex %chin% SEX_LEVELS & age_years %in% AGES_MODEL &
          race_eth %chin% RACES_R3 & educ2 %chin% EDUC_LEVELS
      ]
    },
    key = c("year", "sex", "age_years", "race_eth", "educ2")
  ),
  list(
    id = "birth_exposure_cell",
    basis = build_birth_basis,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & age_years %in% AGES_FEMALE_REPRO &
          race_eth %chin% RACES_R5 & educ2 %chin% EDUC_LEVELS
      ]
    },
    key = c("year", "age_years", "race_eth")
  ),
  list(
    id = "birth_denominator_educ",
    basis = build_birth_basis,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & age_years %in% AGES_FEMALE_REPRO &
          race_eth %chin% RACES_R3 & educ2 %chin% EDUC_LEVELS
      ]
    },
    key = c("year", "age_years", "race_eth", "educ2")
  ),
  list(
    id = "male_fertility_exposure",
    basis = identity,
    scope = function(dt) {
      scoped <- dt[
        year %in% YEARS_MODEL & sex %chin% "male" & age_years %in% AGES_MALE_REPRO &
          race_eth %chin% RACES_R3 & educ2 %chin% EDUC_LEVELS
      ]
      scoped[, father_age_band := FATHER_AGE_BAND_WIDTH * (age_years %/% FATHER_AGE_BAND_WIDTH)]
      scoped
    },
    key = c("year", "race_eth", "educ2", "father_age_band")
  ),
  list(
    id = "education_share_total",
    basis = identity,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & sex %chin% SEX_LEVELS & age_years %in% AGES_EDUC_TRANSITION &
          race_eth %chin% RACES_R3 & educ2 %chin% EDUC_LEVELS
      ]
    },
    key = c("year", "sex", "age_years", "race_eth")
  ),
  list(
    id = "focal_population",
    basis = identity,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & sex %chin% SEX_LEVELS & age_years %in% AGES_MODEL &
          race_eth %chin% RACES_R3
      ]
    },
    key = c("year", "sex", "age_years", "race_eth")
  ),
  list(
    id = "child_population",
    basis = identity,
    scope = function(dt) {
      dt[
        year %in% YEARS_MODEL & sex %chin% SEX_LEVELS & age_years %in% AGES_MINOR &
          race_eth %chin% RACES_R3
      ]
    },
    key = c("year", "age_years", "race_eth")
  )
)

evaluate_support_target <- function(target, panel) {
  scoped <- target$scope(target$basis(panel))
  group_cols <- target$key
  aggregated <- scoped[, lapply(.SD, sum), by = group_cols, .SDcols = WEIGHT_COLUMNS]

  point <- aggregated$n_persons_weighted
  replicate_matrix <- as.matrix(aggregated[, REPLICATE_COLUMNS, with = FALSE])

  # theta* = theta_0 * exp(...) is positive whenever theta_0 is, so a strictly
  # positive point estimate is necessary and sufficient for the draw to stay in
  # support. The signs are shared within a year, so an aggregated cell is
  # perturbed as one quantity and has to be aggregated before it is judged.
  log_deviations <- acs_log_deviations(point, replicate_matrix)
  sigma <- sqrt(ACS_SDR_SCALE * rowSums(log_deviations^2))
  is_perturbed <- aggregated$year >= ACS_REPLICATE_MIN_YEAR

  # Key values other than year, collapsed into one readable label so cells from
  # targets with different keys can share one evidence file.
  label_cols <- setdiff(target$key, "year")
  cell_label <- do.call(
    paste,
    c(
      lapply(label_cols, function(col) paste0(col, "=", aggregated[[col]])),
      list(sep = " ")
    )
  )

  data.table(
    target_id = target$id,
    year = aggregated$year,
    cell = cell_label,
    theta0 = point,
    log_scale_sd = fifelse(is_perturbed, sigma, 0),
    mean_bias = fifelse(is_perturbed, exp(sigma^2 / 2) - 1, 0),
    # Only meaningful where a log is actually taken: before 2005 every replicate
    # is empty by design and the cell is carried through unperturbed.
    n_replicates_floored = fifelse(
      is_perturbed, rowSums(replicate_matrix < ACS_REPLICATE_ZERO_FLOOR), 0L
    ),
    is_perturbed = is_perturbed,
    passes = is.finite(point) & point > 0
  )
}

# ---- Read and reduce the ACS extract ----

acs_files <- find_acs_files(raw_dir)
ddi <- ipumsr::read_ipums_ddi(acs_files$ddi_file, lower_vars = TRUE)
vars_to_read <- select_vars_to_read(ddi)

callback <- ipumsr::IpumsSideEffectCallback$new(function(chunk, index) {
  chunk_dt <- prepare_replicate_chunk(chunk)
  if (is.null(chunk_dt) || !nrow(chunk_dt)) {
    return(invisible(NULL))
  }

  pending_chunks[[length(pending_chunks) + 1L]] <<- aggregate_to_analysis_cells(chunk_dt)
  if (length(pending_chunks) >= CHUNKS_PER_FLUSH) flush_pending_chunks()
  invisible(NULL)
})

ipumsr::read_ipums_micro_chunked(
  ddi = ddi,
  callback = callback,
  chunk_size = chunk_size,
  vars = tidyselect::all_of(vars_to_read),
  data_file = acs_files$data_file,
  verbose = TRUE,
  var_attrs = c("val_labels")
)

flush_pending_chunks()

if (is.null(accumulator) || !nrow(accumulator)) {
  stop("The ACS extract yielded no person records inside the analysis window.", call. = FALSE)
}

replicate_panel <- accumulator
setorderv(replicate_panel, ANALYSIS_KEY)

# ---- Build the panel (written to a temporary name until both checks pass) ----

# Written to a temporary name first: 14 reads this panel directly and checks
# nothing itself, so a panel that fails the checks below must not be left where a
# draw can pick it up. Same rename-into-place as the draw shards.
replicate_panel_tmp <- paste0(replicate_panel_rds, ".tmp")
saveRDS(replicate_panel, file = replicate_panel_tmp)

# ---- Check 1: the point weights reproduce 04's panel ----

equivalence <- compare_point_weights(replicate_panel, readRDS(main_panel_rds))
equivalence_passed <-
  equivalence$n_cells_disagreeing == 0L &&
  equivalence$n_cells_reference == equivalence$n_cells_candidate

# ---- Check 2: every aggregate the rate layers divide by is positive ----

support <- rbindlist(lapply(SUPPORT_TARGETS, evaluate_support_target, panel = replicate_panel))

support_passed <- support[, all(passes)]

if (!equivalence_passed) {
  stop(
    sprintf(
      paste0(
        "Point-weight equivalence gate failed: %d cells disagree with 04 (max weighted ",
        "difference %.6f). The panel stays staged at %s and is not published."
      ),
      equivalence$n_cells_disagreeing, equivalence$max_abs_diff_weighted,
      replicate_panel_tmp
    ),
    call. = FALSE
  )
}

# theta* keeps the sign of theta_0, so an empty required cell is what can still
# put a draw outside the support. Collapsing cells or reseeding would change the
# estimand, so this stops instead of repairing the panel.
if (!support_passed) {
  stop(
    sprintf(
      paste0(
        "Support gate failed: %d required cells have a non-positive point estimate. ",
        "The panel stays staged at %s and is not published."
      ),
      support[, sum(!passes)], replicate_panel_tmp
    ),
    call. = FALSE
  )
}

# ---- Move the panel into place ----

if (!file.rename(replicate_panel_tmp, replicate_panel_rds)) {
  stop(sprintf("Could not move the replicate panel into place at %s.", replicate_panel_rds), call. = FALSE)
}
