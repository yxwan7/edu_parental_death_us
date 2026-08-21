# ==============================================================================
# 15_aggregate_draws.r
#
# Aggregate a verified draw batch into point estimates, confidence intervals
# and MCSE estimates. The reported point estimate is the draw median; the
# deterministic tables are joined in alongside it for comparison.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/15_aggregate_draws.r")
source(here::here("02_codes", "00_functions.r"))
setDTthreads(1L)

# ---- Batch ----

BATCH_ROOT <- Sys.getenv("EDU_PARENT_LOSS_SHARD_ROOT", unset = "")
if (!nzchar(BATCH_ROOT)) {
  stop("15 requires EDU_PARENT_LOSS_SHARD_ROOT to point at a fetched batch.",
       call. = FALSE)
}

batch_info <- readRDS(file.path(BATCH_ROOT, "batch_info.rds"))
if (!identical(batch_info$batch_mode, "random")) {
  stop("Only a random batch can be aggregated into confidence intervals.", call. = FALSE)
}
M_EXPECTED <- as.integer(batch_info$m_expected)
ALLOW_PARTIAL <- nzchar(Sys.getenv("EDU_PARENT_LOSS_ALLOW_PARTIAL", unset = ""))

# ---- Interval definitions ----

CI_PROBS <- c(0.025, 0.975)
# The reported point estimate is the draw median. The deterministic tables are
# still joined in below, so every output keeps the deterministic value beside
# the median for comparison.
MEDIAN_PROB <- 0.5
MCSE_SEED <- 20260720L
MCSE_BOOTSTRAP_B <- 2000L
MEANAGE_MIN_VALID_SHARE <- 0.95

CI_VALUE_COLUMNS <- list(
  incidence_total = c("n_children_losing_parent", "parental_loss_per_1000"),
  incidence_cause6 = c("n_children_losing_parent", "parental_loss_per_1000"),
  meanage_total = "mean_age_at_loss",
  meanage_cause6 = "mean_age_at_loss"
)

KEY_TOTAL <- c("year", "race_eth", "parent_educ2", "parent_type")
KEY_CAUSE6 <- c(KEY_TOTAL, "cause_6")
KEY_GAP <- c("year", "race_eth", "parent_type", "gap_metric")
KEY_MEANAGE_GAP <- c("year", "race_eth", "parent_type")
MAIN_TEXT_KEY_YEARS <- c(2010L, 2019L, 2021L, 2024L)
MAIN_TEXT_PERIOD_LABEL <- "2010-2024"
MAIN_TEXT_COUNT_KEY <- c("period", "race_eth", "parent_type", "parent_educ2")

SHARD_KEYS <- list(
  incidence_total = KEY_TOTAL,
  incidence_cause6 = KEY_CAUSE6,
  meanage_total = KEY_TOTAL,
  meanage_cause6 = KEY_CAUSE6
)

point_estimate_dir <- here("01_data", "02_processed", "edu_parent_loss")
point_estimate_paths <- c(
  incidence_total = file.path(point_estimate_dir, "point_estimate_incidence_total.rds"),
  incidence_cause6 = file.path(point_estimate_dir, "point_estimate_incidence_cause6.rds"),
  meanage_total = file.path(point_estimate_dir, "point_estimate_meanage_total.rds"),
  meanage_cause6 = file.path(point_estimate_dir, "point_estimate_meanage_cause6.rds")
)
if (any(!file.exists(point_estimate_paths))) {
  stop("One or more frozen point-estimate tables are missing.", call. = FALSE)
}

# ---- Draw shards ----

# One row per shard, carrying the metadata the audit trail and the batch gate
# below both read. 14 validates the endpoint tables before it writes a shard.
audit_row <- function(shard, path) {
  meta <- shard$meta
  field <- function(name, default) uncertainty_field(meta, name, default)
  data.table(
    file = basename(path),
    draw_id = as.character(field("draw_id", NA_character_)),
    mode = as.character(field("mode", NA_character_)),
    ledger_index = as.integer(field("ledger_index", NA_integer_)),
    master_seed = as.numeric(field("master_seed", NA_real_)),
    code_version = as.character(field("code_version", NA_character_)),
    batch_id = as.character(field("batch_id", NA_character_)),
    status = as.character(field("status", NA_character_)),
    n_deaths_total = as.numeric(field("n_deaths_total", NA_real_)),
    n_births_total = as.numeric(field("n_births_total", NA_real_)),
    theta_factor_min = as.numeric(field("theta_factor_min", NA_real_)),
    theta_factor_max = as.numeric(field("theta_factor_max", NA_real_)),
    intensity_m_min = as.numeric(field("intensity_m_min", NA_real_)),
    intensity_m_max = as.numeric(field("intensity_m_max", NA_real_)),
    intensity_kappa_min = as.numeric(field("intensity_kappa_min", NA_real_)),
    intensity_kappa_max = as.numeric(field("intensity_kappa_max", NA_real_)),
    intensity_s_birth = as.numeric(field("intensity_s_birth", NA_real_))
  )
}

shard_files <- list.files(
  file.path(BATCH_ROOT, "shards"),
  pattern = "^draw_[0-9]+\\.rds$",
  full.names = TRUE
)
shards <- lapply(shard_files, readRDS)
audit <- rbindlist(Map(audit_row, shards, shard_files))
# Name the shards by ledger position before the audit is reordered, so the two
# are matched by position rather than by file-name order.
names(shards) <- as.character(audit$ledger_index)
setorder(audit, ledger_index)

# ---- Batch gate ----

ledger_positions <- audit$ledger_index
if (!nrow(audit)) {
  stop("The batch contains no shards.", call. = FALSE)
}
if (anyNA(ledger_positions) || anyDuplicated(ledger_positions) ||
    any(!ledger_positions %between% c(1L, M_EXPECTED))) {
  stop("Every shard must hold one distinct ledger position within 1..M.",
       call. = FALSE)
}
if (!all(audit$code_version == batch_info$code_version) ||
    !all(audit$master_seed == batch_info$master_seed)) {
  stop("The batch mixes code versions or master seeds.", call. = FALSE)
}

is_final <- setequal(ledger_positions, seq_len(M_EXPECTED))
if (!is_final && !ALLOW_PARTIAL) {
  stop(sprintf(
    "%d of %d ledger positions are complete. Set EDU_PARENT_LOSS_ALLOW_PARTIAL=1 only for a provisional diagnostic.",
    length(ledger_positions), M_EXPECTED
  ), call. = FALSE)
}

order_ids <- ledger_positions
shard_tables <- lapply(shards[as.character(order_ids)], `[[`, "tables")
n_complete <- length(shard_tables)

# ---- Endpoint functions ----

stack_endpoint <- function(tables, table_name, value_column) {
  key <- SHARD_KEYS[[table_name]]
  ordered <- lapply(tables, function(shard_tables) {
    dt <- as.data.table(shard_tables[[table_name]])
    setorderv(dt, key)
    dt
  })
  reference <- ordered[[1]][, key, with = FALSE]
  values <- vapply(
    ordered,
    function(dt) as.numeric(dt[[value_column]]),
    numeric(nrow(reference))
  )
  list(key = reference, values = matrix(values, nrow = nrow(reference)))
}

mcse_position_cache <- new.env(parent = emptyenv())

bootstrap_positions <- function(n_draws, probs, n_bootstrap, seed) {
  cache_key <- paste(n_draws, paste(probs, collapse = ","), n_bootstrap, seed, sep = ":")
  if (exists(cache_key, envir = mcse_position_cache, inherits = FALSE)) {
    return(get(cache_key, envir = mcse_position_cache, inherits = FALSE))
  }

  set.seed(seed)

  h <- (n_draws - 1) * probs + 1
  lo_rank <- floor(h)
  hi_rank <- pmin(lo_rank + 1L, n_draws)
  positions <- matrix(NA_integer_, nrow = n_bootstrap, ncol = 4L)
  ranks <- c(lo_rank[1], hi_rank[1], lo_rank[2], hi_rank[2])

  for (b in seq_len(n_bootstrap)) {
    counts <- tabulate(
      sample.int(n_draws, n_draws, replace = TRUE),
      nbins = n_draws
    )
    cumulative <- cumsum(counts)
    positions[b, ] <- vapply(
      ranks,
      function(rank) findInterval(rank - 1L, cumulative) + 1L,
      integer(1)
    )
  }

  result <- list(positions = positions, weight = h - lo_rank)
  assign(cache_key, result, envir = mcse_position_cache)
  result
}

bootstrap_quantile_mcse <- function(
    sorted_values, n_draws, probs = CI_PROBS,
    n_bootstrap = MCSE_BOOTSTRAP_B, seed = MCSE_SEED) {
  draws <- bootstrap_positions(n_draws, probs, n_bootstrap, seed)
  position <- draws$positions
  weight <- draws$weight

  lower_replicates <-
    (1 - weight[1]) * sorted_values[, position[, 1], drop = FALSE] +
    weight[1] * sorted_values[, position[, 2], drop = FALSE]
  upper_replicates <-
    (1 - weight[2]) * sorted_values[, position[, 3], drop = FALSE] +
    weight[2] * sorted_values[, position[, 4], drop = FALSE]

  list(
    mcse_lower = apply(lower_replicates, 1L, sd),
    mcse_upper = apply(upper_replicates, 1L, sd)
  )
}

summarise_complete_endpoint <- function(values, label) {
  if (!all(is.finite(values))) {
    stop(sprintf(
      "%s has %d non-finite draw values.", label, sum(!is.finite(values))
    ), call. = FALSE)
  }
  sorted_values <- t(apply(values, 1L, sort))
  n_draws <- ncol(values)
  quantiles <- t(apply(
    sorted_values, 1L, quantile,
    probs = c(CI_PROBS[1], MEDIAN_PROB, CI_PROBS[2]), type = 7, names = FALSE
  ))
  # MCSE stays on the two interval bounds; bootstrap_quantile_mcse() assumes a
  # length-2 probs vector, so the median is taken separately above.
  mcse <- bootstrap_quantile_mcse(sorted_values, n_draws)
  data.table(
    median = quantiles[, 2],
    lower_ci = quantiles[, 1],
    upper_ci = quantiles[, 3],
    mcse_lower = mcse$mcse_lower,
    mcse_upper = mcse$mcse_upper,
    n_valid_draws = n_draws,
    p_zero = 0
  )
}

# `point_estimate` is the deterministic mean age. It no longer supplies the
# reported value -- the median does -- but it still gates which cells are
# summarised at all: a cell the deterministic run leaves undefined has no
# events to average, so the median and the interval are both withheld.
summarise_meanage_endpoint <- function(values, point_estimate, n_complete) {
  out <- data.table(
    median = rep(NA_real_, nrow(values)),
    lower_ci = NA_real_,
    upper_ci = NA_real_,
    mcse_lower = NA_real_,
    mcse_upper = NA_real_,
    n_valid_draws = rowSums(is.finite(values)),
    p_zero = NA_real_
  )
  out[, p_zero := 1 - n_valid_draws / n_complete]

  usable <- which(
    is.finite(point_estimate) &
      out$n_valid_draws >= MEANAGE_MIN_VALID_SHARE * n_complete &
      out$n_valid_draws > 1L
  )
  for (cell in usable) {
    valid <- sort(values[cell, ][is.finite(values[cell, ])])
    quantiles <- quantile(
      valid, probs = c(CI_PROBS[1], MEDIAN_PROB, CI_PROBS[2]),
      type = 7, names = FALSE
    )
    mcse <- bootstrap_quantile_mcse(matrix(valid, nrow = 1L), length(valid))
    set(out, i = cell, j = "median", value = quantiles[2])
    set(out, i = cell, j = "lower_ci", value = quantiles[1])
    set(out, i = cell, j = "upper_ci", value = quantiles[3])
    set(out, i = cell, j = "mcse_lower", value = mcse$mcse_lower)
    set(out, i = cell, j = "mcse_upper", value = mcse$mcse_upper)
  }
  out
}

# Both gap tables carry their point estimate in a column called `value`, so the
# median is named to match rather than left as a bare `median`, which would
# shadow stats::median wherever the table is read downstream.
summarise_gap_endpoint <- function(values, label) {
  out <- summarise_complete_endpoint(values, label)
  setnames(out, "median", "value_median")
  out[]
}

derive_gap_draws <- function(tables) {
  per_draw <- lapply(tables, function(shard_tables) {
    gap <- build_gap_data(as.data.table(shard_tables$incidence_total))
    setorderv(gap, KEY_GAP)
    gap
  })
  reference <- per_draw[[1]][, KEY_GAP, with = FALSE]
  values <- vapply(
    per_draw, function(dt) as.numeric(dt$value), numeric(nrow(reference))
  )
  list(key = reference, values = matrix(values, nrow = nrow(reference)))
}

build_meanage_gap <- function(meanage_total_dt) {
  wide <- dcast(
    meanage_total_dt[parent_educ2 %chin% c("No_BA", "BA_plus")],
    year + race_eth + parent_type ~ parent_educ2,
    value.var = "mean_age_at_loss"
  )
  wide[, .(
    year, race_eth, parent_type,
    value = No_BA - BA_plus
  )]
}

derive_meanage_gap_draws <- function(tables) {
  per_draw <- lapply(tables, function(shard_tables) {
    gap <- build_meanage_gap(as.data.table(shard_tables$meanage_total))
    setorderv(gap, KEY_MEANAGE_GAP)
    gap
  })
  reference <- per_draw[[1]][, KEY_MEANAGE_GAP, with = FALSE]
  values <- vapply(
    per_draw, function(dt) as.numeric(dt$value), numeric(nrow(reference))
  )
  list(key = reference, values = matrix(values, nrow = nrow(reference)))
}

build_bereaved_children_number_table <- function(incidence_total_dt) {
  dt <- copy(as.data.table(incidence_total_dt))[
    year %between% c(2010L, 2024L) &
      parent_educ2 %chin% c("Total", "No_BA", "BA_plus"),
    .(year, race_eth, parent_type, parent_educ2, n_children_losing_parent)
  ]

  annual <- dt[year %in% MAIN_TEXT_KEY_YEARS, .(
    period = as.character(year), race_eth, parent_type, parent_educ2,
    n_children_losing_parent
  )]
  full_period <- dt[, .(
    n_children_losing_parent = sum(n_children_losing_parent)
  ), by = .(race_eth, parent_type, parent_educ2)]
  full_period[, period := MAIN_TEXT_PERIOD_LABEL]
  setcolorder(full_period, names(annual))

  out <- rbindlist(list(annual, full_period), use.names = TRUE)
  totals <- out[parent_educ2 == "Total", .(
    period, race_eth, parent_type,
    total_n_children_losing_parent = n_children_losing_parent
  )]
  out <- merge(
    out, totals,
    by = c("period", "race_eth", "parent_type"),
    all.x = TRUE,
    sort = FALSE
  )
  if (any(!is.finite(out$total_n_children_losing_parent)) ||
      any(out$total_n_children_losing_parent <= 0)) {
    stop("Main-text number table has a missing or non-positive total.",
         call. = FALSE)
  }
  out[, share_percent := fifelse(
    parent_educ2 == "Total",
    NA_real_,
    100 * n_children_losing_parent / total_n_children_losing_parent
  )]
  out[, total_n_children_losing_parent := NULL]
  out[, period_order := match(
    period, c(as.character(MAIN_TEXT_KEY_YEARS), MAIN_TEXT_PERIOD_LABEL)
  )]
  out[, education_order := match(
    parent_educ2, c("Total", "No_BA", "BA_plus")
  )]
  setorder(out, parent_type, race_eth, education_order, period_order)
  out[, c("period_order", "education_order") := NULL]
  out[]
}

derive_bereaved_children_number_draws <- function(tables) {
  per_draw <- lapply(tables, function(shard_tables) {
    out <- build_bereaved_children_number_table(shard_tables$incidence_total)
    setorderv(out, MAIN_TEXT_COUNT_KEY)
    out
  })

  reference <- per_draw[[1]][, MAIN_TEXT_COUNT_KEY, with = FALSE]
  count_values <- vapply(
    per_draw,
    function(dt) as.numeric(dt$n_children_losing_parent),
    numeric(nrow(reference))
  )

  share_reference <- per_draw[[1]][parent_educ2 != "Total"]
  share_key <- share_reference[, MAIN_TEXT_COUNT_KEY, with = FALSE]
  share_values <- vapply(
    per_draw,
    function(dt) as.numeric(dt[parent_educ2 != "Total"]$share_percent),
    numeric(nrow(share_reference))
  )

  list(
    count_key = reference,
    count_values = matrix(count_values, nrow = nrow(reference)),
    share_key = share_key,
    share_values = matrix(share_values, nrow = nrow(share_reference))
  )
}

# ---- Summaries ----

ci_tables <- list()
for (table_name in names(CI_VALUE_COLUMNS)) {
  point_estimate <- as.data.table(readRDS(point_estimate_paths[[table_name]]))
  point_order <- copy(point_estimate)[, SHARD_KEYS[[table_name]], with = FALSE]
  setorderv(point_estimate, SHARD_KEYS[[table_name]])
  summaries <- list()

  for (value_column in CI_VALUE_COLUMNS[[table_name]]) {
    stacked <- stack_endpoint(shard_tables, table_name, value_column)
    summary_dt <- if (identical(value_column, "mean_age_at_loss")) {
      summarise_meanage_endpoint(
        stacked$values, point_estimate[[value_column]], n_complete
      )
    } else {
      summarise_complete_endpoint(
        stacked$values, paste(table_name, value_column, sep = "/")
      )
    }
    setnames(summary_dt, paste(value_column, names(summary_dt), sep = "_"))
    summaries[[value_column]] <- cbind(stacked$key, summary_dt)
  }

  merged <- Reduce(
    function(a, b) merge(a, b, by = SHARD_KEYS[[table_name]]),
    summaries
  )
  joined <- merge(point_estimate, merged, by = SHARD_KEYS[[table_name]])
  if (nrow(joined) != nrow(point_estimate)) {
    stop(table_name, ": the interval join dropped point-estimate rows.", call. = FALSE)
  }
  ci_tables[[table_name]] <- joined[point_order, on = SHARD_KEYS[[table_name]]]
}

gap_stacked <- derive_gap_draws(shard_tables)
gap_point <- build_gap_data(as.data.table(
  readRDS(point_estimate_paths[["incidence_total"]])
))
gap_order <- copy(gap_point)[, KEY_GAP, with = FALSE]
setorderv(gap_point, KEY_GAP)
gap_table <- merge(
  gap_point,
  cbind(
    gap_stacked$key,
    summarise_gap_endpoint(gap_stacked$values, "incidence_gap/value")
  ),
  by = KEY_GAP
)[gap_order, on = KEY_GAP]

meanage_gap_stacked <- derive_meanage_gap_draws(shard_tables)
meanage_gap_point <- build_meanage_gap(as.data.table(
  readRDS(point_estimate_paths[["meanage_total"]])
))
meanage_gap_order <- copy(meanage_gap_point)[, KEY_MEANAGE_GAP, with = FALSE]
setorderv(meanage_gap_point, KEY_MEANAGE_GAP)
meanage_gap_table <- merge(
  meanage_gap_point,
  cbind(
    meanage_gap_stacked$key,
    summarise_gap_endpoint(
      meanage_gap_stacked$values, "meanage_gap/value"
    )
  ),
  by = KEY_MEANAGE_GAP
)[meanage_gap_order, on = KEY_MEANAGE_GAP]

main_number_stacked <- derive_bereaved_children_number_draws(shard_tables)
main_number_point <- build_bereaved_children_number_table(as.data.table(
  readRDS(point_estimate_paths[["incidence_total"]])
))
main_number_order <- copy(main_number_point)[, MAIN_TEXT_COUNT_KEY, with = FALSE]
setorderv(main_number_point, MAIN_TEXT_COUNT_KEY)

main_number_count_summary <- summarise_complete_endpoint(
  main_number_stacked$count_values,
  "main_text_number_table/n_children_losing_parent"
)
setnames(
  main_number_count_summary,
  names(main_number_count_summary),
  paste0("n_children_losing_parent_", names(main_number_count_summary))
)
main_number_table <- merge(
  main_number_point,
  cbind(main_number_stacked$count_key, main_number_count_summary),
  by = MAIN_TEXT_COUNT_KEY,
  all.x = TRUE,
  sort = FALSE
)

main_number_share_summary <- summarise_complete_endpoint(
  main_number_stacked$share_values,
  "main_text_number_table/share_percent"
)
setnames(
  main_number_share_summary,
  names(main_number_share_summary),
  paste0("share_percent_", names(main_number_share_summary))
)
main_number_table <- merge(
  main_number_table,
  cbind(main_number_stacked$share_key, main_number_share_summary),
  by = MAIN_TEXT_COUNT_KEY,
  all.x = TRUE,
  sort = FALSE
)[main_number_order, on = MAIN_TEXT_COUNT_KEY]

# ---- Output bundle ----

batch_summary <- list(
  batch_id = batch_info$batch_id,
  m_complete = n_complete,
  m_expected = M_EXPECTED,
  ledger_positions = order_ids,
  is_final = is_final,
  mode = "random",
  master_seed = batch_info$master_seed,
  code_version = batch_info$code_version,
  point_estimate_rule = "draw_median",
  ci_probs = CI_PROBS,
  mcse_seed = MCSE_SEED,
  mcse_bootstrap_b = MCSE_BOOTSTRAP_B,
  aggregated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

output_root <- Sys.getenv("EDU_PARENT_LOSS_CI_ROOT", unset = "")
if (!nzchar(output_root)) {
  output_root <- here(
    "01_data", "02_processed", "edu_parent_loss", "uncertainty",
    batch_info$batch_id,
    if (is_final) "final" else "provisional"
  )
}
output_root <- normalizePath(output_root, winslash = "/", mustWork = FALSE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

ci_output_paths <- c(
  incidence_total = file.path(output_root, "15_ci_incidence_total.rds"),
  incidence_cause6 = file.path(output_root, "15_ci_incidence_cause6.rds"),
  meanage_total = file.path(output_root, "15_ci_meanage_total.rds"),
  meanage_cause6 = file.path(output_root, "15_ci_meanage_cause6.rds")
)

for (table_name in names(ci_tables)) {
  saveRDS(ci_tables[[table_name]], ci_output_paths[[table_name]])
}
saveRDS(gap_table, file.path(output_root, "15_ci_incidence_gap.rds"))
saveRDS(
  meanage_gap_table,
  file.path(output_root, "15_ci_meanage_gap.rds")
)
saveRDS(
  main_number_table,
  file.path(output_root, "15_ci_bereaved_children_number_table.rds")
)
saveRDS(
  list(
    audit = audit,
    batch_summary = batch_summary,
    batch_info = batch_info
  ),
  file.path(output_root, "15_uncertainty_batch_record.rds")
)

cat(sprintf(
  "Batch %s: %d of %d draws, %s output at %s\n",
  batch_info$batch_id, n_complete, M_EXPECTED,
  if (is_final) "final" else "PROVISIONAL", output_root
))
