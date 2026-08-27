# ==============================================================================
# 02_run_draw.r
#
# Generate one uncertainty draw in a temporary project and save its four
# endpoint tables to the shard directory.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
  library(splines)
})

here::i_am("02_codes/uncertainty/02_run_draw.r")
source(here::here("02_codes", "00_functions.r"))
repo_root <- here::here()

# ---- Draw configuration ----

DRAW_ROOT <- Sys.getenv("EDU_PARENT_LOSS_DRAW_ROOT", unset = "")
BATCH_ROOT <- Sys.getenv("EDU_PARENT_LOSS_SHARD_ROOT", unset = "")
BATCH_ID <- Sys.getenv("EDU_PARENT_LOSS_BATCH_ID", unset = "")
DRAW_ID <- Sys.getenv("EDU_PARENT_LOSS_DRAW_ID", unset = "")
DRAW_MODE <- Sys.getenv("EDU_PARENT_LOSS_DRAW_MODE", unset = "")

batch_info <- readRDS(file.path(BATCH_ROOT, "batch_info.rds"))
seed_ledger <- readRDS(file.path(BATCH_ROOT, "seed_ledger.rds"))
CODE_VERSION <- current_commit()

# A batch is only aggregated when every shard carries one code version, so a
# draw that would write a mixed shard stops before it spends a node-hour.
if (!identical(batch_info$batch_id, BATCH_ID) ||
    !identical(CODE_VERSION, batch_info$code_version)) {
  stop("The batch ID or code version does not match the batch record.",
       call. = FALSE)
}

DRAW_INDEX <- if (identical(DRAW_MODE, "random")) {
  as.integer(DRAW_ID)
} else {
  NA_integer_
}

# Inputs perturbed for each draw.
PERTURBED_INPUTS <- list(
  c("01_data", "02_processed", "deaths", "01_death_main_panel.rds"),
  c("01_data", "02_processed", "births", "01_birth_main_panel.rds"),
  c("01_data", "02_processed", "population", "01_pop_main_panel.rds"),
  c("01_data", "02_processed", "population", "02_pop_denominator_death.rds"),
  c("01_data", "02_processed", "population", "03_pop_denominator_birth.rds")
)

# Year-level intensity layer.
INTENSITY_AGES <- 18:59
INTENSITY_YEARS <- setdiff(1999:2024, c(2020L, 2021L))
INTENSITY_DF <- 4L
INTENSITY_EIGEN_FLOOR <- 0.05
INTENSITY_BIRTH_DF <- 8L
INTENSITY_BIRTH_YEARS <- 1990:2024
# covid19 borrows the residual category, whose behaviour it most resembles.
INTENSITY_COVID_DONOR <- "all_other_causes"

# Fixed inputs used to construct the ACS draw.
replicate_panel_path <- file.path(repo_root, "01_data", "02_processed", "population",
  "04_pop_replicate_panel.rds"
)
main_panel_path <- file.path(repo_root, "01_data", "02_processed", "population",
  "01_pop_main_panel.rds"
)

# Point-estimate stages run in the temporary project.
STAGES <- c(
  "05_mortality_rates.r",
  "06_mortality_rates_by_cause.r",
  "07_fertility_rates.r",
  "08_male_fertility_rates.r",
  "09_education_transitions.r",
  "10_kinship_model.r",
  "11_parental_death.r"
)

RSCRIPT <- file.path(
  R.home("bin"),
  if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
)

setDTthreads(1L)

shard_stem <- if (identical(DRAW_MODE, "random")) {
  sprintf("draw_%04d", DRAW_INDEX)
} else {
  paste0("draw_", DRAW_ID)
}

# ---- Scratch and persistent paths ----

# The submit wrapper creates DRAW_ROOT with mktemp and removes it on exit.
shard_dir <- file.path(BATCH_ROOT, "shards")
failure_dir <- file.path(BATCH_ROOT, "failures")
dir.create(shard_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(failure_dir, recursive = TRUE, showWarnings = FALSE)

analysis_root <- file.path(DRAW_ROOT, "analysis")
analysis_path <- function(...) file.path(analysis_root, ...)
repository_path <- function(parts) {
  do.call(file.path, c(list(repo_root), as.list(parts)))
}
draw_path <- function(parts) {
  do.call(file.path, c(list(analysis_root), as.list(parts)))
}

shard_final_path <- file.path(shard_dir, paste0(shard_stem, ".rds"))
shard_fail_path <- file.path(failure_dir, paste0(shard_stem, ".fail"))

# Leave one small failure record per unresolved draw. Child-process log tails
# are printed so they remain in the Slurm log after scratch is removed.
record_failed_draw <- function(stage, exit_code, message_text,
                               log_files = character(0)) {
  tail_lines <- unlist(lapply(log_files[file.exists(log_files)], function(path) {
    utils::tail(readLines(path, warn = FALSE), 40L)
  }))
  if (length(tail_lines)) {
    message(paste(tail_lines, collapse = "\n"))
  }
  record <- list(
    batch_id = BATCH_ID,
    code_version = CODE_VERSION,
    draw_id = DRAW_ID,
    mode = DRAW_MODE,
    ledger_index = DRAW_INDEX, master_seed = DRAW_MASTER_SEED,
    stage = stage, exit_code = exit_code, message = message_text
  )
  temporary <- tempfile(paste0(shard_stem, "_"), failure_dir, ".tmp")
  saveRDS(record, temporary)
  if (!file.rename(temporary, shard_fail_path)) {
    unlink(temporary)
    warning("Could not commit failure record ", shard_fail_path, call. = FALSE)
  }
}

# Any error from here on leaves a failure record before the job exits, so an
# unresolved draw is never silent in the batch root.
options(error = function() {
  record_failed_draw("draw", NA_integer_, geterrmessage())
  quit(status = 1L, save = "no")
})

dir.create(analysis_path("02_codes"), recursive = TRUE, showWarnings = FALSE)
invisible(file.create(analysis_path(".here")))
invisible(file.copy(
  file.path(repo_root, "02_codes", c("00_functions.r", STAGES)),
  analysis_path("02_codes"),
  overwrite = TRUE
))

# ---- Seeds and the perturbations a draw applies ----

use_rng_stream <- function(stream) {
  RNGkind(DRAW_RNG_KIND)
  assign(".Random.seed", stream, envir = .GlobalEnv)
  invisible(NULL)
}

# Poisson bootstrap of an event count panel, applied row by row: independent
# Poisson draws sum to a Poisson of the summed rate, so the result matches a draw
# taken on any coarser grid downstream. Both parents' events are drawn in one
# pass, which preserves their covariance, and the intensity multiplier scales the
# mean in the same rpois call.
poisson_resample_counts <- function(counts, multiplier = 1) {
  observed <- as.numeric(counts) * as.numeric(multiplier)
  if (anyNA(observed) || any(observed < 0)) {
    stop("Event counts must be non-negative and non-missing before resampling.", call. = FALSE)
  }
  as.integer(rpois(length(observed), observed))
}

# ---- Year-level intensity layer ----

detrended_series <- function(series, df_trend) {
  series <- series[order(series$year)]
  fit <- glm(n ~ ns(year, df = df_trend), family = quasipoisson(), data = series)
  fitted_values <- as.numeric(fitted(fit))
  list(
    observed = as.numeric(series$n),
    fitted = fitted_values,
    rel = (series$n - fitted_values) / fitted_values,
    df_resid = fit$df.residual
  )
}

# Net of the Poisson component that the draw's rpois adds, so it is not counted
# twice. df.residual, not n - 1, charges the spline's parameters.
excess_variance <- function(rel, df_resid, mu, floor = 0) {
  rel <- as.matrix(rel)
  pmax(colSums(rel^2) / df_resid - 1 / mu, floor)
}

# Dispersion and cross-cause correlation from the unperturbed panel; kappa makes
# the implied all-cause dispersion reproduce the directly measured one.
calibrate_death_intensity <- function(death_panel) {
  usable <- death_panel[
    !is.na(age_years) & age_years %in% INTENSITY_AGES &
      race_eth %chin% RACES_R3 & year %in% INTENSITY_YEARS &
      !is.na(cause_6) & cause_6 != "covid19",
    .(n = sum(n_deaths)),
    by = .(race_eth, sex, cause_6, year)
  ]
  causes <- sort(unique(usable$cause_6))
  groups <- unique(usable[, .(race_eth, sex)])

  correlations <- vector("list", nrow(groups))
  dispersions <- vector("list", nrow(groups))
  for (i in seq_len(nrow(groups))) {
    group <- usable[groups[i], on = .(race_eth, sex)]
    fits <- sapply(causes, function(cause) {
      detrended_series(group[cause_6 == cause, .(year, n)], INTENSITY_DF)
    }, simplify = FALSE)
    rel_matrix <- sapply(fits, `[[`, "rel")
    fitted_matrix <- sapply(fits, `[[`, "fitted")
    observed_matrix <- sapply(fits, `[[`, "observed")
    df_resid <- fits[[1]]$df_resid

    covariance <- crossprod(rel_matrix) / df_resid
    diag(covariance) <- excess_variance(rel_matrix, df_resid, colMeans(fitted_matrix),
                                        floor = 1e-10)
    sd_by_cause <- setNames(sqrt(diag(covariance)), causes)
    correlations[[i]] <- covariance / tcrossprod(sd_by_cause)

    # Year-specific weights: the cause mix shifts markedly across the window.
    fitted_totals <- rowSums(fitted_matrix)
    rel_all <- rowSums(observed_matrix - fitted_matrix) / fitted_totals
    dispersions[[i]] <- list(
      race_eth = groups$race_eth[i], sex = groups$sex[i],
      s = sd_by_cause,
      year_weights = fitted_matrix / fitted_totals,
      s_all = sqrt(excess_variance(rel_all, df_resid, mean(fitted_totals)))
    )
  }

  correlation <- Reduce(`+`, correlations) / length(correlations)
  diag(correlation) <- 1

  # The copied row is rank deficient, so the positive-definite floor runs after
  # it is appended -- otherwise chol() fails.
  donor_row <- correlation[INTENSITY_COVID_DONOR, ]
  correlation <- rbind(cbind(correlation, covid19 = donor_row),
                       covid19 = c(donor_row, 1))
  all_causes <- c(causes, "covid19")

  spectrum <- eigen(correlation, symmetric = TRUE)
  if (min(spectrum$values) < INTENSITY_EIGEN_FLOOR) {
    correlation <- spectrum$vectors %*%
      diag(pmax(spectrum$values, INTENSITY_EIGEN_FLOOR)) %*% t(spectrum$vectors)
    correlation <- correlation / tcrossprod(sqrt(diag(correlation)))
  }
  dimnames(correlation) <- list(all_causes, all_causes)

  # Solved on the calibrated causes, with the weights that produced the target.
  calibrated <- correlation[causes, causes]
  sd_by_group <- rbindlist(lapply(dispersions, function(entry) {
    weighted <- sweep(entry$year_weights, 2, entry$s, `*`)
    implied <- sqrt(mean(rowSums((weighted %*% calibrated) * weighted)))
    kappa <- if (implied > 0) entry$s_all / implied else 1
    scaled <- kappa * entry$s
    data.table(race_eth = entry$race_eth, sex = entry$sex, cause_6 = all_causes,
               s = c(scaled, scaled[[INTENSITY_COVID_DONOR]]),
               kappa = kappa, s_all = entry$s_all)
  }))

  list(correlation = correlation, sd = sd_by_group, causes = all_causes)
}

calibrate_birth_intensity <- function(birth_panel) {
  usable <- birth_panel[
    resident_birth_flag & mother_race_eth %chin% RACES_R3 &
      year %in% INTENSITY_BIRTH_YEARS,
    .(n = sum(n_births)),
    by = .(mother_race_eth, year)
  ]
  by_race <- vapply(RACES_R3, function(race) {
    parts <- detrended_series(usable[mother_race_eth == race, .(year, n)],
                              INTENSITY_BIRTH_DF)
    sqrt(excess_variance(parts$rel, parts$df_resid, mean(parts$fitted)))
  }, numeric(1))
  median(by_race)
}

# One correlated vector per year and race x sex group, shared across ages and both
# education legs. meanlog = -sigma^2/2 gives E[m] = 1. Years held out of the
# calibration still get multipliers.
draw_death_multipliers <- function(calibration, years) {
  causes <- calibration$causes
  chol_factor <- chol(calibration$correlation)
  # Fixed order, so a retried draw on the same stream reproduces the multipliers.
  groups <- setorder(unique(calibration$sd[, .(race_eth, sex)]), race_eth, sex)
  years <- sort(years)

  rbindlist(lapply(seq_len(nrow(groups)), function(i) {
    sd_by_cause <- calibration$sd[groups[i], on = .(race_eth, sex)]
    sd_by_cause <- setNames(sd_by_cause$s, sd_by_cause$cause_6)[causes]
    z <- matrix(rnorm(length(causes) * length(years)), nrow = length(causes))
    log_m <- crossprod(chol_factor, z) * sd_by_cause - sd_by_cause^2 / 2
    data.table(year = rep(years, each = length(causes)),
               race_eth = groups$race_eth[i], sex = groups$sex[i],
               cause_6 = causes, m = as.vector(exp(log_m)))
  }))
}

# One Rademacher vector per year, shared across every cell of that year; years
# without published replicate weights get no vector and stay unperturbed.
draw_acs_epsilon <- function(years) {
  perturbed_years <- sort(unique(years[years >= ACS_REPLICATE_MIN_YEAR]))
  epsilon <- lapply(perturbed_years, function(year) {
    sample(c(-1, 1), ACS_REPLICATE_COUNT, replace = TRUE)
  })
  setNames(epsilon, as.character(perturbed_years))
}

# theta* as a per-cell multiplier, so it can be applied to the finer
# point-estimate panel. That keeps the 04_process_population.r schema and its
# raw code columns intact and leaves n_persons_unweighted untouched, while the
# whole ACS bundle derives from one theta*.
acs_draw_factors <- function(replicate_panel, epsilon_by_year,
                             replicate_columns = acs_replicate_columns()) {
  point <- replicate_panel$n_persons_weighted
  replicate_matrix <- as.matrix(replicate_panel[, replicate_columns, with = FALSE])
  factors <- rep(1, length(point))

  for (year_label in names(epsilon_by_year)) {
    rows <- which(replicate_panel$year == as.integer(year_label))
    if (!length(rows)) next
    log_deviations <- acs_log_deviations(point[rows], replicate_matrix[rows, , drop = FALSE])
    factors[rows] <- exp(
      sqrt(ACS_SDR_SCALE) *
        as.vector(log_deviations %*% epsilon_by_year[[year_label]])
    )
  }
  factors
}

# ---- Draw inputs ----

write_identity_inputs <- function() {
  for (parts in PERTURBED_INPUTS) {
    dst <- draw_path(parts)
    dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
    file.copy(repository_path(parts), dst, overwrite = TRUE)
  }
  # Identity draws keep the same summary fields without reading the replicate
  # panel, and fix the intensity multipliers at one.
  list(
    n_epsilon_years = 0L,
    n_deaths_total = NA_real_,
    n_births_total = NA_real_,
    theta_factor_min = 1,
    theta_factor_max = 1,
    intensity_kappa_min = NA_real_,
    intensity_kappa_max = NA_real_,
    intensity_m_min = 1,
    intensity_m_max = 1,
    intensity_s_birth = NA_real_
  )
}

# Poisson event draws and the ACS replicate-weight draw use separate substreams;
# the intensity multipliers come from the same substreams, before the counts.
write_random_inputs <- function(streams) {
  count_totals <- c(n_deaths_total = NA_real_, n_births_total = NA_real_)

  # Deaths: calibrate on the unperturbed panel, then draw intensity x Poisson.
  death_panel <- as.data.table(readRDS(repository_path(PERTURBED_INPUTS[[1]])))
  death_panel[, cause_6 := classify_cause_6_from_cause_atom(cause_atom)]
  death_calibration <- calibrate_death_intensity(death_panel)
  use_rng_stream(streams$death)
  multiplier_years <- sort(unique(death_panel[!is.na(cause_6), year]))
  multipliers <- draw_death_multipliers(death_calibration, multiplier_years)
  death_panel[multipliers, m := i.m, on = .(year, race_eth, sex, cause_6)]
  death_panel[is.na(m), m := 1]
  set(death_panel, j = "n_deaths",
      value = poisson_resample_counts(death_panel[["n_deaths"]], death_panel[["m"]]))
  intensity_summary <- list(
    intensity_kappa_min = min(death_calibration$sd$kappa),
    intensity_kappa_max = max(death_calibration$sd$kappa),
    intensity_m_min = min(multipliers$m),
    intensity_m_max = max(multipliers$m)
  )
  count_totals[["n_deaths_total"]] <- sum(as.numeric(death_panel[["n_deaths"]]))
  set(death_panel, j = c("m", "cause_6"), value = NULL)
  dst <- draw_path(PERTURBED_INPUTS[[1]])
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  saveRDS(as.data.frame(death_panel), dst)
  rm(death_panel, multipliers)
  invisible(gc(verbose = FALSE))

  # Births: one multiplier per year and maternal race.
  birth_panel <- as.data.table(readRDS(repository_path(PERTURBED_INPUTS[[2]])))
  s_birth <- calibrate_birth_intensity(birth_panel)
  use_rng_stream(streams$birth)
  birth_years <- sort(unique(birth_panel$year))
  birth_multipliers <- CJ(year = birth_years, mother_race_eth = RACES_R3)
  birth_multipliers[, m := exp(rnorm(.N, mean = -s_birth^2 / 2, sd = s_birth))]
  birth_panel[birth_multipliers, m := i.m, on = .(year, mother_race_eth)]
  birth_panel[is.na(m), m := 1]
  set(birth_panel, j = "n_births",
      value = poisson_resample_counts(birth_panel[["n_births"]], birth_panel[["m"]]))
  intensity_summary$intensity_s_birth <- s_birth
  count_totals[["n_births_total"]] <- sum(as.numeric(birth_panel[["n_births"]]))
  set(birth_panel, j = "m", value = NULL)
  dst <- draw_path(PERTURBED_INPUTS[[2]])
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  saveRDS(as.data.frame(birth_panel), dst)
  rm(birth_panel, birth_multipliers)
  invisible(gc(verbose = FALSE))

  # Apply a common replicate-weight draw to all population products.
  replicate_panel <- as.data.table(readRDS(replicate_panel_path))
  use_rng_stream(streams$acs)
  epsilon <- draw_acs_epsilon(replicate_panel$year)
  replicate_panel[, theta_factor := acs_draw_factors(replicate_panel, epsilon)]

  main_panel <- as.data.table(readRDS(main_panel_path))
  join_key <- c("source", "year", "sex", "age_years", "race_eth", "educ2")
  main_panel[replicate_panel, theta_factor := i.theta_factor, on = join_key]
  if (anyNA(main_panel$theta_factor)) {
    stop("A population cell has no theta* multiplier; 04 and 13 disagree on the cell grid.", call. = FALSE)
  }
  main_panel[, n_persons_weighted := n_persons_weighted * theta_factor]
  main_panel[, theta_factor := NULL]

  perturbed_main <- as.data.frame(main_panel)
  population_paths <- lapply(PERTURBED_INPUTS[3:5], draw_path)
  lapply(population_paths, function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  })
  saveRDS(perturbed_main, population_paths[[1]])
  saveRDS(derive_death_denominator(perturbed_main), population_paths[[2]])
  saveRDS(derive_birth_denominator(perturbed_main), population_paths[[3]])

  c(
    list(
      n_epsilon_years = length(epsilon),
      n_deaths_total = count_totals[["n_deaths_total"]],
      n_births_total = count_totals[["n_births_total"]],
      theta_factor_min = min(replicate_panel$theta_factor),
      theta_factor_max = max(replicate_panel$theta_factor)
    ),
    intensity_summary
  )
}

DRAW_STREAMS <- if (identical(DRAW_MODE, "random")) {
  draw_component_streams(seed_ledger[[DRAW_INDEX]])
} else {
  NULL
}

draw_prep_summary <- if (identical(DRAW_MODE, "identity")) {
  write_identity_inputs()
} else {
  write_random_inputs(DRAW_STREAMS)
}

# ---- Run the analysis pipeline ----

log_dir <- file.path(DRAW_ROOT, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

for (stage in STAGES) {
  script_path <- analysis_path("02_codes", stage)
  stage_stem <- file.path(log_dir, sub("\\.r$", "", stage))
  log_out <- paste0(stage_stem, ".out")
  log_err <- paste0(stage_stem, ".err")
  previous_wd <- setwd(analysis_root)
  rc <- tryCatch(
    system2(
      RSCRIPT,
      args = c("--vanilla", shQuote(script_path)),
      stdout = log_out,
      stderr = log_err
    ),
    finally = setwd(previous_wd)
  )
  # quit(), not stop(): the record already names the stage and exit code, and
  # raising here would re-enter the handler above and overwrite both.
  if (rc != 0L) {
    message_text <- sprintf("Draw %s stage %s exited %d (see %s).", DRAW_ID, stage, rc, log_dir)
    record_failed_draw(stage, rc, message_text, c(log_out, log_err))
    message(message_text)
    quit(status = 1L, save = "no")
  }
}

# ---- Endpoint shard ----

build_shard_tables <- function() {
  endpoint_dir <- analysis_path("01_data", "02_processed", "edu_parent_loss")
  list(
    incidence_total = order_total(build_total_incidence(
      readRDS(file.path(endpoint_dir, "02_parent_loss_by_cause_atom_annual.rds")),
      EDUC_TABLE_LEVELS
    )),
    incidence_cause_atom = order_cause(
      build_cause_incidence(
        readRDS(file.path(endpoint_dir, "02_parent_loss_by_cause_atom_annual.rds")),
        EDUC_TABLE_LEVELS, "cause_atom", CAUSE_ATOM_LEVELS
      ),
      "cause_atom", CAUSE_ATOM_LEVELS
    ),
    meanage_total = order_total(build_total_mean_age(
      readRDS(file.path(endpoint_dir, "03_parent_loss_by_cause_atom_child_age.rds")),
      EDUC_TABLE_LEVELS
    )),
    meanage_cause_atom = order_cause(
      build_cause_mean_age(
        readRDS(file.path(endpoint_dir, "03_parent_loss_by_cause_atom_child_age.rds")),
        EDUC_TABLE_LEVELS, "cause_atom", CAUSE_ATOM_LEVELS
      ),
      "cause_atom", CAUSE_ATOM_LEVELS
    )
  )
}

tables <- build_shard_tables()
errs <- check_shard_tables(tables)
if (length(errs) > 0L) {
  message_text <- sprintf(
    "Draw %s shard failed validation:\n%s", DRAW_ID, paste0("  - ", errs, collapse = "\n")
  )
  record_failed_draw("shard", 0L, message_text)
  message(message_text)
  quit(status = 1L, save = "no")
}

shard <- list(
  meta = c(
    list(
      draw_id = DRAW_ID,
      batch_id = BATCH_ID,
      mode = DRAW_MODE,
      master_seed = DRAW_MASTER_SEED,
      ledger_index = DRAW_INDEX,
      seed_streams = DRAW_STREAMS,
      code_version = CODE_VERSION,
      status = "complete"
    ),
    draw_prep_summary
  ),
  tables = tables
)

# Write to scratch inside the batch root first, so a job killed mid-write
# cannot leave a truncated shard at the final path.
shard_tmp_path <- tempfile(paste0(shard_stem, "_"), shard_dir, ".tmp")
saveRDS(shard, shard_tmp_path)
invisible(file.rename(shard_tmp_path, shard_final_path))
invisible(unlink(shard_fail_path))

cat(sprintf("Draw %s (%s) complete\n", DRAW_ID, DRAW_MODE))
