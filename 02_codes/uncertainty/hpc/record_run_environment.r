# ==============================================================================
# record_run_environment.r
#
# Create or verify the small set of settings shared by one uncertainty batch.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/uncertainty/hpc/record_run_environment.r")
source(here::here("02_codes", "00_functions.r"))

# ---- Batch settings ----

BATCH_ROOT <- Sys.getenv("EDU_PARENT_LOSS_SHARD_ROOT", unset = "")
BATCH_ID <- Sys.getenv("EDU_PARENT_LOSS_BATCH_ID", unset = "")
BATCH_MODE <- Sys.getenv("EDU_PARENT_LOSS_BATCH_MODE", unset = "")
M_EXPECTED <- suppressWarnings(as.integer(
  Sys.getenv("EDU_PARENT_LOSS_M_EXPECTED", unset = "1000")
))

if (!BATCH_MODE %in% c("random", "smoke")) {
  stop("EDU_PARENT_LOSS_BATCH_MODE must be random or smoke.", call. = FALSE)
}

dir.create(BATCH_ROOT, recursive = TRUE, showWarnings = FALSE)
for (directory in c("shards", "failures", "logs", "accounting")) {
  dir.create(file.path(BATCH_ROOT, directory), recursive = TRUE,
             showWarnings = FALSE)
}

batch_info_path <- file.path(BATCH_ROOT, "batch_info.rds")
batch_text_path <- file.path(BATCH_ROOT, "batch_info.txt")
ledger_path <- file.path(BATCH_ROOT, "seed_ledger.rds")
environment_text_path <- file.path(BATCH_ROOT, "run_environment.txt")

# ---- Environment and fixed inputs ----

PACKAGE_NAMES <- c(
  "data.table", "here", "DemoKin", "Matrix", "tidyselect"
)

FIXED_INPUTS <- list(
  c("01_data", "02_processed", "deaths", "01_death_main_panel.rds"),
  c("01_data", "02_processed", "births", "01_birth_main_panel.rds"),
  c("01_data", "02_processed", "population", "01_pop_main_panel.rds"),
  c("01_data", "02_processed", "population", "02_pop_denominator_death.rds"),
  c("01_data", "02_processed", "population", "03_pop_denominator_birth.rds"),
  c("01_data", "02_processed", "population", "04_pop_replicate_panel.rds"),
  c("01_data", "02_processed", "edu_parent_loss",
    "point_estimate_incidence_total.rds"),
  c("01_data", "02_processed", "edu_parent_loss",
    "point_estimate_incidence_cause6.rds"),
  c("01_data", "02_processed", "edu_parent_loss",
    "point_estimate_meanage_total.rds"),
  c("01_data", "02_processed", "edu_parent_loss",
    "point_estimate_meanage_cause6.rds")
)

package_record <- function(package_name) {
  data.table(
    package = package_name,
    version = as.character(utils::packageVersion(package_name))
  )
}

input_record <- function(parts) {
  path <- do.call(here::here, as.list(parts))
  data.table(
    input = paste(parts, collapse = "/"),
    bytes = file.size(path)
  )
}

packages <- rbindlist(lapply(PACKAGE_NAMES, package_record))
if (!exists(
  "kin_multi_stage_time_variant_2sex",
  envir = asNamespace("DemoKin"),
  inherits = FALSE
)) {
  stop(
    "The installed DemoKin does not export kin_multi_stage_time_variant_2sex().",
    call. = FALSE
  )
}
packages <- as.data.frame(packages)

inputs <- rbindlist(lapply(FIXED_INPUTS, input_record))
if (anyNA(inputs$bytes)) {
  stop(
    "Required batch inputs are missing:\n  ",
    paste(inputs[is.na(bytes), input], collapse = "\n  "),
    call. = FALSE
  )
}
inputs <- as.data.frame(inputs)

code_version <- current_commit()
if (is.na(code_version) || !nzchar(code_version)) {
  stop(
    "Set EDU_PARENT_LOSS_CODE_VERSION, or run from a Git checkout.",
    call. = FALSE
  )
}

# ---- Seeds ----

build_seed_ledger <- function(n_draws, master_seed = DRAW_MASTER_SEED) {
  set.seed(master_seed, kind = DRAW_RNG_KIND)
  ledger <- vector("list", n_draws)
  stream <- get(".Random.seed", envir = .GlobalEnv)
  for (draw_index in seq_len(n_draws)) {
    ledger[[draw_index]] <- stream
    stream <- parallel::nextRNGStream(stream)
  }
  ledger
}

save_rds <- function(object, path) {
  temporary <- tempfile(paste0(basename(path), "_"), dirname(path), ".tmp")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("Could not write ", path, call. = FALSE)
  }
}

ledger <- build_seed_ledger(M_EXPECTED)
if (file.exists(ledger_path)) {
  if (!identical(readRDS(ledger_path), ledger)) {
    stop("The existing seed ledger differs from this batch.", call. = FALSE)
  }
} else {
  save_rds(ledger, ledger_path)
}

# ---- Batch record ----

candidate <- list(
  batch_id = BATCH_ID,
  batch_mode = BATCH_MODE,
  m_expected = M_EXPECTED,
  master_seed = DRAW_MASTER_SEED,
  rng_kind = DRAW_RNG_KIND,
  code_version = code_version,
  r_version = as.character(getRversion()),
  packages = packages,
  inputs = inputs,
  seed_ledger_file = basename(ledger_path)
)

if (file.exists(batch_info_path)) {
  batch_info <- readRDS(batch_info_path)
  if (!identical(batch_info[names(candidate)], candidate)) {
    stop(
      "The existing batch was created with different code, inputs or settings.",
      call. = FALSE
    )
  }
  message("Verified batch: ", BATCH_ID)
} else {
  batch_info <- c(
    candidate,
    list(created_at = format(Sys.time(), tz = "UTC", usetz = TRUE))
  )
  save_rds(batch_info, batch_info_path)
  message("Created batch: ", BATCH_ID)
}

batch_lines <- c(
  paste0("batch_id:      ", batch_info$batch_id),
  paste0("batch_mode:    ", batch_info$batch_mode),
  paste0("created_at:    ", batch_info$created_at),
  paste0("code_version:  ", batch_info$code_version),
  paste0("M:             ", batch_info$m_expected),
  paste0("master_seed:   ", batch_info$master_seed),
  paste0("RNG:           ", batch_info$rng_kind),
  "",
  "Packages:",
  utils::capture.output(print(batch_info$packages)),
  "",
  "Fixed inputs (path and bytes):",
  utils::capture.output(print(batch_info$inputs))
)
writeLines(batch_lines, batch_text_path, useBytes = TRUE)

environment_lines <- c(
  paste0("recorded_at: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("R:           ", R.version.string),
  paste0("platform:    ", R.version$platform),
  paste0("code:        ", code_version),
  "",
  utils::capture.output(utils::sessionInfo())
)
writeLines(environment_lines, environment_text_path, useBytes = TRUE)

cat(sprintf("Batch %s ready: %d draws\n", BATCH_ID, M_EXPECTED))
