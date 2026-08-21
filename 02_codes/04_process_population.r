# ==============================================================================
# 04_process_population.r
#
# Build the ACS population panel and the mortality and fertility denominators.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(here)
  library(data.table)
  library(ipumsr)
})

# ---- ACS processing settings ----

# The ACS extract is a multi-GB fixed-width file, so it is streamed in chunks.
# The chunk size can be adjusted for the available memory.
chunk_size <- {
  raw_value <- Sys.getenv("IPUMS_POP_CHUNK_SIZE", unset = "100000")
  parsed <- suppressWarnings(as.integer(raw_value))
  if (is.na(parsed) || parsed <= 0L) 100000L else parsed
}

# ACS is the sole population source; the national series begins in 2000.
ACS_MIN_YEAR <- 2000L
ACS_MAX_YEAR <- 2024L

# IPUMS-USA variables carrying detailed education and the person weight.
EDUC_VAR <- "educd"
WEIGHT_VAR <- "perwt"

# ---- ACS files ----

here::i_am("02_codes/04_process_population.r")

source(here::here("02_codes", "00_functions.r"))

raw_dir <- here("01_data", "01_raw", "acs")
proc_dir <- here("01_data", "02_processed", "population")

main_panel_rds <- file.path(proc_dir, "01_pop_main_panel.rds")
death_denominator_rds <- file.path(proc_dir, "02_pop_denominator_death.rds")
birth_denominator_rds <- file.path(proc_dir, "03_pop_denominator_birth.rds")

dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)

# ---- ACS reducers ----
# The ACS person-record recodes (sex, age, race/ethnicity, EDUCD education) are
# shared with 13_replicate_weights.r and live in 00_functions.r, so the
# point-weight and replicate panels classify every person identically.

parse_raw_value <- function(x) {
  raw_chr <- trimws(as.character(x))
  raw_chr[raw_chr %in% c("", "NA", "N/A")] <- NA_character_
  raw_chr
}

select_vars_to_read <- function(ddi) {
  available_vars <- tolower(ipumsr::ipums_var_info(ddi)$var_name)

  vars_to_read <- c("year", "sex", "age", "race", "hispan", EDUC_VAR, WEIGHT_VAR)
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

pop_main_panel_columns <- c(
  "source",
  "year",
  "sex_raw",
  "sex",
  "age_raw",
  "age_years",
  "race_raw",
  "hispanic_raw",
  "race_eth",
  "educ_raw",
  "educ2",
  "weight_var",
  "n_persons_unweighted",
  "n_persons_weighted"
)

empty_pop_panel <- function() {
  out <- as.data.frame(setNames(replicate(length(pop_main_panel_columns), logical(0), simplify = FALSE), pop_main_panel_columns))
  character_cols <- c("source", "sex_raw", "sex", "age_raw", "race_raw", "hispanic_raw", "race_eth", "educ_raw", "educ2", "weight_var")
  integer_cols <- c("year", "age_years", "n_persons_unweighted")
  numeric_cols <- c("n_persons_weighted")

  for (col_name in character_cols) out[[col_name]] <- character(0)
  for (col_name in integer_cols) out[[col_name]] <- integer(0)
  for (col_name in numeric_cols) out[[col_name]] <- numeric(0)
  out
}

align_pop_panel_columns <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  missing_cols <- setdiff(pop_main_panel_columns, names(df))
  if (length(missing_cols)) {
    for (col_name in missing_cols) df[[col_name]] <- NA
  }

  df <- df[, pop_main_panel_columns, drop = FALSE]
  df$year <- safe_integer(df$year)
  df$age_years <- safe_integer(df$age_years)
  df$n_persons_unweighted <- safe_integer(df$n_persons_unweighted)
  df$n_persons_weighted <- safe_numeric(df$n_persons_weighted)
  df
}

aggregate_pop_panel <- function(panel_df) {
  if (!nrow(panel_df)) return(empty_pop_panel())

  dt <- data.table::as.data.table(align_pop_panel_columns(panel_df))
  group_cols <- setdiff(pop_main_panel_columns, c("n_persons_unweighted", "n_persons_weighted"))

  dt[is.na(n_persons_unweighted), n_persons_unweighted := 0L]
  dt[is.na(n_persons_weighted), n_persons_weighted := 0]

  out <- dt[, .(
    n_persons_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
    n_persons_weighted = sum(n_persons_weighted, na.rm = TRUE)
  ), by = group_cols]

  data.table::setorderv(out, group_cols)
  as.data.frame(out)
}

# Harmonize one raw chunk into the analysis variables, one row per person.
prepare_population_chunk <- function(chunk_df) {
  panel <- data.frame(
    source = "acs",
    year = safe_integer(chunk_df$year),
    sex_raw = parse_raw_value(chunk_df$sex),
    sex = harmonize_acs_sex(chunk_df$sex),
    age_raw = parse_raw_value(chunk_df$age),
    age_years = harmonize_acs_age(chunk_df$age),
    race_raw = parse_raw_value(chunk_df$race),
    hispanic_raw = parse_raw_value(chunk_df$hispan),
    race_eth = harmonize_acs_race_eth(chunk_df$race, chunk_df$hispan),
    educ_raw = parse_raw_value(chunk_df[[EDUC_VAR]]),
    educ2 = recode_educ2_from_educd(chunk_df[[EDUC_VAR]]),
    weight_var = toupper(WEIGHT_VAR),
    n_persons_unweighted = 1L,
    n_persons_weighted = safe_numeric(chunk_df[[WEIGHT_VAR]]),
    stringsAsFactors = FALSE
  )

  panel <- panel[
    !is.na(panel$year) &
      panel$year >= ACS_MIN_YEAR &
      panel$year <= ACS_MAX_YEAR,
    ,
    drop = FALSE
  ]

  align_pop_panel_columns(panel)
}

# The two denominator derivations are shared with the ACS theta* draw prep and
# live in 00_functions.r and are also used to construct each ACS draw.

# ---- Read and harmonize the ACS extract ----

acs_files <- find_acs_files(raw_dir)
ddi <- ipumsr::read_ipums_ddi(acs_files$ddi_file, lower_vars = TRUE)
vars_to_read <- select_vars_to_read(ddi)

aggregated_chunks <- list()
chunk_counter <- 0L

callback <- ipumsr::IpumsSideEffectCallback$new(function(chunk, index) {
  panel_chunk <- prepare_population_chunk(chunk)

  if (!nrow(panel_chunk)) {
    return(invisible(NULL))
  }

  chunk_counter <<- chunk_counter + 1L
  aggregated_chunks[[chunk_counter]] <<- aggregate_pop_panel(panel_chunk)
  invisible(NULL)
})

ipumsr::read_ipums_micro_chunked(
  ddi = ddi,
  callback = callback,
  chunk_size = chunk_size,
  vars = all_of(vars_to_read),
  data_file = acs_files$data_file,
  verbose = TRUE,
  var_attrs = c("val_labels")
)

main_panel <- if (!length(aggregated_chunks)) {
  empty_pop_panel()
} else {
  aggregate_pop_panel(data.table::rbindlist(aggregated_chunks, use.names = TRUE, fill = TRUE))
}
if (!setequal(unique(main_panel$year), ACS_MIN_YEAR:ACS_MAX_YEAR)) {
  stop(
    sprintf(
      "The ACS panel must contain every year from %d through %d.",
      ACS_MIN_YEAR, ACS_MAX_YEAR
    ),
    call. = FALSE
  )
}

# ---- Derive denominators ----

death_denominator <- derive_death_denominator(main_panel)
birth_denominator <- derive_birth_denominator(main_panel)

# ---- Population panels ----

saveRDS(main_panel, file = main_panel_rds)
saveRDS(death_denominator, file = death_denominator_rds)
saveRDS(birth_denominator, file = birth_denominator_rds)
