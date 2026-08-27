# ==============================================================================
# 03_process_deaths.r
#
# Harmonize NCHS mortality records across coding eras and build the
# cause-classified resident-death count panel.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

here::i_am("02_codes/03_process_deaths.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Mortality processing settings ----

# The multi-GB files are read in chunks; the chunk size may be adjusted for the
# available memory.
chunk_size <- {
  raw_value <- Sys.getenv("NCHS_MORT_CHUNK_SIZE", unset = "100000")
  parsed <- suppressWarnings(as.integer(raw_value))

  if (is.na(parsed) || parsed <= 0L) {
    100000L
  } else {
    parsed
  }
}

save_minimal_extracts <- {
  raw_value <- tolower(trimws(Sys.getenv("NCHS_SAVE_DEATH_MINIMAL_EXTRACTS", unset = "false")))
  raw_value %in% c("true", "t", "1", "yes", "y")
}

# Formal runs rebuild every year. Set false only to resume an interrupted run.
force_reprocess <- {
  raw_value <- tolower(trimws(Sys.getenv("NCHS_DEATH_FORCE_REPROCESS", unset = "true")))
  raw_value %in% c("true", "t", "1", "yes", "y")
}

# ---- Mortality files ----

raw_dir <- here("01_data", "01_raw", "nchs", "mortality")
proc_dir <- here("01_data", "02_processed", "deaths")
extract_dir <- file.path(proc_dir, "minimal_extracts")
per_year_dir <- file.path(proc_dir, "per_year")

main_panel_rds_path <- file.path(proc_dir, "01_death_main_panel.rds")

# ---- Mortality readers and recodes ----

extract_year_from_file_name <- function(file_name) {
  file_name <- basename(file_name)
  year_text <- sub(".*?(\\d{4}).*", "\\1", file_name)

  if (!grepl("^\\d{4}$", year_text)) {
    return(NA_integer_)
  }

  as.integer(year_text)
}

get_zip_member_name <- function(zip_path) {
  zip_listing <- utils::unzip(zip_path, list = TRUE)
  member_names <- zip_listing$Name[!grepl("/$", zip_listing$Name)]

  if (length(member_names) == 0L) {
    stop("No file member was found inside zip archive: ", zip_path, call. = FALSE)
  }

  member_names[1]
}

escape_powershell_single_quoted <- function(x) {
  gsub("'", "''", x, fixed = TRUE)
}

extract_zip_member_with_powershell <- function(zip_path, member_name, temp_dir) {
  ps_exe <- Sys.which("pwsh")

  if (identical(ps_exe, "")) {
    ps_exe <- Sys.which("powershell")
  }

  if (identical(ps_exe, "")) {
    stop(
      "PowerShell was not found, so zip extraction fallback is unavailable.",
      call. = FALSE
    )
  }

  zip_path <- normalizePath(zip_path, winslash = "\\", mustWork = TRUE)
  temp_dir <- normalizePath(temp_dir, winslash = "\\", mustWork = TRUE)
  extracted_path <- file.path(temp_dir, member_name)

  dir.create(dirname(extracted_path), recursive = TRUE, showWarnings = FALSE)

  ps_command <- paste0(
    "Add-Type -AssemblyName System.IO.Compression.FileSystem; ",
    "$zip = [System.IO.Compression.ZipFile]::OpenRead('",
    escape_powershell_single_quoted(zip_path),
    "'); ",
    "$entry = $zip.GetEntry('",
    escape_powershell_single_quoted(member_name),
    "'); ",
    "if ($null -eq $entry) { $zip.Dispose(); throw 'Zip member not found.' }; ",
    "$dest = '",
    escape_powershell_single_quoted(extracted_path),
    "'; ",
    "$destDir = [System.IO.Path]::GetDirectoryName($dest); ",
    "if (-not [string]::IsNullOrEmpty($destDir)) { ",
    "[System.IO.Directory]::CreateDirectory($destDir) | Out-Null ",
    "}; ",
    "[System.IO.Compression.ZipFileExtensions]::ExtractToFile($entry, $dest, $true); ",
    "$zip.Dispose()"
  )

  ps_output <- suppressWarnings(
    system2(
      ps_exe,
      c("-NoProfile", "-Command", ps_command),
      stdout = TRUE,
      stderr = TRUE
    )
  )

  ps_status <- attr(ps_output, "status")

  if (!is.null(ps_status) && ps_status != 0L) {
    stop(
      "PowerShell zip extraction fallback failed for ",
      basename(zip_path),
      ": ",
      paste(ps_output, collapse = "\n"),
      call. = FALSE
    )
  }

  if (!file.exists(extracted_path)) {
    stop(
      "PowerShell zip extraction fallback did not create the expected file: ",
      extracted_path,
      call. = FALSE
    )
  }

  extracted_path
}

extract_zip_member_to_temp <- function(zip_path) {
  member_name <- get_zip_member_name(zip_path)
  temp_dir <- tempfile(pattern = "nchs_mort_")

  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  # `utils::unzip()` intermittently fails on the very large 2020 public-use
  # archive; fall back to PowerShell/.NET extraction.
  extracted_path <- tryCatch(
    suppressWarnings(
      utils::unzip(
        zipfile = zip_path,
        files = member_name,
        exdir = temp_dir,
        overwrite = TRUE
      )
    ),
    error = function(e) character(0)
  )

  extracted_path <- extracted_path[file.exists(extracted_path)]

  if (length(extracted_path) == 0L) {
    extracted_path <- extract_zip_member_with_powershell(
      zip_path = zip_path,
      member_name = member_name,
      temp_dir = temp_dir
    )
  } else {
    extracted_path <- extracted_path[[1]]
  }

  list(
    member_name = member_name,
    extracted_path = extracted_path,
    temp_dir = temp_dir
  )
}

trim_fixed_width_field <- function(x) {
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

pad_record_chunk <- function(record_chunk, expected_length) {
  # Fixed-width mortality files are byte-positioned. Some older files contain
  # non-standard bytes that break character-counting but not byte-counting, so
  # width checks use bytes.
  record_lengths <- nchar(record_chunk, type = "bytes", allowNA = FALSE, keepNA = FALSE)
  short_records <- record_lengths < expected_length

  if (any(short_records)) {
    record_chunk[short_records] <- paste0(
      record_chunk[short_records],
      strrep(" ", expected_length - record_lengths[short_records])
    )
  }

  substr(record_chunk, 1L, expected_length)
}

extract_fields_from_chunk <- function(record_chunk, fields, expected_length) {
  # Slice only the fields declared in layout_registry from the fixed-width text.
  record_chunk <- pad_record_chunk(record_chunk, expected_length)

  extracted <- lapply(
    seq_len(nrow(fields)),
    function(i) {
      trim_fixed_width_field(
        substr(record_chunk, fields$start[i], fields$end[i])
      )
    }
  )

  as.data.frame(
    setNames(extracted, fields$field_name),
    stringsAsFactors = FALSE
  )
}

harmonize_sex <- function(sex_raw) {
  out <- rep(NA_character_, length(sex_raw))
  out[sex_raw %in% c("1", "M", "Male")] <- "male"
  out[sex_raw %in% c("2", "F", "Female")] <- "female"
  out[is.na(out)] <- sex_raw[is.na(out)]
  out
}

convert_age_raw_1990 <- function(age_raw) {
  # NCHS 1990-2002 detail age coding (handled as in Schlüter et al., 2024):
  #   000-199 = completed years; 2/3/4/5/6xx = sub-year units -> 0;
  #   299/399/499/599/699/999 = missing/special codes.
  age_num <- suppressWarnings(as.integer(age_raw))
  age_years <- rep(NA_integer_, length(age_num))

  explicit_missing <- !is.na(age_num) & age_num %in% c(299L, 399L, 499L, 599L, 699L, 999L)
  whole_year_ages <- !is.na(age_num) & age_num >= 0L & age_num <= 199L
  subyear_ages <- !is.na(age_num) & age_num > 199L & !explicit_missing

  age_years[whole_year_ages] <- age_num[whole_year_ages]
  age_years[subyear_ages] <- 0L

  age_years
}

convert_age_raw_2003 <- function(age_raw) {
  # NCHS 2003+ detail age coding (handled as in Schlüter et al., 2024):
  #   1xxx = completed years (leading unit digit 1); 2/4/5/6xxx = sub-year -> 0;
  #   1999/2999/4999/5999/6999/9999 = missing/special codes.
  age_num <- suppressWarnings(as.integer(age_raw))
  age_years <- rep(NA_integer_, length(age_num))

  explicit_missing <- !is.na(age_num) & age_num %in% c(1999L, 2999L, 4999L, 5999L, 6999L, 9999L)
  # Whole-year ages use leading unit digit 1 (1000-1998 -> 0-998 completed
  # years). The lower bound guards against malformed sub-1000 codes producing
  # negative ages.
  whole_year_ages <- !is.na(age_num) & age_num >= 1000L & age_num < 1999L
  subyear_ages <- !is.na(age_num) & age_num > 1999L & !explicit_missing

  age_years[whole_year_ages] <- age_num[whole_year_ages] - 1000L
  age_years[subyear_ages] <- 0L

  age_years
}

harmonize_race_bridged_1989_revision <- function(race_raw) {
  # Bridged single-race buckets (white/black/aian/api) as in Schlüter et al.,
  # 2024: race codes 4-7 and the expanded Asian/Pacific Islander codes
  # (18/28/.../78) map to api; unmapped codes stay NA.
  race_num <- suppressWarnings(as.integer(race_raw))
  out <- rep(NA_character_, length(race_num))

  out[race_num == 1L] <- "white"
  out[race_num == 2L] <- "black"
  out[race_num == 3L] <- "aian"
  out[!is.na(race_num) & race_num >= 4L & race_num <= 7L] <- "api"
  out[race_num %in% c(18L, 28L, 38L, 48L, 58L, 68L, 78L)] <- "api"

  out
}

harmonize_race_eth_1989_revision <- function(hispanic_race_recode_raw, race_raw) {
  # Hispanic origin recode 1-5 and 9 (not stated) -> hispanic, as in Schlüter
  # et al., 2024; otherwise keep the bridged race bucket. Note this folds
  # "not stated" (9) into hispanic for 1990-2020.
  out <- harmonize_race_bridged_1989_revision(race_raw)

  out[hispanic_race_recode_raw %in% c("1", "2", "3", "4", "5", "9")] <- "hispanic"

  out
}

harmonize_race_bridged_2022plus <- function(race_recode_6_raw) {
  # 2022+ public-use files carry single-race recodes instead of bridged race.
  # Collapse Race Recode 6 to white/black/aian/api (Asian and Pacific Islander
  # merged, as in Villaveces et al., 2025); "more than one race" stays NA.
  race_num <- suppressWarnings(as.integer(race_recode_6_raw))
  out <- rep(NA_character_, length(race_num))

  out[race_num == 1L] <- "white"
  out[race_num == 2L] <- "black"
  out[race_num == 3L] <- "aian"
  out[race_num %in% c(4L, 5L)] <- "api"

  out
}

harmonize_race_eth_2022plus <- function(hispanic_race_recode_2022_raw, race_recode_6_raw) {
  # 2022+ Hispanic/race recode: 01-07 and 14 -> hispanic; 08-12 -> non-Hispanic
  # white/black/aian/api; unknown / more-than-one-race stays NA (single-race
  # era, aligned with Villaveces et al., 2025).
  out <- harmonize_race_bridged_2022plus(race_recode_6_raw)

  out[hispanic_race_recode_2022_raw %in% sprintf("%02d", 1:7)] <- "hispanic"
  out[hispanic_race_recode_2022_raw == "14"] <- "hispanic"
  out[hispanic_race_recode_2022_raw == "08"] <- "white"
  out[hispanic_race_recode_2022_raw == "09"] <- "black"
  out[hispanic_race_recode_2022_raw == "10"] <- "aian"
  out[hispanic_race_recode_2022_raw %in% c("11", "12")] <- "api"

  out
}

harmonize_race_bridged_2021 <- function(race_recode_40_raw) {
  # The 2021 public-use mortality file is a transition year: the later Race
  # Recode 6 field used in 2022+ is not populated in the same way, but Race
  # Recode 40 is available. Collapse Race Recode 40 to the project's broad
  # race buckets while leaving multi-racial / other groups as missing.
  race_num <- suppressWarnings(as.integer(race_recode_40_raw))
  out <- rep(NA_character_, length(race_num))

  out[race_num == 1L] <- "white"
  out[race_num == 2L] <- "black"
  out[race_num == 3L] <- "aian"
  out[!is.na(race_num) & race_num >= 4L & race_num <= 14L] <- "api"

  out
}

harmonize_race_eth_2021 <- function(hispanic_raw, race_recode_40_raw) {
  # 2021 still carries the detailed Hispanic-origin recode (the 2-digit
  # Hispanic/race recode used in 2022+ is not yet available). NCHS Hispanic
  # origin codes: 100-199 = non-Hispanic (keep the bridged race); 200-299 =
  # Hispanic; 996-999 = uncodable / deferred / unknown / not stated. Origin
  # that is not stated is set to NA because race/ethnicity
  # cannot be determined; this keeps 2021 consistent with the 2022+ single-race
  # handling (where the unknown code maps to NA).
  out <- harmonize_race_bridged_2021(race_recode_40_raw)
  hispanic_num <- suppressWarnings(as.integer(hispanic_raw))

  out[!is.na(hispanic_num) & hispanic_num >= 200L & hispanic_num <= 299L] <- "hispanic"
  out[!is.na(hispanic_num) & hispanic_num >= 996L] <- NA_character_

  out
}

harmonize_education_2003 <- function(educ_1989_raw, educ_2003_raw, educ_flag_raw) {
  n <- length(educ_flag_raw)
  educ_raw_selected <- rep(NA_character_, n)
  educ_system <- rep(NA_character_, n)
  educ2 <- rep("Unknown", n)

  use_1989 <- educ_flag_raw == "0"
  use_2003 <- educ_flag_raw == "1"
  no_item <- educ_flag_raw == "2"

  educ_raw_selected[use_1989] <- educ_1989_raw[use_1989]
  educ_system[use_1989] <- "1989_revision"
  educ2[use_1989] <- recode_educ2_from_1989(educ_1989_raw[use_1989])

  educ_raw_selected[use_2003] <- educ_2003_raw[use_2003]
  educ_system[use_2003] <- "2003_revision"
  educ2[use_2003] <- recode_educ2_from_2003(educ_2003_raw[use_2003])

  educ_system[no_item] <- "no_education_item"

  fallback_1989 <- is.na(educ_system) & !is.na(educ_1989_raw) & is.na(educ_2003_raw)
  fallback_2003 <- is.na(educ_system) & !is.na(educ_2003_raw)
  fallback_none <- is.na(educ_system)

  educ_raw_selected[fallback_1989] <- educ_1989_raw[fallback_1989]
  educ_system[fallback_1989] <- "1989_revision"
  educ2[fallback_1989] <- recode_educ2_from_1989(educ_1989_raw[fallback_1989])

  educ_raw_selected[fallback_2003] <- educ_2003_raw[fallback_2003]
  educ_system[fallback_2003] <- "2003_revision"
  educ2[fallback_2003] <- recode_educ2_from_2003(educ_2003_raw[fallback_2003])

  educ_system[fallback_none] <- "unknown_revision"

  list(
    educ_raw_selected = educ_raw_selected,
    educ_system = educ_system,
    educ2 = educ2
  )
}

get_optional_field <- function(raw_df, field_name, default = NA_character_) {
  if (field_name %in% names(raw_df)) {
    raw_df[[field_name]]
  } else {
    rep.int(default, nrow(raw_df))
  }
}

derive_mortality_fields <- function(raw_df, year, file_name, member_name) {
  # Each harmonized field keeps its raw source column.
  if (!year %in% 1990L:2024L) {
    stop("No field-derivation logic has been implemented for year ", year, call. = FALSE)
  }

  drug_death <- rep.int(NA_integer_, nrow(raw_df))
  opioid_death <- rep.int(NA_integer_, nrow(raw_df))
  firearm_death <- rep.int(NA_integer_, nrow(raw_df))
  cause_atom <- rep.int(NA_character_, nrow(raw_df))
  age_years <- rep.int(NA_integer_, nrow(raw_df))

  educ_1989_raw <- get_optional_field(raw_df, "educ_1989_raw")
  educ_1989_recode_raw <- get_optional_field(raw_df, "educ_1989_recode_raw")
  educ_2003_raw <- get_optional_field(raw_df, "educ_2003_raw")
  educ_flag_raw <- get_optional_field(raw_df, "educ_flag_raw")
  age_recode_12_raw <- get_optional_field(raw_df, "age_recode_12_raw")
  age_recode_27_raw <- get_optional_field(raw_df, "age_recode_27_raw")
  race_recode_2_raw <- get_optional_field(raw_df, "race_recode_2_raw")
  race_recode_5_raw <- get_optional_field(raw_df, "race_recode_5_raw")
  race_recode_6_raw <- get_optional_field(raw_df, "race_recode_6_raw")
  race_recode_40_raw <- get_optional_field(raw_df, "race_recode_40_raw")
  race_bridged_flag_raw <- get_optional_field(raw_df, "race_bridged_flag_raw")
  race_imputation_flag_raw <- get_optional_field(raw_df, "race_imputation_flag_raw")
  hispanic_race_recode_2022_raw <- get_optional_field(raw_df, "hispanic_race_recode_2022_raw")
  race_raw <- get_optional_field(raw_df, "race_raw")
  race_recode_3_raw <- get_optional_field(raw_df, "race_recode_3_raw")
  hispanic_race_recode_raw <- get_optional_field(raw_df, "hispanic_race_recode_raw")
  race_bridged <- harmonize_race_bridged_1989_revision(race_raw)
  race_eth <- harmonize_race_eth_1989_revision(
    hispanic_race_recode_raw = hispanic_race_recode_raw,
    race_raw = race_raw
  )

  educ_meta <- list(
    educ_raw_selected = educ_1989_raw,
    educ_system = rep.int("1989_revision", nrow(raw_df)),
    educ2 = recode_educ2_from_1989(educ_1989_raw)
  )

  if (year <= 2002L) {
    age_years <- convert_age_raw_1990(raw_df$age_raw)
  } else {
    age_years <- convert_age_raw_2003(raw_df$age_raw)
    educ_meta <- harmonize_education_2003(
      educ_1989_raw = educ_1989_raw,
      educ_2003_raw = educ_2003_raw,
      educ_flag_raw = educ_flag_raw
    )
  }

  if (year == 2021L) {
    race_bridged <- harmonize_race_bridged_2021(race_recode_40_raw)
    race_eth <- harmonize_race_eth_2021(
      hispanic_raw = raw_df$hispanic_raw,
      race_recode_40_raw = race_recode_40_raw
    )
  } else if (year >= 2022L) {
    race_bridged <- harmonize_race_bridged_2022plus(race_recode_6_raw)
    race_eth <- harmonize_race_eth_2022plus(
      hispanic_race_recode_2022_raw = hispanic_race_recode_2022_raw,
      race_recode_6_raw = race_recode_6_raw
    )
  }

  if (year >= 1999L) {
    cause_atom <- classify_cause_atom_icd10(raw_df$ucod_raw)
  }

  data.frame(
    year = rep.int(year, nrow(raw_df)),
    file_name = rep.int(file_name, nrow(raw_df)),
    member_name = rep.int(member_name, nrow(raw_df)),
    resident_status_raw = raw_df$resident_status_raw,
    # Resident status 1-3 are U.S. residents; 4 is a foreign resident. This
    # coding is stable across all supported mortality years.
    is_us_resident = raw_df$resident_status_raw %in% c("1", "2", "3"),
    sex_raw = raw_df$sex_raw,
    sex = harmonize_sex(raw_df$sex_raw),
    age_raw = raw_df$age_raw,
    age_recode_12_raw = age_recode_12_raw,
    age_recode_27_raw = age_recode_27_raw,
    age_years = age_years,
    race_raw = race_raw,
    race_bridged = race_bridged,
    race_recode_3_raw = race_recode_3_raw,
    race_recode_2_raw = race_recode_2_raw,
    race_recode_5_raw = race_recode_5_raw,
    race_recode_6_raw = race_recode_6_raw,
    race_recode_40_raw = race_recode_40_raw,
    race_bridged_flag_raw = race_bridged_flag_raw,
    race_imputation_flag_raw = race_imputation_flag_raw,
    hispanic_raw = raw_df$hispanic_raw,
    hispanic_race_recode_raw = hispanic_race_recode_raw,
    hispanic_race_recode_2022_raw = hispanic_race_recode_2022_raw,
    race_eth = race_eth,
    educ_1989_raw = educ_1989_raw,
    educ_1989_recode_raw = educ_1989_recode_raw,
    educ_2003_raw = educ_2003_raw,
    educ_flag_raw = educ_flag_raw,
    educ_raw_selected = educ_meta$educ_raw_selected,
    educ_system = educ_meta$educ_system,
    educ2 = educ_meta$educ2,
    ucod_raw = raw_df$ucod_raw,
    icd_revision = rep.int(ifelse(year <= 1998L, "ICD9", "ICD10"), nrow(raw_df)),
    cause_atom = cause_atom,
    record_weight = rep.int(1L, nrow(raw_df)),
    stringsAsFactors = FALSE
  )
}

get_main_panel_group_vars <- function() {
  c(
    "year",
    "resident_status_raw",
    "sex_raw",
    "sex",
    "age_raw",
    "age_recode_12_raw",
    "age_recode_27_raw",
    "age_years",
    "race_raw",
    "race_recode_3_raw",
    "race_recode_2_raw",
    "race_recode_5_raw",
    "race_recode_6_raw",
    "race_recode_40_raw",
    "race_bridged_flag_raw",
    "race_imputation_flag_raw",
    "race_bridged",
    "hispanic_raw",
    "hispanic_race_recode_raw",
    "hispanic_race_recode_2022_raw",
    "race_eth",
    "educ_1989_raw",
    "educ_1989_recode_raw",
    "educ_2003_raw",
    "educ_flag_raw",
    "educ_raw_selected",
    "educ_system",
    "educ2",
    "ucod_raw",
    "icd_revision",
    "cause_atom"
  )
}

prepare_main_panel_input <- function(chunk_df) {
  # Resident-only table keeping raw + harmonized fields; age/race aggregation
  # outside education follows Schlüter et al., 2024.
  agg_df <- chunk_df[chunk_df$is_us_resident %in% TRUE, c(
    get_main_panel_group_vars(),
    "record_weight"
  )]

  agg_df$sex[agg_df$sex == ""] <- NA_character_
  agg_df$age_years[!is.na(agg_df$age_years) & agg_df$age_years > 85L] <- 85L
  agg_df$race_bridged[agg_df$race_bridged == ""] <- NA_character_
  agg_df$race_eth[agg_df$race_eth == ""] <- NA_character_
  agg_df$educ2[is.na(agg_df$educ2) | agg_df$educ2 == ""] <- "Unknown"

  agg_df
}

# data.table aggregation keeps NA categories as their own group, so missing
# values need no sentinel encoding.
summarize_chunk_main_panel <- function(chunk_df) {
  agg_input <- prepare_main_panel_input(chunk_df)

  if (nrow(agg_input) == 0L) {
    return(NULL)
  }

  group_vars <- get_main_panel_group_vars()
  agg_dt <- as.data.table(agg_input)

  agg_dt[, .(n_deaths = sum(record_weight)), by = group_vars]
}

combine_aggregate_chunks <- function(chunks, group_vars) {
  # Sum partial death counts (across chunks within a year, or across years) over
  # the shared grouping keys; return a stably ordered data.frame.
  chunks <- chunks[!vapply(chunks, is.null, logical(1))]

  if (length(chunks) == 0L) {
    return(NULL)
  }

  combined <- rbindlist(chunks, use.names = TRUE)
  out <- combined[, .(n_deaths = sum(n_deaths)), by = group_vars]
  setorderv(out, group_vars, na.last = TRUE)

  setDF(out)
  rownames(out) <- NULL
  out
}

process_mortality_zip <- function(zip_path, layout_registry, chunk_size, save_minimal_extracts) {
  file_name <- basename(zip_path)
  year <- extract_year_from_file_name(file_name)

  if (is.na(year)) {
    stop("Could not extract a year from file name: ", file_name, call. = FALSE)
  }

  if (!year %in% as.integer(names(layout_registry))) {
    message(
      "Skipping mortality file with unsupported year ",
      year,
      ": ",
      file_name
    )

    return(NULL)
  }

  year_layout <- layout_registry[[as.character(year)]]
  extraction_target <- extract_zip_member_to_temp(zip_path)
  # Read the extracted text member sequentially, so no raw file is held whole
  # in memory.
  file_con <- file(extraction_target$extracted_path, open = "rb")

  on.exit(close(file_con), add = TRUE)
  on.exit(unlink(extraction_target$temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  minimal_chunks <- list()
  main_panel_chunks <- list()

  n_records_total <- 0L
  n_resident_total <- 0L
  n_missing_educ_raw <- 0L
  n_unknown_educ2 <- 0L
  n_missing_race_raw <- 0L
  n_unknown_race_eth <- 0L
  n_missing_ucod_raw <- 0L
  n_missing_age_years <- 0L
  min_record_length <- Inf
  max_record_length <- -Inf

  repeat {
    record_chunk <- readLines(file_con, n = chunk_size, warn = FALSE, encoding = "bytes")

    if (length(record_chunk) == 0L) {
      break
    }

    n_records_total <- n_records_total + length(record_chunk)

    record_lengths <- nchar(record_chunk, type = "bytes", allowNA = FALSE, keepNA = FALSE)
    min_record_length <- min(min_record_length, min(record_lengths))
    max_record_length <- max(max_record_length, max(record_lengths))

    raw_fields <- extract_fields_from_chunk(
      record_chunk = record_chunk,
      fields = year_layout$fields,
      expected_length = year_layout$record_length
    )

    chunk_df <- derive_mortality_fields(
      raw_df = raw_fields,
      year = year,
      file_name = file_name,
      member_name = extraction_target$member_name
    )

    n_resident_total <- n_resident_total + sum(chunk_df$is_us_resident, na.rm = TRUE)
    n_missing_educ_raw <- n_missing_educ_raw + sum(is.na(chunk_df$educ_raw_selected))
    n_unknown_educ2 <- n_unknown_educ2 + sum(chunk_df$educ2 == "Unknown")
    n_missing_race_raw <- n_missing_race_raw + sum(is.na(chunk_df$race_raw))
    n_unknown_race_eth <- n_unknown_race_eth + sum(is.na(chunk_df$race_eth))
    n_missing_ucod_raw <- n_missing_ucod_raw + sum(is.na(chunk_df$ucod_raw))
    n_missing_age_years <- n_missing_age_years + sum(is.na(chunk_df$age_years))

    main_panel_chunks[[length(main_panel_chunks) + 1L]] <- summarize_chunk_main_panel(chunk_df)

    if (isTRUE(save_minimal_extracts)) {
      minimal_chunks[[length(minimal_chunks) + 1L]] <- chunk_df
    }
  }

  main_panel <- combine_aggregate_chunks(
    chunks = main_panel_chunks,
    group_vars = get_main_panel_group_vars()
  )

  minimal_extract_path <- NA_character_

  if (isTRUE(save_minimal_extracts)) {
    # Row-level extract, for re-auditing the layout work.
    minimal_extract <- do.call(rbind, minimal_chunks)
    rownames(minimal_extract) <- NULL

    minimal_extract_path <- file.path(
      extract_dir,
      sprintf("01_death_minimal_extract_%s.rds", year)
    )

    saveRDS(
      object = minimal_extract,
      file = minimal_extract_path,
      compress = "xz"
    )
  }

  processing_metadata <- data.frame(
    data_type = "mortality",
    year = year,
    file_name = file_name,
    member_name = extraction_target$member_name,
    supported_layout = TRUE,
    record_length_expected = year_layout$record_length,
    min_record_length = if (is.finite(min_record_length)) min_record_length else NA_integer_,
    max_record_length = if (is.finite(max_record_length)) max_record_length else NA_integer_,
    chunk_size = chunk_size,
    n_records_total = n_records_total,
    n_records_us_resident = n_resident_total,
    n_missing_educ_raw = n_missing_educ_raw,
    n_unknown_educ2 = n_unknown_educ2,
    n_missing_race_raw = n_missing_race_raw,
    n_unknown_race_eth = n_unknown_race_eth,
    n_missing_ucod_raw = n_missing_ucod_raw,
    n_missing_age_years = n_missing_age_years,
    minimal_extract_saved = isTRUE(save_minimal_extracts),
    minimal_extract_path = minimal_extract_path,
    stringsAsFactors = FALSE
  )

  list(
    main_panel = main_panel,
    metadata = processing_metadata
  )
}

# ---- Layout registry ----

# Fixed-width field positions for the death-numerator fields: residency, sex,
# detailed age, race / Hispanic recodes, education, and UCOD. The layout changes
# across ICD and certificate eras; verify a year block before extending it.
layout_registry <- list(
  `1990` = list(
    record_length = 440L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_1989_recode_raw",
        "sex_raw",
        "race_raw",
        "race_recode_3_raw",
        "race_recode_2_raw",
        "age_raw",
        "age_recode_12_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw",
        "ucod_raw"
      ),
      start = c(20L, 52L, 54L, 59L, 60L, 62L, 63L, 64L, 71L, 80L, 82L, 146L),
      end = c(20L, 53L, 54L, 59L, 61L, 62L, 63L, 66L, 72L, 81L, 82L, 149L),
      stringsAsFactors = FALSE
    )
  ),
  `1999` = list(
    record_length = 440L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_1989_recode_raw",
        "sex_raw",
        "race_raw",
        "race_recode_3_raw",
        "race_recode_2_raw",
        "age_raw",
        "age_recode_12_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw",
        "ucod_raw"
      ),
      start = c(20L, 52L, 54L, 59L, 60L, 62L, 63L, 64L, 71L, 80L, 82L, 142L),
      end = c(20L, 53L, 54L, 59L, 61L, 62L, 63L, 66L, 72L, 81L, 82L, 145L),
      stringsAsFactors = FALSE
    )
  ),
  `2003` = list(
    record_length = 488L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_raw",
        "race_bridged_flag_raw",
        "race_imputation_flag_raw",
        "race_recode_3_raw",
        "race_recode_5_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw"
      ),
      start = c(20L, 61L, 63L, 64L, 69L, 70L, 77L, 146L, 445L, 447L, 448L, 449L, 450L, 484L, 488L),
      end = c(20L, 62L, 63L, 64L, 69L, 73L, 78L, 149L, 446L, 447L, 448L, 449L, 450L, 486L, 488L),
      stringsAsFactors = FALSE
    )
  ),
  `2018` = list(
    record_length = 490L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_raw",
        "race_bridged_flag_raw",
        "race_imputation_flag_raw",
        "race_recode_3_raw",
        "race_recode_5_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw"
      ),
      start = c(20L, 61L, 63L, 64L, 69L, 70L, 77L, 146L, 445L, 447L, 448L, 449L, 450L, 484L, 488L),
      end = c(20L, 62L, 63L, 64L, 69L, 73L, 78L, 149L, 446L, 447L, 448L, 449L, 450L, 486L, 488L),
      stringsAsFactors = FALSE
    )
  ),
  `2019` = list(
    record_length = 490L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_raw",
        "race_bridged_flag_raw",
        "race_imputation_flag_raw",
        "race_recode_3_raw",
        "race_recode_5_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw"
      ),
      start = c(20L, 61L, 63L, 64L, 69L, 70L, 77L, 146L, 445L, 447L, 448L, 449L, 450L, 484L, 488L),
      end = c(20L, 62L, 63L, 64L, 69L, 73L, 78L, 149L, 446L, 447L, 448L, 449L, 450L, 486L, 488L),
      stringsAsFactors = FALSE
    )
  ),
  `2020` = list(
    record_length = 817L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_raw",
        "race_bridged_flag_raw",
        "race_imputation_flag_raw",
        "race_recode_3_raw",
        "race_recode_5_raw",
        "hispanic_raw",
        "hispanic_race_recode_raw"
      ),
      start = c(20L, 61L, 63L, 64L, 69L, 70L, 77L, 146L, 445L, 447L, 448L, 449L, 450L, 484L, 488L),
      end = c(20L, 62L, 63L, 64L, 69L, 73L, 78L, 149L, 446L, 447L, 448L, 449L, 450L, 486L, 488L),
      stringsAsFactors = FALSE
    )
  ),
  `2021` = list(
    record_length = 817L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_1989_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "hispanic_raw",
        "race_recode_40_raw"
      ),
      start = c(20L, 61L, 63L, 64L, 69L, 70L, 77L, 146L, 484L, 489L),
      end = c(20L, 62L, 63L, 64L, 69L, 73L, 78L, 149L, 486L, 490L),
      stringsAsFactors = FALSE
    )
  ),
  `2022` = list(
    record_length = 817L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_imputation_flag_raw",
        "race_recode_6_raw",
        "hispanic_raw",
        "hispanic_race_recode_2022_raw",
        "race_recode_40_raw"
      ),
      start = c(20L, 63L, 64L, 69L, 70L, 77L, 146L, 448L, 450L, 484L, 487L, 489L),
      end = c(20L, 63L, 64L, 69L, 73L, 78L, 149L, 448L, 450L, 486L, 488L, 490L),
      stringsAsFactors = FALSE
    )
  ),
  `2023` = list(
    record_length = 817L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_imputation_flag_raw",
        "race_recode_6_raw",
        "hispanic_raw",
        "hispanic_race_recode_2022_raw",
        "race_recode_40_raw"
      ),
      start = c(20L, 63L, 64L, 69L, 70L, 77L, 146L, 448L, 450L, 484L, 487L, 489L),
      end = c(20L, 63L, 64L, 69L, 73L, 78L, 149L, 448L, 450L, 486L, 488L, 490L),
      stringsAsFactors = FALSE
    )
  ),
  `2024` = list(
    record_length = 817L,
    fields = data.frame(
      field_name = c(
        "resident_status_raw",
        "educ_2003_raw",
        "educ_flag_raw",
        "sex_raw",
        "age_raw",
        "age_recode_27_raw",
        "ucod_raw",
        "race_imputation_flag_raw",
        "race_recode_6_raw",
        "hispanic_raw",
        "hispanic_race_recode_2022_raw",
        "race_recode_40_raw"
      ),
      start = c(20L, 63L, 64L, 69L, 70L, 77L, 146L, 448L, 450L, 484L, 487L, 489L),
      end = c(20L, 63L, 64L, 69L, 73L, 78L, 149L, 448L, 450L, 486L, 488L, 490L),
      stringsAsFactors = FALSE
    )
  )
)

# Expand the anchor-year layouts to the remaining years sharing the same
# fixed-width structure. Era boundaries follow record length and the certificate
# and recode changes.
for (year in 1991:1998) {
  layout_registry[[as.character(year)]] <- layout_registry[["1990"]]
}

for (year in 2000:2002) {
  layout_registry[[as.character(year)]] <- layout_registry[["1999"]]
}

for (year in 2004:2012) {
  layout_registry[[as.character(year)]] <- layout_registry[["2003"]]
}

for (year in 2013:2017) {
  layout_registry[[as.character(year)]] <- layout_registry[["2018"]]
}

# ---- Process the mortality files, year by year ----

for (output_dir in c(proc_dir, extract_dir, per_year_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

mortality_files <- list.files(
  path = raw_dir,
  pattern = "\\.zip$",
  full.names = TRUE
)

if (length(mortality_files) == 0L) {
  stop("No raw NCHS mortality zip files were found under `01_data/01_raw/nchs/mortality/`.", call. = FALSE)
}

file_years <- vapply(mortality_files, extract_year_from_file_name, integer(1))

missing_years <- setdiff(as.integer(names(layout_registry)), unique(file_years))
if (length(missing_years)) {
  stop(
    "The mortality input is incomplete; missing years: ",
    paste(missing_years, collapse = ", "),
    call. = FALSE
  )
}
duplicate_years <- sort(unique(file_years[duplicated(file_years)]))
if (length(duplicate_years)) {
  stop(
    "More than one mortality archive was found for: ",
    paste(duplicate_years, collapse = ", "),
    call. = FALSE
  )
}

mortality_files <- mortality_files[order(file_years, basename(mortality_files))]

# Each year is processed independently and cached under per_year/. Setting
# NCHS_DEATH_FORCE_REPROCESS=false reuses the years already on disk.
main_panel_results <- list()

for (zip_path in mortality_files) {
  year <- extract_year_from_file_name(basename(zip_path))
  cache_path <- file.path(per_year_dir, sprintf("01_death_year_%s.rds", year))

  if (file.exists(cache_path) && !force_reprocess) {
    message("Using cached year ", year, ": ", basename(cache_path))
    processed <- readRDS(cache_path)
  } else {
    message("Processing mortality file: ", basename(zip_path))
    processed <- process_mortality_zip(
      zip_path = zip_path,
      layout_registry = layout_registry,
      chunk_size = chunk_size,
      save_minimal_extracts = save_minimal_extracts
    )

    if (!is.null(processed)) {
      saveRDS(processed, cache_path)
    }
  }

  if (!is.null(processed)) {
    main_panel_results[[length(main_panel_results) + 1L]] <- processed$main_panel
  }
}

if (length(main_panel_results) == 0L) {
  stop("No supported mortality files were processed. Extend `layout_registry` first.", call. = FALSE)
}

main_panel <- combine_aggregate_chunks(
  chunks = main_panel_results,
  group_vars = get_main_panel_group_vars()
)

# ---- Death panel ----

saveRDS(main_panel, main_panel_rds_path)
