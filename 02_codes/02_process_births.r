# ==============================================================================
# 02_process_births.r
#
# Harmonize NCHS natality records across certificate revisions and build the
# resident-birth count panel.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(here)
  library(data.table)
})

here::i_am("02_codes/02_process_births.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Natality processing settings ----

chunk_size <- {
  raw_value <- Sys.getenv("NCHS_BIRTH_CHUNK_SIZE", unset = "100000")
  parsed <- suppressWarnings(as.integer(raw_value))
  if (is.na(parsed) || parsed <= 0L) 100000L else parsed
}

# Formal runs rebuild every year. Set false only to resume an interrupted run.
force_reprocess <- {
  raw_value <- tolower(trimws(Sys.getenv("NCHS_BIRTH_FORCE_REPROCESS", unset = "true")))
  raw_value %in% c("true", "t", "1", "yes", "y")
}

# ---- Natality files ----

raw_dir <- here("01_data", "01_raw", "nchs", "natality")
proc_dir <- here("01_data", "02_processed", "births")
per_year_dir <- file.path(proc_dir, "per_year")

main_panel_rds_path <- file.path(proc_dir, "01_birth_main_panel.rds")

# ---- Natality readers and recodes ----

extract_year_from_file_name <- function(file_name) {
  year_text <- sub(".*?(\\d{4}).*", "\\1", basename(file_name))
  if (!grepl("^\\d{4}$", year_text)) return(NA_integer_)
  as.integer(year_text)
}

escape_powershell_single_quoted <- function(x) {
  gsub("'", "''", x, fixed = TRUE)
}

get_zip_member_name <- function(zip_path) {
  zip_listing <- utils::unzip(zip_path, list = TRUE)
  member_names <- zip_listing$Name[!grepl("/$", zip_listing$Name)]
  if (!length(member_names)) stop("No file member found inside: ", zip_path, call. = FALSE)
  member_names[1]
}

extract_zip_member_with_powershell <- function(zip_path, member_name, temp_dir) {
  ps_exe <- Sys.which("pwsh")
  if (identical(ps_exe, "")) ps_exe <- Sys.which("powershell")
  if (identical(ps_exe, "")) stop("PowerShell not found.", call. = FALSE)

  zip_path <- normalizePath(zip_path, winslash = "\\", mustWork = TRUE)
  temp_dir <- normalizePath(temp_dir, winslash = "\\", mustWork = TRUE)
  extracted_path <- file.path(temp_dir, member_name)
  dir.create(dirname(extracted_path), recursive = TRUE, showWarnings = FALSE)

  ps_command <- paste0(
    "Add-Type -AssemblyName System.IO.Compression.FileSystem; ",
    "$zip = [System.IO.Compression.ZipFile]::OpenRead('", escape_powershell_single_quoted(zip_path), "'); ",
    "$entry = $zip.GetEntry('", escape_powershell_single_quoted(member_name), "'); ",
    "if ($null -eq $entry) { $zip.Dispose(); throw 'Zip member not found.' }; ",
    "$dest = '", escape_powershell_single_quoted(extracted_path), "'; ",
    "$destDir = [System.IO.Path]::GetDirectoryName($dest); ",
    "if (-not [string]::IsNullOrEmpty($destDir)) { [System.IO.Directory]::CreateDirectory($destDir) | Out-Null }; ",
    "[System.IO.Compression.ZipFileExtensions]::ExtractToFile($entry, $dest, $true); ",
    "$zip.Dispose()"
  )

  ps_output <- suppressWarnings(system2(ps_exe, c("-NoProfile", "-Command", ps_command), stdout = TRUE, stderr = TRUE))
  ps_status <- attr(ps_output, "status")

  if (!is.null(ps_status) && ps_status != 0L) {
    stop("PowerShell unzip fallback failed for ", basename(zip_path), ": ", paste(ps_output, collapse = "\n"), call. = FALSE)
  }

  if (!file.exists(extracted_path)) {
    stop("Expected extracted file was not created: ", extracted_path, call. = FALSE)
  }
  if (is.na(file.info(extracted_path)$size) || file.info(extracted_path)$size <= 0) {
    stop("Extracted file is empty after PowerShell unzip fallback: ", extracted_path, call. = FALSE)
  }

  extracted_path
}

extract_zip_member_with_7zip <- function(zip_path, member_name, temp_dir) {
  seven_zip <- Sys.which("7z")
  if (identical(seven_zip, "")) stop("7z not found.", call. = FALSE)
  ps_exe <- Sys.which("pwsh")
  if (identical(ps_exe, "")) ps_exe <- Sys.which("powershell")
  if (identical(ps_exe, "")) stop("PowerShell not found for 7z fallback.", call. = FALSE)

  zip_path <- normalizePath(zip_path, winslash = "\\", mustWork = TRUE)
  temp_dir <- normalizePath(temp_dir, winslash = "\\", mustWork = TRUE)
  seven_zip <- normalizePath(seven_zip, winslash = "\\", mustWork = TRUE)

  ps_command <- paste0(
    "& '", escape_powershell_single_quoted(seven_zip), "' ",
    "e -y -bd ('-o' + '", escape_powershell_single_quoted(temp_dir), "') ",
    "'", escape_powershell_single_quoted(zip_path), "' ",
    "'", escape_powershell_single_quoted(member_name), "'"
  )

  seven_zip_output <- suppressWarnings(system2(
    ps_exe,
    c("-NoProfile", "-Command", ps_command),
    stdout = TRUE,
    stderr = TRUE
  ))
  seven_zip_status <- attr(seven_zip_output, "status")

  if (!is.null(seven_zip_status) && seven_zip_status != 0L) {
    stop(
      "7z unzip fallback failed for ",
      basename(zip_path),
      ": ",
      paste(seven_zip_output, collapse = "\n"),
      call. = FALSE
    )
  }

  extracted_path <- file.path(temp_dir, basename(member_name))
  if (!file.exists(extracted_path)) {
    stop("Expected 7z extracted file was not created: ", extracted_path, call. = FALSE)
  }
  if (is.na(file.info(extracted_path)$size) || file.info(extracted_path)$size <= 0) {
    stop("Extracted file is empty after 7z unzip fallback: ", extracted_path, call. = FALSE)
  }

  extracted_path
}

extract_zip_member_to_temp <- function(zip_path) {
  member_name <- get_zip_member_name(zip_path)
  temp_dir <- tempfile(pattern = "nchs_birth_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  extracted_path <- tryCatch(
    suppressWarnings(utils::unzip(zipfile = zip_path, files = member_name, exdir = temp_dir, overwrite = TRUE)),
    error = function(e) character(0)
  )
  extracted_path <- extracted_path[file.exists(extracted_path)]
  if (length(extracted_path)) {
    extracted_path <- extracted_path[!is.na(file.info(extracted_path)$size) & file.info(extracted_path)$size > 0]
  }
  if (!length(extracted_path)) {
    extracted_path <- tryCatch(
      extract_zip_member_with_powershell(zip_path, member_name, temp_dir),
      error = function(e) NULL
    )
  } else {
    extracted_path <- extracted_path[[1]]
  }

  if (
    is.null(extracted_path) ||
    !length(extracted_path) ||
    !file.exists(extracted_path) ||
    is.na(file.info(extracted_path)$size) ||
    file.info(extracted_path)$size <= 0
  ) {
    extracted_path <- extract_zip_member_with_7zip(zip_path, member_name, temp_dir)
  } else {
    extracted_path <- extracted_path[[1]]
  }

  list(member_name = member_name, extracted_path = extracted_path, temp_dir = temp_dir)
}

trim_fixed_width_field <- function(x) {
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

pad_record_chunk <- function(record_chunk, min_record_length) {
  record_lengths <- nchar(record_chunk, type = "bytes", allowNA = FALSE, keepNA = FALSE)
  short_records <- record_lengths < min_record_length
  if (any(short_records)) {
    record_chunk[short_records] <- paste0(record_chunk[short_records], strrep(" ", min_record_length - record_lengths[short_records]))
  }
  record_chunk
}

extract_fields_from_chunk <- function(record_chunk, fields, min_record_length) {
  record_chunk <- pad_record_chunk(record_chunk, min_record_length)
  extracted <- lapply(
    seq_len(nrow(fields)),
    function(i) trim_fixed_width_field(substr(record_chunk, fields$start[i], fields$end[i]))
  )
  as.data.frame(setNames(extracted, fields$field_name), stringsAsFactors = FALSE)
}

# -- Field harmonization: education, age, race/ethnicity --

select_parent_education <- function(educ_1989_raw, educ_2003_raw, educ_flag_raw) {
  use_revised <- !is.na(educ_2003_raw) & educ_2003_raw != ""
  use_1989 <- !use_revised & !is.na(educ_1989_raw) & educ_1989_raw != ""
  educ2_1989 <- recode_educ2_from_1989(educ_1989_raw)
  educ2_2003 <- recode_educ2_from_2003(educ_2003_raw)
  out_version <- rep(NA_character_, length(use_revised))
  out_raw <- rep(NA_character_, length(use_revised))
  out_educ2 <- rep("Unknown", length(use_revised))

  # Reporting flags do not fully substitute for certificate version, so the
  # version label follows the source field actually used.
  if (!is.null(educ_flag_raw)) {
    out_version[!is.na(educ_flag_raw) & educ_flag_raw == "1" & !use_1989] <- "2003_revision"
  }

  out_version[use_revised] <- "2003_revision"
  out_version[is.na(out_version) & use_1989] <- "1989_revision"
  out_raw[use_revised] <- educ_2003_raw[use_revised]
  out_raw[use_1989] <- educ_1989_raw[use_1989]
  out_educ2[use_revised] <- educ2_2003[use_revised]
  out_educ2[use_1989] <- educ2_1989[use_1989]

  list(educ_version_used = out_version, educ_raw_selected = out_raw, educ2 = out_educ2)
}

resident_birth_flag <- function(resident_status_raw) {
  status_num <- safe_integer(resident_status_raw)
  !is.na(status_num) & status_num %in% 1:3
}

# Mother age analytic window is 10-54. The effective bounds differ by era with
# the source coding: 1990-2002 carries literal single years 10-54; 2003 is the
# MAGER41 recode (see harmonize_mother_age_transition); 2004+ single-year MAGER
# is bottom-coded at 12 and top-coded at 50 (50 = "50 years and over").
harmonize_mother_age_pre2003 <- function(mother_age_raw) {
  raw_num <- safe_integer(mother_age_raw)
  out <- rep(NA_integer_, length(raw_num))
  out[!is.na(raw_num) & raw_num >= 10L & raw_num <= 54L] <- raw_num[!is.na(raw_num) & raw_num >= 10L & raw_num <= 54L]
  out
}

harmonize_mother_age_2003plus <- function(mother_age_raw) {
  raw_num <- safe_integer(mother_age_raw)
  out <- rep(NA_integer_, length(raw_num))
  out[!is.na(raw_num) & raw_num >= 12L & raw_num <= 49L] <- raw_num[!is.na(raw_num) & raw_num >= 12L & raw_num <= 49L]
  out[!is.na(raw_num) & raw_num == 50L] <- 50L
  out
}

harmonize_mother_age_transition <- function(mother_age_raw) {
  raw_num <- safe_integer(mother_age_raw)
  out <- rep(NA_integer_, length(raw_num))

  # The 2003 public natality file stores mother age in MAGER41: 01 = under 15,
  # 02 = 15, ..., 41 = 54. Mapped back to single years, with the under-15
  # bucket taken as age 14 and codes 02:41 as 15:54.
  valid_idx <- !is.na(raw_num) & raw_num >= 1L & raw_num <= 41L
  out[valid_idx] <- raw_num[valid_idx] + 13L
  out
}

harmonize_father_age <- function(father_age_raw) {
  raw_num <- safe_integer(father_age_raw)
  out <- rep(NA_integer_, length(raw_num))

  # In the natality public-use files used here, father age is stored as a
  # direct single-year age in the raw age field. Blank / 99 represent unknown.
  valid_idx <- !is.na(raw_num) & raw_num >= 10L & raw_num <= 98L
  out[valid_idx] <- raw_num[valid_idx]
  out
}

# Bridged-era Hispanic origin recode: 1-5 = Hispanic subtypes, 9 = origin
# not stated. Origin-unknown (9) is folded into Hispanic to match the death
# numerator (Schlüter et al., 2024: hspanicr %in% c(1:5, 9) ~ "hispanic").
is_hispanic_old <- function(hispanic_raw) {
  raw_num <- safe_integer(hispanic_raw)
  !is.na(raw_num) & raw_num %in% c(1L, 2L, 3L, 4L, 5L, 9L)
}

# Bridged-era (2003-2020) Hispanic origin recode, matching the death numerator
# (Schlüter et al., 2024: 1-5 plus not-stated 9 -> hispanic). The 2021+
# single-race era treats not-stated origin differently (see
# harmonize_race_eth_2021plus), so it does not use this helper.
is_hispanic_2003plus <- function(hispanic_raw) {
  raw_num <- safe_integer(hispanic_raw)
  !is.na(raw_num) & raw_num %in% c(1L, 2L, 3L, 4L, 5L, 9L)
}

# Origin not stated, coded 9 in every era. The death numerator folds it into
# Hispanic; the father path marks it not stated instead, because unstated origin
# runs 15-16% on the father side against 3-5% on the mother side.
is_origin_not_stated <- function(hispanic_raw) {
  raw_num <- safe_integer(hispanic_raw)
  !is.na(raw_num) & raw_num == 9L
}

map_old_race_to_group <- function(race_raw) {
  raw_num <- safe_integer(race_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num == 1L] <- "white"
  out[!is.na(raw_num) & raw_num == 2L] <- "black"
  out[!is.na(raw_num) & raw_num == 3L] <- "aian"
  # API codes match the death numerator: single Asian/PI races 4-7 plus
  # combined-race codes 18..78.
  out[!is.na(raw_num) & raw_num %in% c(4L, 5L, 6L, 7L, 18L, 28L, 38L, 48L, 58L, 68L, 78L)] <- "api"
  out
}

map_race31_to_group <- function(race31_raw) {
  raw_num <- safe_integer(race31_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num == 1L] <- "white"
  out[!is.na(raw_num) & raw_num == 2L] <- "black"
  out[!is.na(raw_num) & raw_num == 3L] <- "aian"
  out[!is.na(raw_num) & raw_num %in% c(4L, 5L)] <- "api"
  out
}

map_race6_to_group <- function(race6_raw) {
  raw_num <- safe_integer(race6_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num == 1L] <- "white"
  out[!is.na(raw_num) & raw_num == 2L] <- "black"
  out[!is.na(raw_num) & raw_num == 3L] <- "aian"
  out[!is.na(raw_num) & raw_num %in% c(4L, 5L)] <- "api"
  out
}

map_public_race_eth_recode <- function(race_eth_public_raw) {
  raw_num <- safe_integer(race_eth_public_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num == 1L] <- "white"
  out[!is.na(raw_num) & raw_num == 2L] <- "black"
  out[!is.na(raw_num) & raw_num == 3L] <- "aian"
  out[!is.na(raw_num) & raw_num %in% c(4L, 5L)] <- "api"
  out[!is.na(raw_num) & raw_num == 7L] <- "hispanic"
  out
}

map_pre2003_public_race_eth_recode <- function(race_eth_public_raw) {
  raw_num <- safe_integer(race_eth_public_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num == 6L] <- "white"
  out[!is.na(raw_num) & raw_num == 7L] <- "black"
  out[!is.na(raw_num) & raw_num %in% c(1L, 2L, 3L, 4L, 5L, 11L, 22L, 33L, 44L, 55L)] <- "hispanic"
  out
}

map_transition_public_race_eth_recode <- function(race_eth_public_raw) {
  raw_num <- safe_integer(race_eth_public_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num %in% c(1L, 2L, 3L, 4L, 5L)] <- "hispanic"
  out[!is.na(raw_num) & raw_num == 6L] <- "white"
  out[!is.na(raw_num) & raw_num == 7L] <- "black"
  out
}

map_transition_race_to_group <- function(race_raw) {
  raw_num <- safe_integer(race_raw)
  out <- rep(NA_character_, length(raw_num))
  out[!is.na(raw_num) & raw_num %in% c(1L, 11L)] <- "white"
  out[!is.na(raw_num) & raw_num %in% c(2L, 22L)] <- "black"
  out[!is.na(raw_num) & raw_num %in% c(3L, 33L)] <- "aian"
  out[!is.na(raw_num) & raw_num %in% c(4L, 44L, 54L, 64L, 74L, 84L)] <- "api"
  out
}

harmonize_race_eth_pre2003 <- function(race_eth_public_raw, hispanic_raw, race_raw,
                                       origin_not_stated_to_na = FALSE) {
  out <- map_pre2003_public_race_eth_recode(race_eth_public_raw)
  unresolved <- is.na(out)
  if (any(unresolved)) {
    fallback <- map_old_race_to_group(race_raw)
    fallback[is_hispanic_old(hispanic_raw)] <- "hispanic"
    out[unresolved] <- fallback[unresolved]
  }
  if (origin_not_stated_to_na) out[is_origin_not_stated(hispanic_raw)] <- NA_character_
  out
}

harmonize_race_eth_transition <- function(race_eth_public_raw, hispanic_raw, race_raw,
                                          origin_not_stated_to_na = FALSE) {
  out <- map_transition_public_race_eth_recode(race_eth_public_raw)
  unresolved <- is.na(out)
  if (any(unresolved)) {
    fallback <- map_transition_race_to_group(race_raw)
    fallback[is_hispanic_2003plus(hispanic_raw)] <- "hispanic"
    out[unresolved] <- fallback[unresolved]
  }
  out[is_hispanic_2003plus(hispanic_raw)] <- "hispanic"
  if (origin_not_stated_to_na) out[is_origin_not_stated(hispanic_raw)] <- NA_character_
  out
}

harmonize_race_eth_2016plus <- function(hispanic_raw, race31_raw, race6_raw = NULL, race_eth_public_raw = NULL,
                                        origin_not_stated_to_na = FALSE) {
  out <- if (is.null(race_eth_public_raw)) {
    map_race31_to_group(race31_raw)
  } else {
    mapped <- map_public_race_eth_recode(race_eth_public_raw)
    unresolved <- is.na(mapped)
    if (any(unresolved)) {
      race31_group <- map_race31_to_group(race31_raw)
      mapped[unresolved] <- race31_group[unresolved]
      still_unresolved <- is.na(mapped) & !is.null(race6_raw)
      if (any(still_unresolved)) {
        mapped[still_unresolved] <- map_race6_to_group(race6_raw)[still_unresolved]
      }
    }
    mapped
  }
  out[is_hispanic_2003plus(hispanic_raw)] <- "hispanic"
  if (origin_not_stated_to_na) out[is_origin_not_stated(hispanic_raw)] <- NA_character_
  out
}

# 2021-2024 single-race era: race comes from the single-race recode (recode 15,
# multi-race -> NA), falling back to recode 6, with Hispanic origin taking
# precedence. Origin not stated (recode 9) -> NA, as in the death numerator
# (Villaveces et al., 2025).
harmonize_race_eth_2021plus <- function(hispanic_recode_raw, race15_raw, race6_raw = NULL, race_eth_public_raw = NULL,
                                        origin_not_stated_to_na = FALSE) {
  out <- map_race31_to_group(race15_raw)
  if (!is.null(race6_raw)) {
    unresolved <- is.na(out)
    if (any(unresolved)) out[unresolved] <- map_race6_to_group(race6_raw)[unresolved]
  }
  hispanic_num <- safe_integer(hispanic_recode_raw)
  out[!is.na(hispanic_num) & hispanic_num %in% c(1L, 2L, 3L, 4L, 5L)] <- "hispanic"
  # A no-op for this era, which already marks origin not stated for both parents.
  out[!is.na(hispanic_num) & hispanic_num == 9L] <- NA_character_
  out
}

birth_main_panel_columns <- c(
  "year", "resident_status_raw", "resident_birth_flag",
  "mother_age_raw", "mother_age_years",
  "father_age_raw", "father_age_years",
  "mother_race_raw", "mother_race_recode6_raw", "mother_race_recode15_raw",
  "mother_hispanic_raw", "mother_hispanic_recode_raw", "mother_race_eth_public_raw", "mother_race_eth",
  "mother_educ_flag_raw", "mother_educ_1989_raw", "mother_educ_1989_recode_raw",
  "mother_educ_2003_raw", "mother_educ_version_used", "mother_educ_raw_selected", "mother_educ2",
  "father_race_raw", "father_race_recode6_raw", "father_race_recode15_raw",
  "father_hispanic_raw", "father_hispanic_recode_raw", "father_race_eth_public_raw", "father_race_eth",
  "father_educ_flag_raw", "father_educ_1989_raw", "father_educ_2003_raw",
  "father_educ_version_used", "father_educ_raw_selected", "father_educ2",
  "record_weight", "n_births"
)

empty_birth_panel <- function() {
  out <- as.data.frame(setNames(vector("list", length(birth_main_panel_columns)), birth_main_panel_columns), stringsAsFactors = FALSE)
  out$year <- integer(0)
  out$resident_birth_flag <- logical(0)
  out$mother_age_years <- integer(0)
  out$father_age_years <- integer(0)
  out$record_weight <- integer(0)
  out$n_births <- integer(0)
  out
}

align_birth_panel_columns <- function(df) {
  missing_cols <- setdiff(birth_main_panel_columns, names(df))
  if (length(missing_cols)) {
    for (col_name in missing_cols) df[[col_name]] <- rep(NA, nrow(df))
  }
  df <- df[, birth_main_panel_columns, drop = FALSE]
  df$year <- safe_integer(df$year)
  df$resident_birth_flag <- as.logical(df$resident_birth_flag)
  df$mother_age_years <- safe_integer(df$mother_age_years)
  df$father_age_years <- safe_integer(df$father_age_years)
  df$record_weight <- safe_integer(df$record_weight)
  df$n_births <- safe_integer(df$n_births)
  df
}

# -- Birth chunk harmonizers and per-year panel assembly --

harmonize_birth_chunk_pre2003 <- function(raw_chunk, year) {
  mother_education <- select_parent_education(raw_chunk$mother_educ_1989_raw, rep(NA_character_, nrow(raw_chunk)), rep(NA_character_, nrow(raw_chunk)))
  father_education <- select_parent_education(raw_chunk$father_educ_1989_raw, rep(NA_character_, nrow(raw_chunk)), rep(NA_character_, nrow(raw_chunk)))

  out <- data.frame(
    year = year,
    resident_status_raw = raw_chunk$resident_status_raw,
    resident_birth_flag = resident_birth_flag(raw_chunk$resident_status_raw),
    mother_age_raw = raw_chunk$mother_age_raw,
    mother_age_years = harmonize_mother_age_pre2003(raw_chunk$mother_age_raw),
    father_age_raw = raw_chunk$father_age_raw,
    father_age_years = harmonize_father_age(raw_chunk$father_age_raw),
    mother_race_raw = raw_chunk$mother_race_raw,
    mother_race_recode6_raw = NA_character_,
    mother_race_recode15_raw = NA_character_,
    mother_hispanic_raw = raw_chunk$mother_hispanic_raw,
    mother_hispanic_recode_raw = NA_character_,
    mother_race_eth_public_raw = raw_chunk$mother_race_eth_public_raw,
    mother_race_eth = harmonize_race_eth_pre2003(raw_chunk$mother_race_eth_public_raw, raw_chunk$mother_hispanic_raw, raw_chunk$mother_race_raw),
    mother_educ_flag_raw = NA_character_,
    mother_educ_1989_raw = raw_chunk$mother_educ_1989_raw,
    mother_educ_1989_recode_raw = raw_chunk$mother_educ_1989_recode_raw,
    mother_educ_2003_raw = NA_character_,
    mother_educ_version_used = mother_education$educ_version_used,
    mother_educ_raw_selected = mother_education$educ_raw_selected,
    mother_educ2 = mother_education$educ2,
    father_race_raw = raw_chunk$father_race_raw,
    father_race_recode6_raw = NA_character_,
    father_race_recode15_raw = NA_character_,
    father_hispanic_raw = raw_chunk$father_hispanic_raw,
    father_hispanic_recode_raw = NA_character_,
    father_race_eth_public_raw = raw_chunk$father_race_eth_public_raw,
    father_race_eth = harmonize_race_eth_pre2003(raw_chunk$father_race_eth_public_raw, raw_chunk$father_hispanic_raw, raw_chunk$father_race_raw, origin_not_stated_to_na = TRUE),
    father_educ_flag_raw = NA_character_,
    father_educ_1989_raw = raw_chunk$father_educ_1989_raw,
    father_educ_2003_raw = NA_character_,
    father_educ_version_used = father_education$educ_version_used,
    father_educ_raw_selected = father_education$educ_raw_selected,
    father_educ2 = father_education$educ2,
    record_weight = 1L,
    n_births = 1L,
    stringsAsFactors = FALSE
  )

  align_birth_panel_columns(out)
}

harmonize_birth_chunk_transition <- function(raw_chunk, year, age_harmonizer = harmonize_mother_age_transition) {
  mother_education <- select_parent_education(raw_chunk$mother_educ_1989_raw, raw_chunk$mother_educ_2003_raw, raw_chunk$mother_educ_flag_raw)
  father_education <- select_parent_education(rep(NA_character_, nrow(raw_chunk)), raw_chunk$father_educ_2003_raw, raw_chunk$father_educ_flag_raw)
  father_race_selected <- ifelse(
    !is.na(raw_chunk$father_race_raw_revised) & raw_chunk$father_race_raw_revised != "",
    raw_chunk$father_race_raw_revised,
    raw_chunk$father_race_raw_old
  )
  # The combined father age field is blank for every record of the 2012 and 2013
  # files, which carry the father's age only in the revised-certificate field.
  # It is never blank in 2003-2011, so the fallback fires only in those two years.
  father_age_selected <- ifelse(
    !is.na(raw_chunk$father_age_raw) & raw_chunk$father_age_raw != "",
    raw_chunk$father_age_raw,
    raw_chunk$father_age_alt_raw
  )

  out <- data.frame(
    year = year,
    resident_status_raw = raw_chunk$resident_status_raw,
    resident_birth_flag = resident_birth_flag(raw_chunk$resident_status_raw),
    mother_age_raw = raw_chunk$mother_age_raw,
    mother_age_years = age_harmonizer(raw_chunk$mother_age_raw),
    father_age_raw = father_age_selected,
    father_age_years = harmonize_father_age(father_age_selected),
    mother_race_raw = raw_chunk$mother_race_raw,
    mother_race_recode6_raw = NA_character_,
    mother_race_recode15_raw = NA_character_,
    mother_hispanic_raw = raw_chunk$mother_hispanic_raw,
    mother_hispanic_recode_raw = NA_character_,
    mother_race_eth_public_raw = raw_chunk$mother_race_eth_public_raw,
    mother_race_eth = harmonize_race_eth_transition(raw_chunk$mother_race_eth_public_raw, raw_chunk$mother_hispanic_raw, raw_chunk$mother_race_raw),
    mother_educ_flag_raw = raw_chunk$mother_educ_flag_raw,
    mother_educ_1989_raw = raw_chunk$mother_educ_1989_raw,
    mother_educ_1989_recode_raw = raw_chunk$mother_educ_1989_recode_raw,
    mother_educ_2003_raw = raw_chunk$mother_educ_2003_raw,
    mother_educ_version_used = mother_education$educ_version_used,
    mother_educ_raw_selected = mother_education$educ_raw_selected,
    mother_educ2 = mother_education$educ2,
    father_race_raw = father_race_selected,
    father_race_recode6_raw = NA_character_,
    father_race_recode15_raw = NA_character_,
    father_hispanic_raw = raw_chunk$father_hispanic_raw,
    father_hispanic_recode_raw = NA_character_,
    father_race_eth_public_raw = raw_chunk$father_race_eth_public_raw,
    father_race_eth = harmonize_race_eth_transition(raw_chunk$father_race_eth_public_raw, raw_chunk$father_hispanic_raw, father_race_selected, origin_not_stated_to_na = TRUE),
    father_educ_flag_raw = raw_chunk$father_educ_flag_raw,
    father_educ_1989_raw = NA_character_,
    father_educ_2003_raw = raw_chunk$father_educ_2003_raw,
    father_educ_version_used = father_education$educ_version_used,
    father_educ_raw_selected = father_education$educ_raw_selected,
    father_educ2 = father_education$educ2,
    record_weight = 1L,
    n_births = 1L,
    stringsAsFactors = FALSE
  )

  align_birth_panel_columns(out)
}

harmonize_birth_chunk_transition_direct_age <- function(raw_chunk, year) {
  harmonize_birth_chunk_transition(raw_chunk, year, age_harmonizer = harmonize_mother_age_2003plus)
}

harmonize_birth_chunk_2016plus <- function(raw_chunk, year, race_eth_harmonizer = harmonize_race_eth_2016plus) {
  mother_education <- select_parent_education(rep(NA_character_, nrow(raw_chunk)), raw_chunk$mother_educ_2003_raw, raw_chunk$mother_educ_flag_raw)
  father_education <- select_parent_education(rep(NA_character_, nrow(raw_chunk)), raw_chunk$father_educ_2003_raw, raw_chunk$father_educ_flag_raw)

  out <- data.frame(
    year = year,
    resident_status_raw = raw_chunk$resident_status_raw,
    resident_birth_flag = resident_birth_flag(raw_chunk$resident_status_raw),
    mother_age_raw = raw_chunk$mother_age_raw,
    mother_age_years = harmonize_mother_age_2003plus(raw_chunk$mother_age_raw),
    father_age_raw = raw_chunk$father_age_raw,
    father_age_years = harmonize_father_age(raw_chunk$father_age_raw),
    mother_race_raw = raw_chunk$mother_race_raw,
    mother_race_recode6_raw = raw_chunk$mother_race_recode6_raw,
    mother_race_recode15_raw = raw_chunk$mother_race_recode15_raw,
    mother_hispanic_raw = raw_chunk$mother_hispanic_raw,
    mother_hispanic_recode_raw = raw_chunk$mother_hispanic_recode_raw,
    mother_race_eth_public_raw = raw_chunk$mother_race_eth_public_raw,
    mother_race_eth = race_eth_harmonizer(
      raw_chunk$mother_hispanic_recode_raw,
      raw_chunk$mother_race_recode15_raw,
      raw_chunk$mother_race_recode6_raw,
      raw_chunk$mother_race_eth_public_raw
    ),
    mother_educ_flag_raw = raw_chunk$mother_educ_flag_raw,
    mother_educ_1989_raw = NA_character_,
    mother_educ_1989_recode_raw = NA_character_,
    mother_educ_2003_raw = raw_chunk$mother_educ_2003_raw,
    mother_educ_version_used = mother_education$educ_version_used,
    mother_educ_raw_selected = mother_education$educ_raw_selected,
    mother_educ2 = mother_education$educ2,
    father_race_raw = raw_chunk$father_race_raw,
    father_race_recode6_raw = raw_chunk$father_race_recode6_raw,
    father_race_recode15_raw = raw_chunk$father_race_recode15_raw,
    father_hispanic_raw = raw_chunk$father_hispanic_raw,
    father_hispanic_recode_raw = raw_chunk$father_hispanic_recode_raw,
    father_race_eth_public_raw = raw_chunk$father_race_eth_public_raw,
    father_race_eth = race_eth_harmonizer(
      raw_chunk$father_hispanic_recode_raw,
      raw_chunk$father_race_recode15_raw,
      raw_chunk$father_race_recode6_raw,
      raw_chunk$father_race_eth_public_raw,
      origin_not_stated_to_na = TRUE
    ),
    father_educ_flag_raw = raw_chunk$father_educ_flag_raw,
    father_educ_1989_raw = NA_character_,
    father_educ_2003_raw = raw_chunk$father_educ_2003_raw,
    father_educ_version_used = father_education$educ_version_used,
    father_educ_raw_selected = father_education$educ_raw_selected,
    father_educ2 = father_education$educ2,
    record_weight = 1L,
    n_births = 1L,
    stringsAsFactors = FALSE
  )

  align_birth_panel_columns(out)
}

harmonize_birth_chunk_2021plus <- function(raw_chunk, year) {
  harmonize_birth_chunk_2016plus(raw_chunk, year, race_eth_harmonizer = harmonize_race_eth_2021plus)
}

get_birth_group_vars <- function() {
  setdiff(birth_main_panel_columns, c("record_weight", "n_births"))
}

# data.table aggregation keeps NA categories as their own group, so missing
# values need no sentinel encoding.
aggregate_resident_births <- function(resident_df) {
  if (nrow(resident_df) == 0L) {
    return(NULL)
  }
  group_vars <- get_birth_group_vars()
  agg_dt <- as.data.table(resident_df)
  agg_dt[, .(n_births = sum(n_births)), by = group_vars]
}

combine_aggregate_chunks <- function(chunks, group_vars) {
  # Sum partial birth counts (across chunks within a year, or across years) over
  # the shared grouping keys; return a stably ordered, type-aligned panel.
  chunks <- chunks[!vapply(chunks, is.null, logical(1))]
  if (length(chunks) == 0L) {
    return(empty_birth_panel())
  }

  combined <- rbindlist(chunks, use.names = TRUE, fill = TRUE)
  out <- combined[, .(n_births = sum(n_births)), by = group_vars]
  setorderv(out, group_vars, na.last = TRUE)

  setDF(out)
  rownames(out) <- NULL
  out$record_weight <- 1L
  align_birth_panel_columns(out)
}

summarize_birth_processing <- function(year, zip_path, layout_name, min_record_length,
                                       observed_record_length_min, observed_record_length_max,
                                       n_records_total, processing_counts) {
  data.frame(
    year = year,
    source_file = basename(zip_path),
    layout_name = layout_name,
    min_required_record_length = min_record_length,
    observed_record_length_min = observed_record_length_min,
    observed_record_length_max = observed_record_length_max,
    n_records_total = n_records_total,
    n_resident_births = processing_counts$n_resident,
    n_missing_mother_age_years = processing_counts$n_missing_mother_age,
    n_missing_father_age_years = processing_counts$n_missing_father_age,
    n_unknown_mother_educ2 = processing_counts$n_unknown_mother_educ2,
    n_unknown_father_educ2 = processing_counts$n_unknown_father_educ2,
    n_missing_mother_race_eth = processing_counts$n_missing_mother_race_eth,
    n_missing_father_race_eth = processing_counts$n_missing_father_race_eth,
    n_nonmissing_father_race_raw = processing_counts$n_nonmissing_father_race_raw,
    n_nonmissing_father_educ_2003_raw = processing_counts$n_nonmissing_father_educ_2003_raw,
    stringsAsFactors = FALSE
  )
}

process_birth_year_file <- function(zip_path, year, layout, chunk_size) {
  extraction <- extract_zip_member_to_temp(zip_path)
  on.exit(unlink(extraction$temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  file_con <- file(extraction$extracted_path, open = "rb")
  on.exit(close(file_con), add = TRUE)

  aggregated_chunks <- list()
  n_records_total <- 0L
  observed_record_length_min <- Inf
  observed_record_length_max <- -Inf
  processing_counts <- list(
    n_resident = 0L,
    n_missing_mother_age = 0L,
    n_missing_father_age = 0L,
    n_unknown_mother_educ2 = 0L,
    n_unknown_father_educ2 = 0L,
    n_missing_mother_race_eth = 0L,
    n_missing_father_race_eth = 0L,
    n_nonmissing_father_race_raw = 0L,
    n_nonmissing_father_educ_2003_raw = 0L
  )

  repeat {
    record_chunk <- readLines(file_con, n = chunk_size, warn = FALSE, encoding = "bytes")
    if (!length(record_chunk)) break
    chunk_lengths <- nchar(record_chunk, type = "bytes")
    n_records_total <- n_records_total + length(chunk_lengths)
    observed_record_length_min <- min(observed_record_length_min, min(chunk_lengths))
    observed_record_length_max <- max(observed_record_length_max, max(chunk_lengths))

    raw_fields <- extract_fields_from_chunk(record_chunk, layout$fields, layout$min_record_length)
    harmonized_chunk <- layout$harmonizer(raw_fields, year)
    resident_chunk <- align_birth_panel_columns(
      harmonized_chunk[harmonized_chunk$resident_birth_flag %in% TRUE, , drop = FALSE]
    )

    if (nrow(resident_chunk)) {
      # Missingness counts accumulate per record, before aggregation.
      processing_counts$n_resident <- processing_counts$n_resident + nrow(resident_chunk)
      processing_counts$n_missing_mother_age <- processing_counts$n_missing_mother_age + sum(is.na(resident_chunk$mother_age_years))
      processing_counts$n_missing_father_age <- processing_counts$n_missing_father_age + sum(is.na(resident_chunk$father_age_years))
      processing_counts$n_unknown_mother_educ2 <- processing_counts$n_unknown_mother_educ2 + sum(resident_chunk$mother_educ2 == "Unknown", na.rm = TRUE)
      processing_counts$n_unknown_father_educ2 <- processing_counts$n_unknown_father_educ2 + sum(resident_chunk$father_educ2 == "Unknown", na.rm = TRUE)
      processing_counts$n_missing_mother_race_eth <- processing_counts$n_missing_mother_race_eth + sum(is.na(resident_chunk$mother_race_eth))
      processing_counts$n_missing_father_race_eth <- processing_counts$n_missing_father_race_eth + sum(is.na(resident_chunk$father_race_eth))
      processing_counts$n_nonmissing_father_race_raw <- processing_counts$n_nonmissing_father_race_raw + sum(!is.na(resident_chunk$father_race_raw))
      processing_counts$n_nonmissing_father_educ_2003_raw <- processing_counts$n_nonmissing_father_educ_2003_raw + sum(!is.na(resident_chunk$father_educ_2003_raw))

      aggregated_chunks[[length(aggregated_chunks) + 1L]] <- aggregate_resident_births(resident_chunk)
    }
  }

  year_panel <- combine_aggregate_chunks(aggregated_chunks, get_birth_group_vars())

  if (!is.finite(observed_record_length_min)) observed_record_length_min <- NA_integer_
  if (!is.finite(observed_record_length_max)) observed_record_length_max <- NA_integer_

  year_metadata <- summarize_birth_processing(
    year = year,
    zip_path = zip_path,
    layout_name = layout$layout_name,
    min_record_length = layout$min_record_length,
    observed_record_length_min = observed_record_length_min,
    observed_record_length_max = observed_record_length_max,
    n_records_total = n_records_total,
    processing_counts = processing_counts
  )
  list(panel = year_panel, metadata = year_metadata)
}

# -- Fixed-width field positions and per-year layout registry --

pre2003_fields <- data.frame(
  field_name = c(
    "resident_status_raw", "mother_age_raw", "mother_hispanic_raw", "mother_race_eth_public_raw",
    "mother_race_raw", "mother_race_recode_raw", "mother_educ_1989_raw", "mother_educ_1989_recode_raw",
    "father_age_raw",
    "father_hispanic_raw", "father_race_eth_public_raw", "father_race_raw", "father_race_recode_raw",
    "father_educ_1989_raw"
  ),
  start = c(6L, 70L, 77L, 78L, 80L, 82L, 83L, 85L, 154L, 158L, 159L, 160L, 162L, 163L),
  end = c(6L, 71L, 77L, 78L, 81L, 82L, 84L, 85L, 155L, 158L, 159L, 161L, 162L, 164L),
  stringsAsFactors = FALSE
)

transition_fields <- data.frame(
  field_name = c(
    "resident_status_raw", "mother_age_raw",
    "mother_educ_2003_raw", "mother_educ_flag_raw",
    "mother_race_raw", "mother_hispanic_raw", "mother_race_eth_public_raw",
    "mother_educ_1989_raw", "mother_educ_1989_recode_raw",
    "father_age_raw", "father_age_alt_raw",
    "father_race_raw_revised", "father_hispanic_raw", "father_race_eth_public_raw", "father_race_raw_old",
    "father_educ_2003_raw", "father_educ_flag_raw"
  ),
  start = c(138L, 89L, 155L, 571L, 142L, 148L, 149L, 156L, 158L, 184L, 182L, 188L, 195L, 196L, 199L, 197L, 572L),
  end = c(138L, 90L, 155L, 571L, 143L, 148L, 149L, 157L, 158L, 185L, 183L, 189L, 195L, 196L, 200L, 197L, 572L),
  stringsAsFactors = FALSE
)

revised_fields <- data.frame(
  field_name = c(
    "resident_status_raw", "mother_age_raw", "mother_race_raw", "mother_race_recode6_raw",
    "mother_race_recode15_raw", "mother_hispanic_raw", "mother_hispanic_recode_raw",
    "mother_race_eth_public_raw", "mother_educ_2003_raw", "mother_educ_flag_raw",
    "father_age_raw",
    "father_race_raw", "father_race_recode6_raw", "father_race_recode15_raw",
    "father_hispanic_raw", "father_hispanic_recode_raw", "father_race_eth_public_raw",
    "father_educ_2003_raw", "father_educ_flag_raw"
  ),
  start = c(104L, 75L, 105L, 107L, 108L, 112L, 115L, 117L, 124L, 126L, 147L, 151L, 153L, 154L, 159L, 160L, 162L, 163L, 165L),
  end = c(104L, 76L, 106L, 107L, 109L, 112L, 115L, 117L, 124L, 126L, 148L, 152L, 153L, 155L, 159L, 160L, 162L, 163L, 165L),
  stringsAsFactors = FALSE
)

layout_registry <- c(
  setNames(
    replicate(
      length(1990:2002),
      list(layout_name = "1989_revision_pre2003", min_record_length = 164L, fields = pre2003_fields, harmonizer = harmonize_birth_chunk_pre2003),
      simplify = FALSE
    ),
    as.character(1990:2002)
  ),
  list(`2003` = list(layout_name = "transition_2003", min_record_length = 572L, fields = transition_fields, harmonizer = harmonize_birth_chunk_transition)),
  setNames(
    replicate(
      length(2004:2013),
      list(layout_name = "transition_2004_2013", min_record_length = 572L, fields = transition_fields, harmonizer = harmonize_birth_chunk_transition_direct_age),
      simplify = FALSE
    ),
    as.character(2004:2013)
  ),
  # 2014-2020 keep the bridged-compatible race/eth handling; 2021-2024 switch
  # to the single-race era (see harmonize_race_eth_2021plus).
  setNames(
    replicate(
      length(2014:2020),
      list(layout_name = "revised_2014_2020_bridged", min_record_length = 165L, fields = revised_fields, harmonizer = harmonize_birth_chunk_2016plus),
      simplify = FALSE
    ),
    as.character(2014:2020)
  ),
  setNames(
    replicate(
      length(2021:2024),
      list(layout_name = "revised_2021_2024_single_race", min_record_length = 165L, fields = revised_fields, harmonizer = harmonize_birth_chunk_2021plus),
      simplify = FALSE
    ),
    as.character(2021:2024)
  )
)

# ---- Process the natality files, year by year ----

for (output_dir in c(proc_dir, per_year_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

zip_files <- list.files(path = raw_dir, pattern = "\\.zip$", full.names = TRUE)
if (!length(zip_files)) stop("No natality zip files were found in: ", raw_dir, call. = FALSE)

zip_years <- vapply(zip_files, extract_year_from_file_name, integer(1), USE.NAMES = FALSE)
keep_idx <- !is.na(zip_years)
zip_files <- zip_files[keep_idx]
zip_years <- zip_years[keep_idx]

missing_years <- setdiff(as.integer(names(layout_registry)), unique(zip_years))
if (length(missing_years)) {
  stop(
    "The natality input is incomplete; missing years: ",
    paste(missing_years, collapse = ", "),
    call. = FALSE
  )
}
duplicate_years <- sort(unique(zip_years[duplicated(zip_years)]))
if (length(duplicate_years)) {
  stop(
    "More than one natality archive was found for: ",
    paste(duplicate_years, collapse = ", "),
    call. = FALSE
  )
}

ordered_idx <- order(zip_years, zip_files)
zip_files <- zip_files[ordered_idx]
zip_years <- zip_years[ordered_idx]

# Each year is processed independently and cached under per_year/. Setting
# NCHS_BIRTH_FORCE_REPROCESS=false reuses the years already on disk.
panel_results <- list()

for (i in seq_along(zip_files)) {
  year <- zip_years[i]
  layout <- layout_registry[[as.character(year)]]
  if (is.null(layout)) stop("No layout registry entry was found for year: ", year, call. = FALSE)

  cache_path <- file.path(per_year_dir, sprintf("01_birth_year_%s.rds", year))

  if (file.exists(cache_path) && !force_reprocess) {
    message("Using cached birth year ", year, ": ", basename(cache_path))
    processed <- readRDS(cache_path)
  } else {
    message("Processing birth file for year ", year, ": ", basename(zip_files[i]))
    processed <- process_birth_year_file(zip_files[i], year, layout, chunk_size)
    saveRDS(processed, cache_path)
  }

  panel_results[[length(panel_results) + 1L]] <- processed$panel
}

main_panel <- combine_aggregate_chunks(panel_results, get_birth_group_vars())

# ---- Birth panel ----

saveRDS(main_panel, main_panel_rds_path)
