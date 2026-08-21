# ==============================================================================
# 01_download_nchs_files.r
#
# Download public-use NCHS natality and mortality files and record the download
# inventory.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(here)
  library(rvest)
})

here::i_am("02_codes/01_download_nchs_files.r")

# ---- Download settings ----

PORTAL_URL <- "https://www.cdc.gov/nchs/data_access/vitalstatsonline.htm"

start_year <- 1990L
end_year <- NA_integer_
overwrite_existing <- FALSE
run_download <- FALSE
max_retries <- 3L
retry_sleep_seconds <- 2L
use_inventory_fallback <- TRUE
validate_downloads <- TRUE

# NCHS public-use zips can be large, so raise R's default download timeout.
options(timeout = max(600, getOption("timeout")))

# ---- Download directories ----

raw_dir      <- here("01_data", "01_raw")
nchs_dir     <- file.path(raw_dir, "nchs")
birth_dir    <- file.path(nchs_dir, "natality")
death_dir    <- file.path(nchs_dir, "mortality")
metadata_dir <- file.path(nchs_dir, "metadata")

# ---- Portal and download functions ----

get_env_integer <- function(var_name, default_value = NA_integer_) {
  raw_value <- Sys.getenv(var_name, unset = "")

  if (identical(raw_value, "")) {
    return(default_value)
  }

  parsed_value <- suppressWarnings(as.integer(raw_value))

  if (is.na(parsed_value)) {
    stop("Environment variable `", var_name, "` must be an integer.", call. = FALSE)
  }

  parsed_value
}

get_env_flag <- function(var_name, default_value = FALSE) {
  raw_value <- tolower(trimws(Sys.getenv(var_name, unset = "")))

  if (identical(raw_value, "")) {
    return(default_value)
  }

  if (raw_value %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }

  if (raw_value %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }

  stop(
    "Environment variable `",
    var_name,
    "` must be one of TRUE/FALSE/1/0/yes/no.",
    call. = FALSE
  )
}

write_download_csv <- function(df, output_path) {
  utils::write.csv(
    x = df,
    file = output_path,
    row.names = FALSE
  )
}

normalize_href <- function(href) {
  href <- trimws(href)
  href[href == ""] <- NA_character_

  is_relative <- !is.na(href) & !grepl("^(https?|ftp)://", href)
  href[is_relative] <- paste0("https://www.cdc.gov", href[is_relative])

  href
}

get_output_dir_for_data_type <- function(data_type) {
  if (identical(data_type, "natality")) {
    return(birth_dir)
  }

  if (identical(data_type, "mortality")) {
    return(death_dir)
  }

  stop("Unsupported `data_type` in inventory: ", data_type, call. = FALSE)
}

normalize_inventory <- function(inventory) {
  required_columns <- c("data_type", "year", "href")
  missing_columns <- setdiff(required_columns, names(inventory))

  if (length(missing_columns) > 0L) {
    stop(
      "Inventory is missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  inventory$data_type <- trimws(as.character(inventory$data_type))
  inventory$year <- as.integer(inventory$year)
  inventory$href <- trimws(as.character(inventory$href))

  if (!"file_name" %in% names(inventory)) {
    inventory$file_name <- basename(inventory$href)
  } else {
    inventory$file_name <- trimws(as.character(inventory$file_name))
    missing_file_name <- is.na(inventory$file_name) | inventory$file_name == ""
    inventory$file_name[missing_file_name] <- basename(inventory$href[missing_file_name])
  }

  output_dirs <- vapply(
    inventory$data_type,
    get_output_dir_for_data_type,
    character(1)
  )

  inventory$dest_path <- file.path(output_dirs, inventory$file_name)
  inventory <- inventory[, c("data_type", "year", "file_name", "href", "dest_path"), drop = FALSE]
  inventory <- inventory[order(inventory$data_type, inventory$year), , drop = FALSE]
  rownames(inventory) <- NULL

  inventory
}

get_portal_links <- function(url) {
  portal_html <- rvest::read_html(url)
  link_nodes  <- rvest::html_elements(portal_html, "a[href]")
  link_text   <- trimws(rvest::html_text2(link_nodes))
  link_href   <- normalize_href(rvest::html_attr(link_nodes, "href"))

  link_table <- data.frame(
    label = link_text,
    href  = link_href,
    stringsAsFactors = FALSE
  )

  link_table <- link_table[!is.na(link_table$href), , drop = FALSE]
  rownames(link_table) <- NULL

  link_table
}

build_series_inventory <- function(
  link_table,
  data_type,
  start_year = 1990L,
  end_year = NA_integer_
) {
  if (identical(data_type, "natality")) {
    # Historical natality files use mixed naming conventions, e.g. `Nat1990.zip`
    # for early years and `Nat2024us.zip` for later public-use releases.
    file_pattern <- "^nat\\d{4}(us)?\\.zip$"
  } else if (identical(data_type, "mortality")) {
    file_pattern <- "^mort\\d{4}us\\.zip$"
  } else {
    stop("`data_type` must be either 'natality' or 'mortality'.", call. = FALSE)
  }

  output_dir <- get_output_dir_for_data_type(data_type)

  file_names <- basename(link_table$href)

  inventory <- link_table[
    grepl(file_pattern, file_names, ignore.case = TRUE),
    ,
    drop = FALSE
  ]

  inventory <- inventory[!duplicated(inventory$href), , drop = FALSE]
  inventory$year <- as.integer(sub(".*?(\\d{4}).*", "\\1", basename(inventory$href)))
  inventory <- inventory[!is.na(inventory$year), , drop = FALSE]
  inventory <- inventory[inventory$year >= start_year, , drop = FALSE]

  if (!is.na(end_year)) {
    inventory <- inventory[inventory$year <= end_year, , drop = FALSE]
  }

  inventory <- inventory[order(inventory$year), , drop = FALSE]

  if (nrow(inventory) == 0L) {
    stop(
      "No ",
      data_type,
      " files were found on the CDC portal for years >= ",
      start_year,
      ".",
      call. = FALSE
    )
  }

  inventory$data_type <- data_type
  inventory$file_name <- basename(inventory$href)
  inventory$dest_path <- file.path(output_dir, inventory$file_name)

  inventory <- inventory[
    ,
    c("data_type", "year", "file_name", "href", "dest_path"),
    drop = FALSE
  ]

  rownames(inventory) <- NULL

  inventory
}

build_download_inventory <- function(start_year = 1990L, end_year = NA_integer_) {
  portal_links <- get_portal_links(PORTAL_URL)

  inventory <- rbind(
    build_series_inventory(
      link_table = portal_links,
      data_type  = "natality",
      start_year = start_year,
      end_year   = end_year
    ),
    build_series_inventory(
      link_table = portal_links,
      data_type  = "mortality",
      start_year = start_year,
      end_year   = end_year
    )
  )

  normalize_inventory(inventory)
}

build_download_inventory_or_cached <- function(
  start_year,
  end_year,
  inventory_path,
  use_inventory_fallback = TRUE
) {
  inventory <- tryCatch(
    build_download_inventory(start_year = start_year, end_year = end_year),
    error = function(e) {
      message("Live portal inventory build failed: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(inventory)) {
    return(inventory)
  }

  if (isTRUE(use_inventory_fallback) && file.exists(inventory_path)) {
    message("Falling back to saved inventory: ", inventory_path)

    inventory <- utils::read.csv(
      file = inventory_path,
      stringsAsFactors = FALSE
    )

    inventory <- inventory[inventory$year >= start_year, , drop = FALSE]

    if (!is.na(end_year)) {
      inventory <- inventory[inventory$year <= end_year, , drop = FALSE]
    }

    if (nrow(inventory) == 0L) {
      stop("Saved inventory exists but has no rows in the requested year range.", call. = FALSE)
    }

    return(normalize_inventory(inventory))
  }

  stop(
    "Failed to build the live inventory and no usable fallback inventory was available.",
    call. = FALSE
  )
}

validate_downloaded_zip <- function(dest_path) {
  if (!file.exists(dest_path)) {
    stop("Downloaded file does not exist: ", dest_path, call. = FALSE)
  }

  file_size <- file.info(dest_path)$size

  if (is.na(file_size) || file_size <= 0) {
    stop("Downloaded file is empty: ", dest_path, call. = FALSE)
  }

  zip_listing <- tryCatch(
    utils::unzip(dest_path, list = TRUE),
    error = function(e) {
      stop("Downloaded zip failed validation: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (nrow(zip_listing) == 0L) {
    stop("Downloaded zip contains no members: ", dest_path, call. = FALSE)
  }

  TRUE
}

download_one_file <- function(
  file_url,
  dest_path,
  overwrite_existing = FALSE,
  max_retries = 3L,
  retry_sleep_seconds = 2L,
  validate_downloads = TRUE
) {
  result <- list(
    status = NA_character_,
    attempts = 0L,
    error_message = NA_character_,
    file_size_bytes = NA_real_,
    validated_zip = NA
  )

  if (file.exists(dest_path) && !overwrite_existing) {
    if (isTRUE(validate_downloads)) {
      existing_file_is_valid <- tryCatch(
        {
          validate_downloaded_zip(dest_path)
          TRUE
        },
        error = function(e) {
          result$error_message <<- conditionMessage(e)
          FALSE
        }
      )

      if (!isTRUE(existing_file_is_valid)) {
        message(
          "Existing file failed validation and will be re-downloaded: ",
          basename(dest_path)
        )

        overwrite_existing <- TRUE
      }
    }

    if (!isTRUE(overwrite_existing)) {
      result$status <- "skipped_existing"
      result$file_size_bytes <- unname(file.info(dest_path)$size)
      result$validated_zip <- isTRUE(validate_downloads)
      return(result)
    }
  }

  temp_path <- paste0(dest_path, ".part")

  if (file.exists(temp_path)) {
    file.remove(temp_path)
  }

  for (attempt in seq_len(max_retries)) {
    result$attempts <- attempt

    if (file.exists(temp_path)) {
      file.remove(temp_path)
    }

    attempt_result <- tryCatch(
      {
        utils::download.file(
          url      = file_url,
          destfile = temp_path,
          mode     = "wb",
          method   = "libcurl",
          quiet    = FALSE
        )

        if (file.exists(dest_path) && overwrite_existing) {
          file.remove(dest_path)
        }

        rename_success <- file.rename(from = temp_path, to = dest_path)

        if (!isTRUE(rename_success)) {
          stop("Downloaded file could not be moved into place: ", dest_path, call. = FALSE)
        }

        if (isTRUE(validate_downloads)) {
          validate_downloaded_zip(dest_path)
        }

        result$status <- "downloaded"
        result$file_size_bytes <- unname(file.info(dest_path)$size)
        result$validated_zip <- isTRUE(validate_downloads)
        result$error_message <- NA_character_

        result
      },
      error = function(e) {
        if (file.exists(temp_path)) {
          file.remove(temp_path)
        }

        list(
          status = "failed",
          attempts = attempt,
          error_message = conditionMessage(e),
          file_size_bytes = if (file.exists(dest_path)) unname(file.info(dest_path)$size) else NA_real_,
          validated_zip = FALSE
        )
      }
    )

    if (identical(attempt_result$status, "downloaded")) {
      return(attempt_result)
    }

    result <- attempt_result

    if (attempt < max_retries) {
      message(
        "Attempt ",
        attempt,
        " failed for ",
        basename(dest_path),
        ". Retrying after ",
        retry_sleep_seconds * attempt,
        " seconds."
      )

      Sys.sleep(retry_sleep_seconds * attempt)
    }
  }

  result
}

download_nchs_files <- function(
  inventory,
  download_log_path,
  overwrite_existing = FALSE,
  max_retries = 3L,
  retry_sleep_seconds = 2L,
  validate_downloads = TRUE
) {
  results <- inventory
  results$status <- "pending"
  results$attempts <- 0L
  results$error_message <- NA_character_
  results$file_size_bytes <- NA_real_
  results$validated_zip <- NA
  results$download_timestamp_utc <- NA_character_

  write_download_csv(results, download_log_path)

  for (i in seq_len(nrow(results))) {
    message(
      "Processing ",
      results$data_type[i],
      " ",
      results$year[i],
      ": ",
      results$file_name[i]
    )

    attempt_result <- tryCatch(
      download_one_file(
        file_url = results$href[i],
        dest_path = results$dest_path[i],
        overwrite_existing = overwrite_existing,
        max_retries = max_retries,
        retry_sleep_seconds = retry_sleep_seconds,
        validate_downloads = validate_downloads
      ),
      error = function(e) {
        list(
          status = "failed",
          attempts = max_retries,
          error_message = conditionMessage(e),
          file_size_bytes = if (file.exists(results$dest_path[i])) {
            unname(file.info(results$dest_path[i])$size)
          } else {
            NA_real_
          },
          validated_zip = FALSE
        )
      }
    )

    results$status[i] <- attempt_result$status
    results$attempts[i] <- attempt_result$attempts
    results$error_message[i] <- attempt_result$error_message
    results$file_size_bytes[i] <- attempt_result$file_size_bytes
    results$validated_zip[i] <- attempt_result$validated_zip
    results$download_timestamp_utc[i] <- format(
      Sys.time(),
      tz = "UTC",
      usetz = TRUE
    )

    write_download_csv(results, download_log_path)
  }

  results
}

# ---- Build the download inventory ----

# Settings can be overridden per run via environment variables; the inventory is
# rebuilt from the live portal each run so new years appear without hard-coding.
start_year <- get_env_integer(
  var_name = "NCHS_DOWNLOAD_START_YEAR",
  default_value = start_year
)
end_year <- get_env_integer(
  var_name = "NCHS_DOWNLOAD_END_YEAR",
  default_value = end_year
)
overwrite_existing <- get_env_flag(
  var_name = "NCHS_OVERWRITE_EXISTING",
  default_value = overwrite_existing
)
run_download <- get_env_flag(
  var_name = "NCHS_RUN_DOWNLOAD",
  default_value = run_download
)
max_retries <- get_env_integer(
  var_name = "NCHS_MAX_RETRIES",
  default_value = max_retries
)
retry_sleep_seconds <- get_env_integer(
  var_name = "NCHS_RETRY_SLEEP_SECONDS",
  default_value = retry_sleep_seconds
)
use_inventory_fallback <- get_env_flag(
  var_name = "NCHS_USE_INVENTORY_FALLBACK",
  default_value = use_inventory_fallback
)
validate_downloads <- get_env_flag(
  var_name = "NCHS_VALIDATE_DOWNLOADS",
  default_value = validate_downloads
)

for (output_dir in c(raw_dir, nchs_dir, birth_dir, death_dir, metadata_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

inventory_path <- file.path(metadata_dir, "00_nchs_download_inventory.csv")
download_log_path <- file.path(metadata_dir, "00_nchs_download_log.csv")

download_inventory <- build_download_inventory_or_cached(
  start_year = start_year,
  end_year = end_year,
  inventory_path = inventory_path,
  use_inventory_fallback = use_inventory_fallback
)

write_download_csv(download_inventory, inventory_path)

latest_birth_year <- max(
  download_inventory$year[download_inventory$data_type == "natality"]
)
latest_death_year <- max(
  download_inventory$year[download_inventory$data_type == "mortality"]
)

message(
  "Inventory: natality through ", latest_birth_year,
  ", mortality through ", latest_death_year,
  " (", nrow(download_inventory), " files from ", start_year, ")."
)

# ---- Download raw files ----

if (isTRUE(run_download)) {
  download_log <- download_nchs_files(
    inventory = download_inventory,
    download_log_path = download_log_path,
    overwrite_existing = overwrite_existing,
    max_retries = max_retries,
    retry_sleep_seconds = retry_sleep_seconds,
    validate_downloads = validate_downloads
  )

  n_downloaded <- sum(download_log$status == "downloaded")
  n_skipped <- sum(download_log$status == "skipped_existing")
  n_failed <- sum(download_log$status == "failed")

  if (n_failed > 0L) {
    warning(
      n_failed,
      " file(s) failed. Progress was preserved in the download log and can be resumed.",
      call. = FALSE
    )
  }
} else {
  message("Download skipped. Set NCHS_RUN_DOWNLOAD=true to fetch the zips.")
}
