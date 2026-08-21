# ==============================================================================
# 12_point_estimates.r
#
# Build the deterministic incidence and mean-age tables from 11. They are the
# comparison values that 15 joins alongside the draw medians, so they have to
# exist before the uncertainty layer runs.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

here::i_am("02_codes/12_point_estimates.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Paths ----

processed_dir <- here::here("01_data", "02_processed", "edu_parent_loss")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

incidence_total_input_path <- file.path(
  processed_dir, "03b_parent_loss_by_3_cause_annual.rds"
)
incidence_cause6_input_path <- file.path(
  processed_dir, "04b_parent_loss_by_6_15_cause_annual.rds"
)
meanage_total_input_path <- file.path(
  processed_dir, "03c_parent_loss_by_3_cause_child_age.rds"
)
meanage_cause6_input_path <- file.path(
  processed_dir, "04c_parent_loss_by_6_15_cause_child_age.rds"
)

write_point_estimate <- function(dt, stem) {
  # Written straight after the table is built. A subset in between would leave a
  # data.table index on the object and change the bytes, which the batch record
  # in record_run_environment.r checks when a draw is resubmitted.
  saveRDS(dt, file.path(processed_dir, paste0("point_estimate_", stem, ".rds")))
}

# ---- Incidence ----

incidence_total <- order_total(
  build_total_incidence(readRDS(incidence_total_input_path), EDUC_TABLE_LEVELS)
)
incidence_cause6 <- order_cause6(
  build_cause6_incidence(readRDS(incidence_cause6_input_path), EDUC_TABLE_LEVELS)
)

write_point_estimate(incidence_total, "incidence_total")
write_point_estimate(incidence_cause6, "incidence_cause6")

# ---- Mean age at loss ----

meanage_total <- order_total(
  build_total_mean_age(readRDS(meanage_total_input_path), EDUC_TABLE_LEVELS)
)
meanage_cause6 <- order_cause6(
  build_cause6_mean_age(readRDS(meanage_cause6_input_path), EDUC_TABLE_LEVELS)
)

write_point_estimate(meanage_total, "meanage_total")
write_point_estimate(meanage_cause6, "meanage_cause6")
