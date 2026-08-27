# ==============================================================================
# 12_point_estimates.r
#
# Build the deterministic incidence and mean-age tables from 11_parental_death.r.
# They are the comparison values that uncertainty/03_aggregate_draws.r joins
# alongside the draw medians, so they have to exist before the uncertainty layer
# runs.
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

# Both endpoints come off the cause_atom panels: the all-cause total is the
# cause_6 partition summed, so there is no separate all-cause input.
incidence_input_path <- file.path(
  processed_dir, "02_parent_loss_by_cause_atom_annual.rds"
)
meanage_input_path <- file.path(
  processed_dir, "03_parent_loss_by_cause_atom_child_age.rds"
)

write_point_estimate <- function(dt, stem) {
  # Written straight after the table is built. A subset in between would leave a
  # data.table index on the object and change the bytes, which the batch record
  # in record_run_environment.r checks when a draw is resubmitted.
  saveRDS(dt, file.path(processed_dir, paste0("point_estimate_", stem, ".rds")))
}

# ---- Incidence ----

incidence_raw <- readRDS(incidence_input_path)
incidence_total <- order_total(build_total_incidence(incidence_raw, EDUC_TABLE_LEVELS))
incidence_cause6 <- order_cause(
  build_cause_incidence(incidence_raw, EDUC_TABLE_LEVELS, "cause_6", CAUSE_6_LEVELS),
  "cause_6", CAUSE_6_LEVELS
)

write_point_estimate(incidence_total, "incidence_total")
write_point_estimate(incidence_cause6, "incidence_cause6")

# ---- Mean age at loss ----

meanage_raw <- readRDS(meanage_input_path)
meanage_total <- order_total(build_total_mean_age(meanage_raw, EDUC_TABLE_LEVELS))
meanage_cause6 <- order_cause(
  build_cause_mean_age(meanage_raw, EDUC_TABLE_LEVELS, "cause_6", CAUSE_6_LEVELS),
  "cause_6", CAUSE_6_LEVELS
)

write_point_estimate(meanage_total, "meanage_total")
write_point_estimate(meanage_cause6, "meanage_cause6")
