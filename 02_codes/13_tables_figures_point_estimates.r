# ==============================================================================
# 13_tables_figures_point_estimates.r
#
# Build the main-text figures and table from the deterministic endpoints in
# 12_point_estimates.r. These need no draw batch, so they run without a cluster.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(scales)
  library(patchwork)
  library(ragg)
})

here::i_am("02_codes/13_tables_figures_point_estimates.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Paths ----

processed_dir <- here::here("01_data", "02_processed", "edu_parent_loss")
output_dir <- here::here("03_outputs", "00_point_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_point <- function(name) as.data.table(readRDS(file.path(processed_dir, name)))

incidence_total <- read_point("point_estimate_incidence_total.rds")
incidence_cause6 <- read_point("point_estimate_incidence_cause6.rds")
meanage_total <- read_point("point_estimate_meanage_total.rds")
meanage_cause6 <- read_point("point_estimate_meanage_cause6.rds")

TOTAL_KEY_LABEL <- "Total gaps"

# ---- Cause decompositions ----

# Both decompositions close exactly on their total here, so the total comes off
# the contributions themselves. The marker still carries information a stack
# does not: with contributions on both sides of zero, the net gap is not the
# top of the stack.

incidence_decomposition <- local({
  wide <- dcast(
    incidence_cause6[parent_educ2 %chin% EDUC_LEVELS & cause_6 %chin% CAUSE_6_LEVELS],
    year + race_eth + parent_type + cause_6 ~ parent_educ2,
    value.var = "parental_loss_per_1000"
  )
  wide[, rd_contribution_per_1000 := No_BA - BA_plus]
  totals <- wide[, .(total_mid = sum(rd_contribution_per_1000, na.rm = TRUE)),
                 by = .(year, race_eth, parent_type)]
  out <- merge(wide, totals, by = c("year", "race_eth", "parent_type"), sort = FALSE)
  setnames(out, c("No_BA", "BA_plus"),
           c("no_ba_rate_per_1000", "ba_plus_rate_per_1000"))
  out[, `:=`(race_order = match(race_eth, RACES_R3),
             parent_order = match(parent_type, PARENT_TYPE_LEVELS),
             cause_order = match(cause_6, CAUSE_6_LEVELS))]
  setorder(out, parent_order, race_order, year, cause_order)
  out[, c("race_order", "parent_order", "cause_order") := NULL]
  out[]
})

meanage_decomposition <- local({
  frame <- merge(
    meanage_cause6[, .(year, race_eth, parent_educ2, parent_type, cause_6,
                       mean_age_at_loss)],
    incidence_cause6[, .(year, race_eth, parent_educ2, parent_type, cause_6,
                         n_children_losing_parent)],
    by = c("year", "race_eth", "parent_educ2", "parent_type", "cause_6")
  )
  out <- build_meanage_decomposition(frame)
  out[, total_mid := mean_age_gap_total]
  out[]
})

# ---- Write figures ----

# ---- Figure 1: incidence rate by race, education and parent ----

incidence_total_figure <- label_common(copy(incidence_total[parent_educ2 %chin% EDUC_LEVELS]))
incidence_total_figure[, mid := as.numeric(parental_loss_per_1000)]

save_figure(
  plot_total(incidence_total_figure,
             function(w) sprintf("Incidence rate of children experiencing %s death, by %s education", w, w),
             AXIS_RATE, 0.1),
  file.path(output_dir, "Fig_1_incidence_rates_by_race_education_parent.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 2: education gap in the incidence rate, by cause ----

save_figure(
  plot_incidence_decomposition(incidence_decomposition, TOTAL_KEY_LABEL),
  file.path(output_dir, "Fig_2_incidence_rate_education_gap_by_race_parent_cause.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 3: mean age at loss by race, education and parent ----

meanage_total_figure <- label_common(copy(meanage_total[parent_educ2 %chin% EDUC_LEVELS]))
meanage_total_figure[, mid := as.numeric(mean_age_at_loss)]

save_figure(
  plot_total_mean_age(meanage_total_figure),
  file.path(output_dir, "Fig_3_mean_age_by_race_education_parent.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 4: education gap in mean age at loss, by cause ----

save_figure(
  plot_meanage_decomposition(meanage_decomposition, TOTAL_KEY_LABEL),
  file.path(output_dir, "Fig_4_mean_age_education_gap_by_race_parent_cause.png"),
  width = 10.5, height = 8.4
)

# ---- Figure data ----

# The table each figure is drawn from. The level tables also carry the pooled
# Total education row, which the figures do not plot.

write_figure_data <- function(dt, stem) {
  fwrite(dt, file.path(output_dir, paste0(stem, "_data.csv")))
}

write_figure_data(incidence_total, "Fig_1_incidence_rates_by_race_education_parent")
write_figure_data(
  incidence_decomposition[, .(
    year, race_eth, parent_type, cause_6,
    no_ba_rate_per_1000, ba_plus_rate_per_1000,
    rd_contribution_per_1000, total_gap_per_1000 = total_mid
  )],
  "Fig_2_incidence_rate_education_gap_by_race_parent_cause"
)
write_figure_data(meanage_total, "Fig_3_mean_age_by_race_education_parent")
write_figure_data(
  meanage_decomposition[, .(
    year, race_eth, parent_type, cause_6,
    share_No_BA, share_BA_plus, mean_age_No_BA, mean_age_BA_plus,
    shapley_contribution, total_gap = mean_age_gap_total
  )],
  "Fig_4_mean_age_education_gap_by_race_parent_cause"
)

# ---- Table 1: bereaved children in key years ----

table1_period_levels <- c("2010", "2019", "2021", "2024", MAIN_TEXT_PERIOD_LABEL)
table1_row_levels <- c("Total", "No BA", "BA or higher")

table1 <- build_bereaved_children_number_table(incidence_total)
table1[, panel := PUBLICATION_PANEL_LABELS[parent_type]]
table1[, race := PUBLICATION_RACE_LABELS[race_eth]]
table1[, row_label := factor(
  fcase(
    parent_educ2 == "Total", "Total",
    parent_educ2 == "No_BA", "No BA",
    parent_educ2 == "BA_plus", "BA or higher"
  ),
  levels = table1_row_levels
)]
table1[, csv_value := fcase(
  parent_educ2 == "Total", comma(round(n_children_losing_parent), accuracy = 1),
  default = sprintf("%.1f%%", share_percent)
)]
table1[, `:=`(
  parent_order = match(parent_type, PARENT_TYPE_LEVELS),
  race_order = match(race_eth, names(PUBLICATION_RACE_LABELS))
)]

table1_csv <- dcast(
  table1,
  panel + race + row_label + parent_order + race_order ~ period,
  value.var = "csv_value"
)
setorder(table1_csv, parent_order, race_order, row_label)
setcolorder(table1_csv, c(
  "panel", "race", "row_label", "parent_order", "race_order", table1_period_levels
))
table1_csv[, c("parent_order", "race_order") := NULL]
setnames(table1_csv, "row_label", "parental_education")

fwrite(
  table1_csv,
  file.path(output_dir, "Table_1_bereaved_children_number_key_years.csv")
)

cat("Point-estimate table and figures written\n")
