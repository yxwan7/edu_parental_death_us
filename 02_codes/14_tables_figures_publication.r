# ==============================================================================
# 14_tables_figures_publication.r
#
# Build the reported tables and figures from the aggregated draws written by
# uncertainty/03_aggregate_draws.r.
# Reported values are draw medians with 95% confidence intervals; the
# deterministic values travel alongside as comparison columns.
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

here::i_am("02_codes/14_tables_figures_publication.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Paths ----

table_dir <- here::here("03_outputs", "01_tables")
figure_dir <- here::here("03_outputs", "02_figures")
appendix_table_dir <- here::here("03_outputs", "03_appendix", "tables")
appendix_figure_dir <- here::here("03_outputs", "03_appendix", "figures")
for (.dir in c(table_dir, figure_dir, appendix_table_dir, appendix_figure_dir)) {
  dir.create(.dir, recursive = TRUE, showWarnings = FALSE)
}

ci_root <- Sys.getenv("EDU_PARENT_LOSS_CI_ROOT", unset = "")
if (!nzchar(ci_root)) {
  stop("Set EDU_PARENT_LOSS_CI_ROOT to the aggregated bundle written by 15.",
       call. = FALSE)
}

batch_record <- readRDS(file.path(ci_root, "15_uncertainty_batch_record.rds"))
batch_summary <- batch_record$batch_summary
ALLOW_PARTIAL <- nzchar(Sys.getenv("EDU_PARENT_LOSS_ALLOW_PARTIAL", unset = ""))
if (!isTRUE(batch_summary$is_final) && !ALLOW_PARTIAL) {
  stop("The figures need a final 1..M bundle unless EDU_PARENT_LOSS_ALLOW_PARTIAL=1.",
       call. = FALSE)
}

read_ci <- function(name) as.data.table(readRDS(file.path(ci_root, name)))

incidence_total <- read_ci("15_ci_incidence_total.rds")
incidence_cause6 <- read_ci("15_ci_incidence_cause6.rds")
incidence_gap <- read_ci("15_ci_incidence_gap.rds")
meanage_total <- read_ci("15_ci_meanage_total.rds")
meanage_cause6 <- read_ci("15_ci_meanage_cause6.rds")
meanage_gap <- read_ci("15_ci_meanage_gap.rds")
main_number <- read_ci("15_ci_bereaved_children_number_table.rds")

# ---- Cause decompositions on the draw-median scale ----

# Every input is a draw median, so the bands and the total close exactly. That
# total differs slightly from value_median in the gap table -- a median is not
# additive -- so the marker and its interval are drawn as a separate layer.

incidence_decomposition <- local({
  wide <- dcast(
    incidence_cause6[parent_educ2 %chin% EDUC_LEVELS & cause_6 %chin% CAUSE_6_LEVELS],
    year + race_eth + parent_type + cause_6 ~ parent_educ2,
    value.var = "parental_loss_per_1000_median"
  )
  wide[, rd_contribution_per_1000 := No_BA - BA_plus]
  totals <- wide[, .(stacked_total_per_1000 = sum(rd_contribution_per_1000, na.rm = TRUE)),
                 by = .(year, race_eth, parent_type)]
  out <- merge(wide, totals, by = c("year", "race_eth", "parent_type"), sort = FALSE)
  out <- merge(
    out,
    incidence_gap[gap_metric == "RD", .(year, race_eth, parent_type,
                                        total_mid = value_median,
                                        total_lo = lower_ci,
                                        total_hi = upper_ci)],
    by = c("year", "race_eth", "parent_type"), sort = FALSE
  )
  setnames(out, c("No_BA", "BA_plus"),
           c("no_ba_rate_per_1000_median", "ba_plus_rate_per_1000_median"))
  out[, `:=`(race_order = match(race_eth, RACES_R3),
             parent_order = match(parent_type, PARENT_TYPE_LEVELS),
             cause_order = match(cause_6, CAUSE_6_LEVELS))]
  setorder(out, parent_order, race_order, year, cause_order)
  out[, c("race_order", "parent_order", "cause_order") := NULL]
  out[]
})

# The Shapley weights are event counts, which live in the incidence tables. Both
# come from the same aggregation, so taking the median counts alongside the
# median ages keeps the decomposition on one scale and closing exactly.
meanage_decomposition <- local({
  counts <- incidence_cause6[
    , .(year, race_eth, parent_educ2, parent_type, cause_6,
        n_children_losing_parent = n_children_losing_parent_median)
  ]
  ages <- meanage_cause6[
    , .(year, race_eth, parent_educ2, parent_type, cause_6,
        mean_age_at_loss = mean_age_at_loss_median)
  ]
  median_scale <- merge(
    ages, counts,
    by = c("year", "race_eth", "parent_educ2", "parent_type", "cause_6")
  )
  out <- build_meanage_decomposition(median_scale)
  merge(
    out,
    meanage_gap[, .(year, race_eth, parent_type,
                    total_mid = value_median,
                    total_lo = lower_ci,
                    total_hi = upper_ci)],
    by = c("year", "race_eth", "parent_type"), sort = FALSE
  )
})

# ---- Shared plotting pieces ----

TOTAL_KEY_LABEL <- "Total gaps and 95% CIs"
CAUSE_TAGS <- c("A", "B", "C", "D", "E", "F")

# One incidence level frame per reported value, for Figure 1 and Figure S1.
total_figure_data <- function(value_stem) {
  dt <- label_common(copy(incidence_total[parent_educ2 %chin% EDUC_LEVELS]))
  dt[, `:=`(
    mid = as.numeric(get(paste0(value_stem, "_median"))),
    lo = as.numeric(get(paste0(value_stem, "_lower_ci"))),
    hi = as.numeric(get(paste0(value_stem, "_upper_ci")))
  )]
  dt[]
}




# ---- Incidence figures ----

plot_cause6 <- function(parent_type_value, metric) {
  value_stem <- if (identical(metric, "rate")) {
    "parental_loss_per_1000"
  } else {
    "n_children_losing_parent"
  }
  dt <- label_common(copy(incidence_cause6[
    parent_type == parent_type_value &
      parent_educ2 %chin% EDUC_LEVELS &
      cause_6 %chin% CAUSE_6_LEVELS
  ]))
  dt[, `:=`(
    mid = as.numeric(get(paste0(value_stem, "_median"))),
    lo = as.numeric(get(paste0(value_stem, "_lower_ci"))),
    hi = as.numeric(get(paste0(value_stem, "_upper_ci")))
  )]
  dt <- dt[is.finite(mid)]

  if (identical(metric, "rate")) {
    y_label <- AXIS_RATE
    y_labeler <- label_number(accuracy = 0.01)
  } else {
    y_label <- AXIS_NUMBER
    y_labeler <- label_comma(accuracy = 1)
  }

  panels <- lapply(seq_along(CAUSE_6_LEVELS), function(i) {
    cause_value <- CAUSE_6_LEVELS[i]
    tag <- CAUSE_TAGS[i]
    is_last <- i > length(CAUSE_6_LEVELS) - 2L
    ggplot(dt[cause_6 == cause_value],
           aes(x = year, y = mid, color = parent_educ_label,
               fill = parent_educ_label, group = parent_educ_label)) +
      geom_ribbon(aes(ymin = lo, ymax = hi), alpha = RIBBON_ALPHA,
                  color = NA, show.legend = FALSE) +
      geom_line(linewidth = LINE_MAIN, lineend = "round") +
      geom_point(size = POINT_MAIN, stroke = 0) +
      facet_wrap(~race_label, nrow = 1) +
      common_x_scale() +
      scale_y_continuous(
        labels = y_labeler,
        expand = expansion(mult = c(0.05, if (identical(tag, "A")) 0.34 else 0.14))
      ) +
      educ_scales() +
      labs(title = panel_title(tag, CAUSE_6_LABELS[[cause_value]]),
           x = if (is_last) "Year" else NULL, y = NULL) +
      theme_parent_loss(axis_text_size = 7.8) +
      (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
  })

  compose_causes(panels, y_label)
}

# ---- Mean-age figures ----

plot_cause6_mean_age <- function(parent_type_value) {
  dt <- label_common(copy(meanage_cause6[
    parent_type == parent_type_value &
      parent_educ2 %chin% EDUC_LEVELS &
      cause_6 %chin% CAUSE_6_LEVELS
  ]))
  dt <- dt[is.finite(mean_age_at_loss_median)]

  panels <- lapply(seq_along(CAUSE_6_LEVELS), function(i) {
    cause_value <- CAUSE_6_LEVELS[i]
    tag <- CAUSE_TAGS[i]
    is_last <- i > length(CAUSE_6_LEVELS) - 2L
    ggplot(dt[cause_6 == cause_value],
           aes(x = year, y = mean_age_at_loss_median, color = parent_educ_label,
               fill = parent_educ_label, group = parent_educ_label)) +
      geom_ribbon(aes(ymin = mean_age_at_loss_lower_ci,
                      ymax = mean_age_at_loss_upper_ci),
                  alpha = RIBBON_ALPHA, color = NA, show.legend = FALSE) +
      geom_line(linewidth = LINE_MAIN, lineend = "round") +
      geom_point(size = POINT_MAIN, stroke = 0) +
      facet_wrap(~race_label, nrow = 1) +
      common_x_scale() +
      scale_y_continuous(
        labels = label_number(accuracy = 0.1),
        expand = expansion(mult = c(0.08, if (identical(tag, "A")) 0.36 else 0.16))
      ) +
      educ_scales() +
      labs(title = panel_title(tag, CAUSE_6_LABELS[[cause_value]]),
           x = if (is_last) "Year" else NULL, y = NULL) +
      theme_parent_loss(axis_text_size = 7.8) +
      (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
  })

  compose_causes(panels, AXIS_MEAN_AGE)
}

# ---- Write figures ----

# ---- Figure 1: incidence rate by race, education and parent ----

save_figure(
  plot_total(total_figure_data("parental_loss_per_1000"),
             function(w) sprintf("Incidence rate of children experiencing %s death, by %s education", w, w),
             AXIS_RATE, 0.1),
  file.path(figure_dir, "Fig_1_incidence_rates_by_race_education_parent.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 2: education gap in the incidence rate, by cause ----

save_figure(
  plot_incidence_decomposition(incidence_decomposition, TOTAL_KEY_LABEL),
  file.path(figure_dir, "Fig_2_incidence_rate_education_gap_by_race_parent_cause.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 3: mean age at loss by race, education and parent ----

meanage_total_figure <- label_common(copy(meanage_total[parent_educ2 %chin% EDUC_LEVELS]))
meanage_total_figure[, `:=`(
  mid = mean_age_at_loss_median,
  lo = mean_age_at_loss_lower_ci,
  hi = mean_age_at_loss_upper_ci
)]

save_figure(
  plot_total_mean_age(meanage_total_figure),
  file.path(figure_dir, "Fig_3_mean_age_by_race_education_parent.png"),
  width = 10.5, height = 8.4
)

# ---- Figure 4: education gap in mean age at loss, by cause ----

save_figure(
  plot_meanage_decomposition(meanage_decomposition, TOTAL_KEY_LABEL),
  file.path(figure_dir, "Fig_4_mean_age_education_gap_by_race_parent_cause.png"),
  width = 10.5, height = 8.4
)

# ---- Figure S1: number of bereaved children ----

save_figure(
  plot_total(total_figure_data("n_children_losing_parent"),
             function(w) sprintf("Number of children experiencing %s death, by %s education", w, w),
             AXIS_NUMBER, NULL),
  file.path(appendix_figure_dir, "Fig_S1_bereaved_children_number_by_race_education_parent.png"),
  width = 10.5, height = 8.4
)

# ---- Figure S2: maternal deaths by cause, number of children ----

save_figure(
  plot_cause6("mother", "count"),
  file.path(appendix_figure_dir, "Fig_S2_bereaved_children_number_maternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure S3: paternal deaths by cause, number of children ----

save_figure(
  plot_cause6("father", "count"),
  file.path(appendix_figure_dir, "Fig_S3_bereaved_children_number_paternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure S4: maternal deaths by cause, incidence rate ----

save_figure(
  plot_cause6("mother", "rate"),
  file.path(appendix_figure_dir, "Fig_S4_incidence_rates_maternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure S5: paternal deaths by cause, incidence rate ----

save_figure(
  plot_cause6("father", "rate"),
  file.path(appendix_figure_dir, "Fig_S5_incidence_rates_paternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure S6: maternal deaths by cause, mean age at loss ----

save_figure(
  plot_cause6_mean_age("mother"),
  file.path(appendix_figure_dir, "Fig_S6_mean_age_maternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure S7: paternal deaths by cause, mean age at loss ----

save_figure(
  plot_cause6_mean_age("father"),
  file.path(appendix_figure_dir, "Fig_S7_mean_age_paternal_by_race_education_cause.png"),
  width = 14, height = 10.5
)

# ---- Figure data ----

# The table each figure is drawn from. Several figures share one table, so a
# file can carry more rows than its figure plots. Aggregation diagnostics stay
# in the bundle.

figure_data_dir <- here::here("03_outputs", "04_figure_data")
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)

write_figure_data <- function(dt, stem) {
  diagnostics <- grep("_mcse_|_n_valid_draws$|_p_zero$", names(dt), value = TRUE)
  out <- if (length(diagnostics)) dt[, !diagnostics, with = FALSE] else copy(dt)
  fwrite(out, file.path(figure_data_dir, paste0(stem, "_data.csv")))
}

write_figure_data(incidence_total, "Fig_1_S1_incidence_by_race_education_parent")
write_figure_data(
  incidence_decomposition[, .(
    year, race_eth, parent_type, cause_6,
    no_ba_rate_per_1000_median, ba_plus_rate_per_1000_median,
    rd_contribution_per_1000, stacked_total_per_1000,
    total_gap_median = total_mid,
    total_gap_lower_ci = total_lo, total_gap_upper_ci = total_hi
  )],
  "Fig_2_incidence_rate_education_gap_by_race_parent_cause"
)
write_figure_data(meanage_total, "Fig_3_mean_age_by_race_education_parent")
write_figure_data(
  meanage_decomposition[, .(
    year, race_eth, parent_type, cause_6,
    share_No_BA, share_BA_plus, mean_age_No_BA, mean_age_BA_plus,
    shapley_contribution, stacked_total = mean_age_gap_total,
    total_gap_median = total_mid,
    total_gap_lower_ci = total_lo, total_gap_upper_ci = total_hi
  )],
  "Fig_4_mean_age_education_gap_by_race_parent_cause"
)
write_figure_data(incidence_cause6, "Fig_S2_S5_incidence_by_race_education_cause")
write_figure_data(meanage_cause6, "Fig_S6_S7_mean_age_by_race_education_cause")

# ---- Publication tables ----

draw_label <- if (isTRUE(batch_summary$is_final)) {
  sprintf("M = %d draws", batch_summary$m_complete)
} else {
  sprintf("Provisional: M = %d of %d draws",
          batch_summary$m_complete, batch_summary$m_expected)
}

# ---- Labels, formatting and table builders ----

PUBLICATION_KEY_YEARS <- c(2010L, 2019L, 2021L, 2024L)
PUBLICATION_EDUCATION_LABELS <- c(
  No_BA = "No BA", BA_plus = "BA or higher"
)

# Compact cause labels for the SI tables; CAUSE_6_LABELS keeps the fuller
# wording the figures have room for.
PUBLICATION_CAUSE_LABELS <- c(
  drug_poisoning = "Drug poisoning",
  firearm_deaths = "Firearm",
  covid19 = "COVID-19",
  cardiovascular_disease_and_diabetes = "CVD and diabetes",
  cancer = "Cancer",
  all_other_causes = "All other"
)
PUBLICATION_GAP_ROW_LABELS <- c(RD = "Gap", RR = "Ratio")

format_publication_ci <- function(value, lower, upper, digits = 2L) {
  rounding_threshold <- 0.5 * 10^(-digits)
  value[is.finite(value) & abs(value) < rounding_threshold] <- 0
  lower[is.finite(lower) & abs(lower) < rounding_threshold] <- 0
  upper[is.finite(upper) & abs(upper) < rounding_threshold] <- 0
  number_format <- paste0("%.", digits, "f")
  sprintf(
    paste0(number_format, " (", number_format, ", ", number_format, ")"),
    value, lower, upper
  )
}

format_publication_point <- function(value, digits = 2L) {
  rounding_threshold <- 0.5 * 10^(-digits)
  value[is.finite(value) & abs(value) < rounding_threshold] <- 0
  sprintf(paste0("%.", digits, "f"), value)
}

build_level_gap_publication_table <- function(
    level_dt, gap_dt, level_columns, gap_columns, digits = 2L,
    gap_label_column = NULL, gap_row_levels = "Gap") {
  levels <- as.data.table(level_dt)[
    year %in% PUBLICATION_KEY_YEARS &
      parent_educ2 %in% names(PUBLICATION_EDUCATION_LABELS),
    .(
      year, race_eth, parent_type, parent_educ2,
      estimate = get(level_columns[["estimate"]]),
      lower = get(level_columns[["lower"]]),
      upper = get(level_columns[["upper"]])
    )
  ]
  levels[, row_label := PUBLICATION_EDUCATION_LABELS[parent_educ2]]

  gaps <- as.data.table(gap_dt)[
    year %in% PUBLICATION_KEY_YEARS,
    .(
      year, race_eth, parent_type,
      estimate = get(gap_columns[["estimate"]]),
      lower = get(gap_columns[["lower"]]),
      upper = get(gap_columns[["upper"]]),
      row_label = if (is.null(gap_label_column)) "Gap" else
        as.character(get(gap_label_column))
    )
  ]
  gaps[, parent_educ2 := NA_character_]

  out <- rbindlist(list(levels, gaps), use.names = TRUE)
  row_levels <- c(unname(PUBLICATION_EDUCATION_LABELS), gap_row_levels)
  out[, `:=`(
    panel = PUBLICATION_PANEL_LABELS[parent_type],
    race = PUBLICATION_RACE_LABELS[race_eth],
    row_order = match(row_label, row_levels),
    race_order = match(race_eth, RACES_R3),
    parent_order = match(parent_type, PARENT_TYPE_LEVELS),
    display_value = format_publication_ci(estimate, lower, upper, digits)
  )]
  setorder(out, parent_order, race_order, row_order, year)

  display <- dcast(
    out,
    panel + race + row_label + parent_order + race_order + row_order ~ year,
    value.var = "display_value"
  )
  setorder(display, parent_order, race_order, row_order)
  setcolorder(display, c(
    "panel", "race", "row_label", "parent_order", "race_order", "row_order",
    as.character(PUBLICATION_KEY_YEARS)
  ))
  display[, c("parent_order", "race_order", "row_order") := NULL]
  display[]
}

build_decomposition_publication_table <- function(
    decomposition_dt, contribution_column, total_columns, digits = 2L) {
  dt <- as.data.table(decomposition_dt)[
    year %in% PUBLICATION_KEY_YEARS & cause_6 %in% CAUSE_6_LEVELS
  ]
  totals <- unique(dt[, .(
    year, race_eth, parent_type,
    estimate = get(total_columns[["estimate"]]),
    lower = get(total_columns[["lower"]]),
    upper = get(total_columns[["upper"]])
  )])
  totals[, `:=`(cause_6 = NA_character_, row_label = "Total gap")]

  contributions <- dt[, .(
    year, race_eth, parent_type, cause_6,
    estimate = get(contribution_column),
    lower = NA_real_,
    upper = NA_real_
  )]
  contributions[, row_label := PUBLICATION_CAUSE_LABELS[cause_6]]

  out <- rbindlist(list(totals, contributions), use.names = TRUE)
  row_levels <- c(
    "Total gap", unname(PUBLICATION_CAUSE_LABELS[CAUSE_6_LEVELS])
  )
  out[, `:=`(
    panel = PUBLICATION_PANEL_LABELS[parent_type],
    race = PUBLICATION_RACE_LABELS[race_eth],
    row_order = match(row_label, row_levels),
    race_order = match(race_eth, RACES_R3),
    parent_order = match(parent_type, PARENT_TYPE_LEVELS)
  )]
  out[, display_value := fifelse(
    row_label == "Total gap",
    format_publication_ci(estimate, lower, upper, digits),
    format_publication_point(estimate, digits)
  )]
  setorder(out, parent_order, race_order, row_order, year)

  display <- dcast(
    out,
    panel + race + row_label + parent_order + race_order + row_order ~ year,
    value.var = "display_value"
  )
  setorder(display, parent_order, race_order, row_order)
  setcolorder(display, c(
    "panel", "race", "row_label", "parent_order", "race_order", "row_order",
    as.character(PUBLICATION_KEY_YEARS)
  ))
  display[, c("parent_order", "race_order", "row_order") := NULL]
  display[]
}

write_publication_panel_table_tex <- function(
    display_dt, path, second_column_label, year_columns, panel_levels,
    race_levels, row_levels, notes, source = NULL, bold_rows = character(),
    vertically_center_rows = character(), thin_rule_after_rows = character(),
    panel_separator = "\\addlinespace", notes_size = "\\footnotesize") {
  display <- as.data.frame(display_dt, stringsAsFactors = FALSE)
  n_columns <- length(year_columns) + 2L
  tabular_columns <- paste0(
    "ll", paste(rep("l", length(year_columns)), collapse = "")
  )
  latex_years <- gsub("-", "--", year_columns, fixed = TRUE)
  row_ending <- "\\\\"
  latex_lines <- c(
    "\\begin{threeparttable}",
    sprintf("\\begin{tabular}{%s}", tabular_columns),
    "\\toprule",
    paste0(
      "Race and ethnicity & ", second_column_label, " & ",
      paste(latex_years, collapse = " & "), " ", row_ending
    ),
    "\\midrule"
  )

  for (panel_index in seq_along(panel_levels)) {
    panel_value <- panel_levels[[panel_index]]
    latex_lines <- c(
      latex_lines,
      paste0(
        sprintf(
          "\\multicolumn{%d}{l}{\\textbf{%s}}",
          n_columns, panel_value
        ),
        " ", row_ending
      )
    )
    for (race_index in seq_along(race_levels)) {
      race_value <- race_levels[[race_index]]
      rows <- display[
        display$panel == panel_value & display$race == race_value,
        , drop = FALSE
      ]
      rows <- rows[match(row_levels, rows$row_label), , drop = FALSE]
      for (row_index in seq_along(row_levels)) {
        row_value <- row_levels[[row_index]]
        row_label <- if (row_value %in% bold_rows) {
          sprintf("\\textbf{%s}", row_value)
        } else {
          row_value
        }
        if (row_value %in% vertically_center_rows) {
          row_label <- sprintf(
            "\\raisebox{0.5\\height}{%s}", row_label
          )
        }
        race_label <- if (row_index == 1L) {
          sprintf("\\multirow{%d}{*}{%s}", length(row_levels), race_value)
        } else {
          ""
        }
        values <- unlist(
          rows[row_index, year_columns, drop = FALSE],
          use.names = FALSE
        )
        values <- gsub(
          "%", "\\%", as.character(values), fixed = TRUE
        )
        latex_lines <- c(
          latex_lines,
          paste0(
            race_label, " & ", row_label, " & ",
            paste(values, collapse = " & "), " ", row_ending
          )
        )
        if (row_value %in% thin_rule_after_rows) {
          latex_lines <- c(
            latex_lines, sprintf("\\cmidrule(lr){2-%d}", n_columns)
          )
        }
      }
      if (race_index < length(race_levels)) {
        latex_lines <- c(
          latex_lines, sprintf("\\cmidrule(lr){1-%d}", n_columns)
        )
      }
    }
    if (panel_index < length(panel_levels)) {
      latex_lines <- c(latex_lines, panel_separator)
    }
  }

  latex_lines <- c(
    latex_lines,
    "\\bottomrule",
    "\\end{tabular}",
    "\\begin{tablenotes}[flushleft]",
    notes_size,
    paste0("\\item \\textbf{Notes:} ", paste(notes, collapse = "; ")),
    if (!is.null(source)) paste0("\\item \\textbf{Source:} ", source),
    "\\end{tablenotes}",
    "\\end{threeparttable}"
  )
  writeLines(latex_lines, path, useBytes = TRUE)
}

# ---- Table 1: bereaved children in key years ----

table1_period_levels <- c("2010", "2019", "2021", "2024", "2010-2024")
table1_row_levels <- c("Total", "No BA", "BA or higher")
table1_numeric <- copy(main_number)
table1_numeric[, panel := PUBLICATION_PANEL_LABELS[parent_type]]
table1_numeric[, race := PUBLICATION_RACE_LABELS[race_eth]]
table1_numeric[, row_label := factor(
  fcase(
    parent_educ2 == "Total", "Total",
    parent_educ2 == "No_BA", "No BA",
    parent_educ2 == "BA_plus", "BA or higher"
  ),
  levels = table1_row_levels
)]
table1_numeric[, parent_order := match(parent_type, PARENT_TYPE_LEVELS)]
table1_numeric[, race_order := match(race_eth, names(PUBLICATION_RACE_LABELS))]
table1_numeric[, period_order := match(period, table1_period_levels)]
setorder(table1_numeric, parent_order, race_order, row_label, period_order)

format_count_ci <- function(value, lower, upper) {
  format_count <- function(x) {
    comma(round(x), accuracy = 1)
  }
  sprintf(
    "\\shortstack[l]{%s \\\\ (%s, %s)}",
    format_count(value),
    format_count(lower),
    format_count(upper)
  )
}

format_count_ci_csv <- function(value, lower, upper) {
  format_count <- function(x) comma(round(x), accuracy = 1)
  sprintf(
    "%s (%s-%s)",
    format_count(value),
    format_count(lower),
    format_count(upper)
  )
}

table1_numeric[, display_value := fcase(
  parent_educ2 == "Total",
  format_count_ci(
    n_children_losing_parent_median,
    n_children_losing_parent_lower_ci,
    n_children_losing_parent_upper_ci
  ),
  default = sprintf("%.1f%%", share_percent_median)
)]

table1_numeric[, csv_value := fcase(
  parent_educ2 == "Total",
  format_count_ci_csv(
    n_children_losing_parent_median,
    n_children_losing_parent_lower_ci,
    n_children_losing_parent_upper_ci
  ),
  default = sprintf("%.1f%%", share_percent_median)
)]

table1_display <- dcast(
  table1_numeric,
  panel + race + row_label + parent_order + race_order ~ period,
  value.var = "display_value"
)
setorder(table1_display, parent_order, race_order, row_label)
setcolorder(table1_display, c(
  "panel", "race", "row_label", "parent_order", "race_order",
  table1_period_levels
))
table1_display[, c("parent_order", "race_order") := NULL]

table1_stem <- "Table_1_bereaved_children_number_key_years"
table1_numeric_output <- copy(table1_numeric)
table1_numeric_output[, c(
  "parent_order", "race_order", "period_order", "display_value", "csv_value"
) := NULL]
fwrite(
  table1_numeric_output,
  file.path(table_dir, paste0(table1_stem, "_data.csv"))
)
table1_csv <- dcast(
  table1_numeric,
  panel + race + row_label + parent_order + race_order ~ period,
  value.var = "csv_value"
)
setorder(table1_csv, parent_order, race_order, row_label)
setcolorder(table1_csv, c(
  "panel", "race", "row_label", "parent_order", "race_order",
  table1_period_levels
))
table1_csv[, c("parent_order", "race_order") := NULL]
setnames(table1_csv, "row_label", "parental_education")
fwrite(table1_csv, file.path(table_dir, paste0(table1_stem, ".csv")))
write_publication_panel_table_tex(
  display_dt = table1_display,
  path = file.path(table_dir, paste0(table1_stem, ".tex")),
  second_column_label = "Parental education",
  year_columns = table1_period_levels,
  panel_levels = unname(PUBLICATION_PANEL_LABELS[PARENT_TYPE_LEVELS]),
  race_levels = unname(PUBLICATION_RACE_LABELS),
  row_levels = table1_row_levels,
  bold_rows = "Total",
  vertically_center_rows = "Total",
  thin_rule_after_rows = "Total",
  panel_separator = "\\midrule",
  notes_size = "\\small",
  notes = paste0(
    "BA = bachelor's degree; CI = confidence interval; NH = non-Hispanic; ",
    "the 2010--2024 column is the sum of annual numbers over the 15-year ",
    "period; values in parentheses are 95\\% CIs from 1,000 Monte Carlo ",
    "replicates taking the 2.5th and 97.5th percentiles. No BA and BA or ",
    "higher rows are the percentages of each total attributed to the ",
    "corresponding parental education group, and percentages may not sum ",
    "to 100 because of rounding. All groups rows pool the three racial and ",
    "ethnic groups, with the pooling applied within each replicate."
  )
)

# ---- Shared pieces for the SI tables ----

publication_year_columns <- as.character(PUBLICATION_KEY_YEARS)
publication_panel_levels <- unname(
  PUBLICATION_PANEL_LABELS[PARENT_TYPE_LEVELS]
)
publication_race_levels <- unname(PUBLICATION_RACE_LABELS[RACES_R3])
publication_gap_rows <- unname(PUBLICATION_GAP_ROW_LABELS)

# ---- Table S4: incidence rates and education gaps ----

incidence_gap_rows <- incidence_gap[
  gap_metric %in% names(PUBLICATION_GAP_ROW_LABELS)
]
incidence_gap_rows[, gap_metric := PUBLICATION_GAP_ROW_LABELS[gap_metric]]

table_s4 <- build_level_gap_publication_table(
  level_dt = incidence_total,
  gap_dt = incidence_gap_rows,
  level_columns = c(
    estimate = "parental_loss_per_1000_median",
    lower = "parental_loss_per_1000_lower_ci",
    upper = "parental_loss_per_1000_upper_ci"
  ),
  gap_columns = c(
    estimate = "value_median", lower = "lower_ci", upper = "upper_ci"
  ),
  gap_label_column = "gap_metric",
  gap_row_levels = publication_gap_rows
)
table_s4_csv <- copy(table_s4)
setnames(
  table_s4_csv,
  c("race", "row_label"),
  c("race_ethnicity", "parental_education")
)
table_s4_stem <- "Table_S4_incidence_rates_education_gaps_key_years"
fwrite(
  table_s4_csv,
  file.path(appendix_table_dir, paste0(table_s4_stem, ".csv"))
)
write_publication_panel_table_tex(
  display_dt = table_s4,
  path = file.path(appendix_table_dir, paste0(table_s4_stem, ".tex")),
  second_column_label = "Parental education",
  year_columns = publication_year_columns,
  panel_levels = publication_panel_levels,
  race_levels = publication_race_levels,
  row_levels = c(unname(PUBLICATION_EDUCATION_LABELS), publication_gap_rows),
  bold_rows = publication_gap_rows,
  panel_separator = "\\midrule",
  notes_size = "\\small",
  notes = c(
    "BA = bachelor's degree",
    "CI = confidence interval",
    "NH = non-Hispanic",
    "Gap = No BA $-$ BA or higher",
    "Ratio = No BA / BA or higher",
    paste0(
      "values in parentheses are 95\\% CIs from 1,000 Monte Carlo replicates ",
      "taking the 2.5th and 97.5th percentiles."
    )
  )
)

# ---- Table S5: mean age and education gaps ----

table_s5 <- build_level_gap_publication_table(
  level_dt = meanage_total,
  gap_dt = meanage_gap,
  level_columns = c(
    estimate = "mean_age_at_loss_median",
    lower = "mean_age_at_loss_lower_ci",
    upper = "mean_age_at_loss_upper_ci"
  ),
  gap_columns = c(
    estimate = "value_median", lower = "lower_ci", upper = "upper_ci"
  )
)
table_s5_csv <- copy(table_s5)
setnames(
  table_s5_csv,
  c("race", "row_label"),
  c("race_ethnicity", "parental_education")
)
table_s5_stem <- "Table_S5_mean_age_education_gaps_key_years"
fwrite(
  table_s5_csv,
  file.path(appendix_table_dir, paste0(table_s5_stem, ".csv"))
)
write_publication_panel_table_tex(
  display_dt = table_s5,
  path = file.path(appendix_table_dir, paste0(table_s5_stem, ".tex")),
  second_column_label = "Parental education",
  year_columns = publication_year_columns,
  panel_levels = publication_panel_levels,
  race_levels = publication_race_levels,
  row_levels = c(unname(PUBLICATION_EDUCATION_LABELS), "Gap"),
  bold_rows = "Gap",
  panel_separator = "\\midrule",
  notes_size = "\\small",
  notes = c(
    "BA = bachelor's degree",
    "CI = confidence interval",
    "NH = non-Hispanic",
    "Gap = No BA $-$ BA or higher",
    paste0(
      "values in parentheses are 95\\% CIs from 1,000 Monte Carlo replicates ",
      "taking the 2.5th and 97.5th percentiles."
    )
  )
)

# ---- Table S6: incidence gap cause decomposition ----

table_s6 <- build_decomposition_publication_table(
  decomposition_dt = incidence_decomposition,
  contribution_column = "rd_contribution_per_1000",
  total_columns = c(
    estimate = "total_mid",
    lower = "total_lo",
    upper = "total_hi"
  )
)
table_s6_csv <- copy(table_s6)
setnames(
  table_s6_csv,
  c("race", "row_label"),
  c("race_ethnicity", "cause_of_death")
)
table_s6_stem <- "Table_S6_incidence_gap_cause_decomposition_key_years"
fwrite(
  table_s6_csv,
  file.path(appendix_table_dir, paste0(table_s6_stem, ".csv"))
)
write_publication_panel_table_tex(
  display_dt = table_s6,
  path = file.path(appendix_table_dir, paste0(table_s6_stem, ".tex")),
  second_column_label = "Cause of death",
  year_columns = publication_year_columns,
  panel_levels = publication_panel_levels,
  race_levels = publication_race_levels,
  row_levels = c(
    "Total gap", unname(PUBLICATION_CAUSE_LABELS[CAUSE_6_LEVELS])
  ),
  bold_rows = "Total gap",
  panel_separator = "\\midrule",
  notes_size = "\\small",
  notes = c(
    "BA = bachelor's degree",
    "CI = confidence interval",
    "NH = non-Hispanic",
    "Total gap = No BA $-$ BA or higher for all causes",
    paste0(
      "values in parentheses are 95\\% CIs from 1,000 Monte Carlo replicates ",
      "taking the 2.5th and 97.5th percentiles."
    )
  )
)

# ---- Table S7: mean-age gap cause decomposition ----

table_s7 <- build_decomposition_publication_table(
  decomposition_dt = meanage_decomposition,
  contribution_column = "shapley_contribution",
  total_columns = c(
    estimate = "total_mid",
    lower = "total_lo",
    upper = "total_hi"
  )
)
table_s7_csv <- copy(table_s7)
setnames(
  table_s7_csv,
  c("race", "row_label"),
  c("race_ethnicity", "cause_of_death")
)
table_s7_stem <- "Table_S7_mean_age_gap_cause_decomposition_key_years"
fwrite(
  table_s7_csv,
  file.path(appendix_table_dir, paste0(table_s7_stem, ".csv"))
)
write_publication_panel_table_tex(
  display_dt = table_s7,
  path = file.path(appendix_table_dir, paste0(table_s7_stem, ".tex")),
  second_column_label = "Cause of death",
  year_columns = publication_year_columns,
  panel_levels = publication_panel_levels,
  race_levels = publication_race_levels,
  row_levels = c(
    "Total gap", unname(PUBLICATION_CAUSE_LABELS[CAUSE_6_LEVELS])
  ),
  bold_rows = "Total gap",
  panel_separator = "\\midrule",
  notes_size = "\\small",
  notes = c(
    "BA = bachelor's degree",
    "CI = confidence interval",
    "NH = non-Hispanic",
    "Total gap = No BA $-$ BA or higher for all causes",
    paste0(
      "values in parentheses are 95\\% CIs from 1,000 Monte Carlo replicates ",
      "taking the 2.5th and 97.5th percentiles."
    )
  )
)

cat(sprintf("Tables and figures written for batch %s\n", batch_summary$batch_id))
