# ==============================================================================
# 16_tables_figures.r
#
# Build the reported tables and figures from the aggregated draws written by 15.
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

here::i_am("02_codes/16_tables_figures.r")

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

# ==============================================================================
# Cause decompositions on the draw-median scale
# ==============================================================================

# Every input is a draw median, so the coloured bands and the total come from
# one set of numbers and close exactly. That total therefore differs slightly
# from value_median in the gap table -- a median is not additive -- which is why
# the marker and its interval are drawn as a separate layer.

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
                                        total_rd_median = value_median,
                                        total_rd_lower_ci = lower_ci,
                                        total_rd_upper_ci = upper_ci)],
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
                    total_gap_median = value_median,
                    total_gap_lower_ci = lower_ci,
                    total_gap_upper_ci = upper_ci)],
    by = c("year", "race_eth", "parent_type"), sort = FALSE
  )
})

# ==============================================================================
# Shared plotting pieces
# ==============================================================================

RIBBON_ALPHA <- 0.16
LINE_MAIN <- 0.75
POINT_MAIN <- 2.1
TOTAL_KEY_LABEL <- "Total gaps and 95% CIs"
CAUSE_TAGS <- c("A", "B", "C", "D", "E", "F")

label_common <- function(dt) {
  dt[, `:=`(
    race_label = factor(RACE_LABELS[race_eth], levels = RACE_LABELS[RACES_R3]),
    parent_type_label = factor(PARENT_LABELS[parent_type],
                               levels = PARENT_LABELS[PARENT_TYPE_LEVELS])
  )]
  if ("parent_educ2" %chin% names(dt)) {
    dt[, parent_educ_label := factor(EDUC_LABELS[parent_educ2],
                                     levels = EDUC_LABELS[EDUC_LEVELS])]
  }
  if ("cause_6" %chin% names(dt)) {
    dt[, cause_6_label := factor(CAUSE_6_LABELS[cause_6],
                                 levels = CAUSE_6_LABELS[CAUSE_6_LEVELS])]
  }
  dt[]
}

educ_scales <- list(
  scale_color_manual(
    values = setNames(EDUC_PALETTE[EDUC_LEVELS], EDUC_LABELS[EDUC_LEVELS]),
    name = "Parent education"
  ),
  scale_fill_manual(
    values = setNames(EDUC_PALETTE[EDUC_LEVELS], EDUC_LABELS[EDUC_LEVELS]),
    guide = "none"
  )
)

as_stack_factor <- function(dt) {
  dt[, cause_6_label := factor(CAUSE_6_LABELS[cause_6],
                               levels = CAUSE_6_LABELS[CAUSE_6_STACK_LEVELS])]
  dt[]
}

# ==============================================================================
# Incidence figures
# ==============================================================================

incidence_total_panel <- function(dt, parent_type_value, tag, description,
                                  y_breaks, y_labeler, y_label, show_x_title) {
  ggplot(dt[parent_type == parent_type_value],
         aes(x = year, y = mid, color = parent_educ_label,
             fill = parent_educ_label, group = parent_educ_label)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = RIBBON_ALPHA,
                color = NA, show.legend = FALSE) +
    geom_line(linewidth = LINE_MAIN, lineend = "round") +
    geom_point(size = POINT_MAIN, stroke = 0) +
    facet_wrap(~race_label, nrow = 1) +
    common_x_scale() +
    scale_y_continuous(breaks = y_breaks, limits = c(0, max(y_breaks)),
                       labels = y_labeler, expand = expansion(mult = c(0, 0.05))) +
    educ_scales +
    labs(title = panel_title(tag, description),
         x = if (show_x_title) "Year" else NULL, y = y_label) +
    theme_parent_loss(axis_text_size = 7.8) +
    (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
}

plot_total <- function(value_stem, describe, y_label, y_accuracy) {
  dt <- label_common(copy(incidence_total[parent_educ2 %chin% EDUC_LEVELS]))
  dt[, `:=`(
    mid = as.numeric(get(paste0(value_stem, "_median"))),
    lo = as.numeric(get(paste0(value_stem, "_lower_ci"))),
    hi = as.numeric(get(paste0(value_stem, "_upper_ci")))
  )]
  y_breaks <- make_zero_based_axis_breaks(dt[, max(hi, na.rm = TRUE)])
  y_labeler <- if (is.null(y_accuracy)) {
    label_number(big.mark = ",")
  } else {
    label_number(accuracy = y_accuracy, big.mark = ",")
  }
  compose_ab(
    incidence_total_panel(dt, "mother", "A", describe("maternal"),
                          y_breaks, y_labeler, y_label, FALSE),
    incidence_total_panel(dt, "father", "B", describe("paternal"),
                          y_breaks, y_labeler, y_label, TRUE)
  )
}

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
      educ_scales +
      labs(title = panel_title(tag, CAUSE_6_LABELS[[cause_value]]),
           x = if (is_last) "Year" else NULL, y = NULL) +
      theme_parent_loss(axis_text_size = 7.8) +
      (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
  })

  compose_causes(panels, y_label)
}

incidence_decomposition_panel <- function(part_dt, total_dt, tag, description, y_axis,
                                          show_x_title, legend_here) {
  ggplot(part_dt, aes(x = year, y = rd_contribution_per_1000, fill = cause_6_label)) +
    geom_col(width = 0.82, colour = NA) +
    geom_linerange(data = total_dt,
                   aes(x = year, ymin = total_rd_lower_ci, ymax = total_rd_upper_ci,
                       colour = TOTAL_KEY_LABEL),
                   inherit.aes = FALSE, linewidth = 0.45) +
    geom_point(data = total_dt,
               aes(x = year, y = total_rd_median, colour = TOTAL_KEY_LABEL),
               inherit.aes = FALSE, shape = 16, size = 1.5) +
    scale_colour_manual(values = setNames("grey15", TOTAL_KEY_LABEL), name = NULL) +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    facet_wrap(~race_label, nrow = 1) +
    common_x_scale(pad = 1.2) +
    scale_y_continuous(breaks = y_axis$breaks, limits = y_axis$limits,
                       labels = label_number(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      values = setNames(CAUSE_6_PALETTE[CAUSE_6_STACK_LEVELS],
                        CAUSE_6_LABELS[CAUSE_6_STACK_LEVELS]),
      name = "Cause of death"
    ) +
    labs(title = panel_title(tag, description),
         x = if (show_x_title) "Year" else NULL, y = AXIS_RATE_DIFFERENCE) +
    theme_parent_loss(axis_text_size = 7.8) +
    (if (legend_here) {
      inside_legend("top") + paired_legend_box()
    } else {
      theme(legend.position = "none")
    })
}

plot_incidence_decomposition <- function() {
  parts <- as_stack_factor(label_common(copy(incidence_decomposition)))
  totals <- label_common(unique(incidence_decomposition[, .(
    year, race_eth, parent_type,
    total_rd_median, total_rd_lower_ci, total_rd_upper_ci
  )]))

  stack_bounds <- parts[, .(
    positive = sum(pmax(rd_contribution_per_1000, 0), na.rm = TRUE),
    negative = sum(pmin(rd_contribution_per_1000, 0), na.rm = TRUE)
  ), by = .(year, race_eth, parent_type)]
  y_axis <- make_stacked_axis(
    min(stack_bounds$negative, totals$total_rd_lower_ci, na.rm = TRUE),
    max(stack_bounds$positive, totals$total_rd_upper_ci, na.rm = TRUE)
  )

  describe <- function(word) {
    sprintf("Cause contributions to the educational gap in the incidence rate of children experiencing %s death",
            word)
  }
  build <- function(parent_type_value, tag, is_last) {
    incidence_decomposition_panel(
      parts[parent_type == parent_type_value],
      totals[parent_type == parent_type_value],
      tag, describe(PARENT_WORDS[[parent_type_value]]), y_axis, is_last,
      identical(parent_type_value, "mother")
    )
  }
  compose_ab(build("mother", "A", FALSE), build("father", "B", TRUE))
}

# ==============================================================================
# Mean-age figures
# ==============================================================================

mean_age_panel <- function(dt, parent_type_value, tag, description, y_limits,
                           show_x_title) {
  ggplot(dt[parent_type == parent_type_value],
         aes(x = year, y = mean_age_at_loss_median, color = parent_educ_label,
             fill = parent_educ_label, group = parent_educ_label)) +
    geom_ribbon(aes(ymin = mean_age_at_loss_lower_ci, ymax = mean_age_at_loss_upper_ci),
                alpha = RIBBON_ALPHA, color = NA, show.legend = FALSE) +
    geom_line(linewidth = LINE_MAIN, lineend = "round") +
    geom_point(size = POINT_MAIN, stroke = 0) +
    facet_wrap(~race_label, nrow = 1) +
    common_x_scale() +
    scale_y_continuous(limits = y_limits,
                       breaks = seq(y_limits[1], y_limits[2], by = 0.5),
                       labels = label_number(accuracy = 0.1),
                       expand = expansion(mult = c(0.02, 0.02))) +
    educ_scales +
    labs(title = panel_title(tag, description),
         x = if (show_x_title) "Year" else NULL, y = AXIS_MEAN_AGE) +
    theme_parent_loss(axis_text_size = 7.8) +
    (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
}

plot_total_mean_age <- function() {
  dt <- label_common(copy(meanage_total[parent_educ2 %chin% EDUC_LEVELS]))
  dt <- dt[is.finite(mean_age_at_loss_median)]
  y_limits <- c(
    min(9, floor(dt[, min(mean_age_at_loss_lower_ci, na.rm = TRUE)] * 2) / 2),
    max(12, ceiling(dt[, max(mean_age_at_loss_upper_ci, na.rm = TRUE)] * 2) / 2)
  )
  describe <- function(word) {
    sprintf("Mean age of children experiencing %s death, by %s education", word, word)
  }
  compose_ab(
    mean_age_panel(dt, "mother", "A", describe("maternal"), y_limits, FALSE),
    mean_age_panel(dt, "father", "B", describe("paternal"), y_limits, TRUE)
  )
}

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
      educ_scales +
      labs(title = panel_title(tag, CAUSE_6_LABELS[[cause_value]]),
           x = if (is_last) "Year" else NULL, y = NULL) +
      theme_parent_loss(axis_text_size = 7.8) +
      (if (identical(tag, "A")) inside_legend("top") else theme(legend.position = "none"))
  })

  compose_causes(panels, AXIS_MEAN_AGE)
}

meanage_decomposition_panel <- function(part_dt, total_dt, tag, description, y_axis,
                                        show_x_title, legend_here) {
  ggplot(part_dt, aes(x = year, y = shapley_contribution, fill = cause_6_label)) +
    geom_col(width = 0.82, colour = NA) +
    geom_linerange(data = total_dt,
                   aes(x = year, ymin = total_gap_lower_ci, ymax = total_gap_upper_ci,
                       colour = TOTAL_KEY_LABEL),
                   inherit.aes = FALSE, linewidth = 0.45) +
    geom_point(data = total_dt,
               aes(x = year, y = total_gap_median, colour = TOTAL_KEY_LABEL),
               inherit.aes = FALSE, shape = 16, size = 1.5) +
    scale_colour_manual(values = setNames("grey15", TOTAL_KEY_LABEL), name = NULL) +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2)) +
    facet_wrap(~race_label, nrow = 1) +
    common_x_scale(pad = 1.2) +
    scale_y_continuous(breaks = y_axis$breaks, limits = y_axis$limits,
                       labels = label_number(accuracy = 0.1),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      values = setNames(CAUSE_6_PALETTE[CAUSE_6_STACK_LEVELS],
                        CAUSE_6_LABELS[CAUSE_6_STACK_LEVELS]),
      name = "Cause of death"
    ) +
    labs(title = panel_title(tag, description),
         x = if (show_x_title) "Year" else NULL, y = AXIS_MEAN_AGE_DIFFERENCE) +
    theme_parent_loss(axis_text_size = 7.8) +
    (if (legend_here) {
      inside_legend("bottom") + paired_legend_box()
    } else {
      theme(legend.position = "none")
    })
}

plot_meanage_decomposition <- function() {
  parts <- as_stack_factor(label_common(copy(meanage_decomposition)))
  totals <- label_common(unique(meanage_decomposition[, .(
    year, race_eth, parent_type,
    total_gap_median, total_gap_lower_ci, total_gap_upper_ci
  )]))

  stack_bounds <- parts[, .(
    positive = sum(pmax(shapley_contribution, 0), na.rm = TRUE),
    negative = sum(pmin(shapley_contribution, 0), na.rm = TRUE)
  ), by = .(year, race_eth, parent_type)]
  y_axis <- make_stacked_axis(
    min(stack_bounds$negative, totals$total_gap_lower_ci, na.rm = TRUE),
    max(stack_bounds$positive, totals$total_gap_upper_ci, na.rm = TRUE),
    target_intervals = 6L, slack = 1.5
  )

  describe <- function(word) {
    sprintf("Cause contributions to the educational gap in mean age of children experiencing %s death",
            word)
  }
  build <- function(parent_type_value, tag, is_last) {
    meanage_decomposition_panel(
      parts[parent_type == parent_type_value],
      totals[parent_type == parent_type_value],
      tag, describe(PARENT_WORDS[[parent_type_value]]), y_axis, is_last,
      identical(parent_type_value, "mother")
    )
  }
  compose_ab(build("mother", "A", FALSE), build("father", "B", TRUE))
}

# ==============================================================================
# Write figures
# ==============================================================================

figures <- list(
  list(figure_dir, "Fig_1_incidence_rates_by_race_education_parent.png", 10.5, 8.4,
       plot_total("parental_loss_per_1000",
                  function(w) sprintf("Incidence rate of children experiencing %s death, by %s education", w, w),
                  AXIS_RATE, 0.1)),
  list(figure_dir, "Fig_2_incidence_rate_education_gap_by_race_parent_cause.png", 10.5, 8.4,
       plot_incidence_decomposition()),
  list(figure_dir, "Fig_3_mean_age_by_race_education_parent.png", 10.5, 8.4,
       plot_total_mean_age()),
  list(figure_dir, "Fig_4_mean_age_education_gap_by_race_parent_cause.png", 10.5, 8.4,
       plot_meanage_decomposition()),
  list(appendix_figure_dir, "Fig_S1_bereaved_children_number_by_race_education_parent.png", 10.5, 8.4,
       plot_total("n_children_losing_parent",
                  function(w) sprintf("Number of children experiencing %s death, by %s education", w, w),
                  AXIS_NUMBER, NULL)),
  list(appendix_figure_dir, "Fig_S2_bereaved_children_number_maternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6("mother", "count")),
  list(appendix_figure_dir, "Fig_S3_bereaved_children_number_paternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6("father", "count")),
  list(appendix_figure_dir, "Fig_S4_incidence_rates_maternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6("mother", "rate")),
  list(appendix_figure_dir, "Fig_S5_incidence_rates_paternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6("father", "rate")),
  list(appendix_figure_dir, "Fig_S6_mean_age_maternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6_mean_age("mother")),
  list(appendix_figure_dir, "Fig_S7_mean_age_paternal_by_race_education_cause.png", 14, 10.5,
       plot_cause6_mean_age("father"))
)

for (entry in figures) {
  save_figure(entry[[5]], file.path(entry[[1]], entry[[2]]),
              width = entry[[3]], height = entry[[4]])
}

# ==============================================================================
# Publication tables
# ==============================================================================

draw_label <- if (isTRUE(batch_summary$is_final)) {
  sprintf("M = %d draws", batch_summary$m_complete)
} else {
  sprintf("Provisional: M = %d of %d draws",
          batch_summary$m_complete, batch_summary$m_expected)
}

PUBLICATION_KEY_YEARS <- c(2010L, 2019L, 2021L, 2024L)
PUBLICATION_RACE_LABELS <- c(
  white = "NH White", black = "NH Black", hispanic = "Hispanic"
)
PUBLICATION_PANEL_LABELS <- c(
  mother = "Panel A: maternal death",
  father = "Panel B: paternal death"
)
PUBLICATION_EDUCATION_LABELS <- c(
  No_BA = "No BA", BA_plus = "BA or higher"
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
  contributions[, row_label := CAUSE_6_LABELS[cause_6]]

  out <- rbindlist(list(totals, contributions), use.names = TRUE)
  row_levels <- c(
    "Total gap", unname(CAUSE_6_LABELS[CAUSE_6_LEVELS])
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
table1_numeric[, race_order := match(race_eth, RACES_R3)]
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
  race_levels = unname(PUBLICATION_RACE_LABELS[RACES_R3]),
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
    "to 100 because of rounding."
  )
)

publication_year_columns <- as.character(PUBLICATION_KEY_YEARS)
publication_panel_levels <- unname(
  PUBLICATION_PANEL_LABELS[PARENT_TYPE_LEVELS]
)
publication_race_levels <- unname(PUBLICATION_RACE_LABELS[RACES_R3])
publication_gap_rows <- unname(PUBLICATION_GAP_ROW_LABELS)

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

table_s6 <- build_decomposition_publication_table(
  decomposition_dt = incidence_decomposition,
  contribution_column = "rd_contribution_per_1000",
  total_columns = c(
    estimate = "total_rd_median",
    lower = "total_rd_lower_ci",
    upper = "total_rd_upper_ci"
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
    "Total gap", unname(CAUSE_6_LABELS[CAUSE_6_LEVELS])
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

table_s7 <- build_decomposition_publication_table(
  decomposition_dt = meanage_decomposition,
  contribution_column = "shapley_contribution",
  total_columns = c(
    estimate = "total_gap_median",
    lower = "total_gap_lower_ci",
    upper = "total_gap_upper_ci"
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
    "Total gap", unname(CAUSE_6_LABELS[CAUSE_6_LEVELS])
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

cat(sprintf("Tables and %d figures written for batch %s\n",
            length(figures), batch_summary$batch_id))
