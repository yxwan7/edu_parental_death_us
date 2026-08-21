# ==============================================================================
# 01_male_fertility.r
#
# Summarize annual father-side RD/RR for the S1 and S2 sensitivity analyses
# and generate Appendix Figures S9 and S10.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(patchwork)
  library(scales)
})

here::i_am("02_codes/sensitivity/01_male_fertility.r")
source(here::here("02_codes", "00_functions.r"))

# ---- Sensitivity legs and files ----

leg_specs <- as.data.table(copy(MALE_SENSITIVITY_LEGS))
annual_input_name <- "03b_parent_loss_by_3_cause_annual.rds"

processed_root <- here(
  "01_data", "02_processed", "edu_parent_loss", "sensitivity"
)
output_dir <- here("03_outputs", "03_appendix", "sensitivity")
figure_dir <- here("03_outputs", "03_appendix", "figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

s1_annual_path <- file.path(
  output_dir,
  "s1_father_side_rd_rr_by_scenario_annual.csv"
)
s2_annual_path <- file.path(
  output_dir,
  "s2_father_side_rd_rr_by_scenario_annual.csv"
)
s1_figure_path <- file.path(
  figure_dir,
  "Fig_S9_male_fertility_gradient_sensitivity.png"
)
s2_figure_path <- file.path(
  figure_dir,
  "Fig_S10_father_education_missingness_sensitivity.png"
)

# ---- Read annual sensitivity incidence ----

read_leg_incidence <- function(leg, directory) {
  input_path <- file.path(processed_root, directory, annual_input_name)
  dt <- as.data.table(readRDS(input_path))
  out <- dt[
    ,
    .(
      n_children_losing_parent = sum(n_children_losing_parent),
      n_children_with_alive_parent_under18 =
        first(n_children_with_alive_parent_under18)
    ),
    by = .(year, race_eth, parent_educ2, parent_type)
  ]
  out[
    ,
    parental_loss_per_1000 := 1000 *
      n_children_losing_parent / n_children_with_alive_parent_under18
  ]
  out[, leg := leg]
  out[]
}

# ---- Annual RD and RR tables ----

annual_incidence <- rbindlist(lapply(seq_len(nrow(leg_specs)), function(i) {
  read_leg_incidence(leg_specs$leg[i], leg_specs$directory[i])
}))

annual_gap <- dcast(
  annual_incidence[
    year %in% YEARS_REPORT &
      parent_type == "father" &
      parent_educ2 %chin% EDUC_LEVELS,
    .(year, race_eth, leg, parent_educ2, parental_loss_per_1000)
  ],
  year + race_eth + leg ~ parent_educ2,
  value.var = "parental_loss_per_1000"
)
annual_gap[, `:=`(
  rd_per_1000 = No_BA - BA_plus,
  rr = safe_divide(No_BA, BA_plus)
)]

s1_annual_specs <- data.table(
  leg = c("primary", "female_borrow", "zero_gradient"),
  scenario = c("S1", "S2", "S3"),
  scenario_label = c(
    "S1: Primary",
    "S2: Female gradient",
    "S3: No gradient"
  ),
  definition = c(
    "Father-observed education gradient; missing education allocated under MAR",
    "Male proxy borrows the female education gradient",
    "No male fertility education gradient"
  )
)

s2_annual_specs <- data.table(
  leg = c("primary", "complete_case", "mnar_all_ba", "mnar_all_noba"),
  scenario = c("S1", "S2", "S3", "S4"),
  scenario_label = c(
    "S1: Primary",
    "S2: Complete case",
    "S3: All missing to BA",
    "S4: All missing to No BA"
  ),
  definition = c(
    "Missing paternal education allocated under MAR",
    "Births missing paternal education are excluded",
    "All missing paternal education is assigned to BA+",
    "All missing paternal education is assigned to No BA"
  )
)

build_annual_gap_table <- function(specs, analysis_name) {
  out <- merge(
    copy(annual_gap[leg %chin% specs$leg]),
    specs,
    by = "leg",
    all.x = TRUE,
    sort = FALSE
  )
  out <- melt(
    out,
    id.vars = c(
      "year", "race_eth", "leg", "scenario", "scenario_label", "definition",
      "BA_plus", "No_BA"
    ),
    measure.vars = c("rd_per_1000", "rr"),
    variable.name = "gap_metric",
    value.name = "value",
    variable.factor = FALSE
  )
  out[, analysis := analysis_name]
  setorder(out, race_eth, gap_metric, scenario, year)
  out[
    ,
    .(
      analysis, year, race_eth, scenario, scenario_label, leg, definition,
      gap_metric,
      BA_plus_rate_per_1000 = BA_plus,
      No_BA_rate_per_1000 = No_BA,
      value
    )
  ]
}

s1_annual_gap <- build_annual_gap_table(
  s1_annual_specs,
  "male_fertility_education_gradient"
)
s2_annual_gap <- build_annual_gap_table(
  s2_annual_specs,
  "father_education_missingness"
)

# ---- Export annual RD and RR tables ----

fwrite(s1_annual_gap, s1_annual_path)
fwrite(s2_annual_gap, s2_annual_path)

# ---- Annual RD and RR figures ----

s1_plot_specs <- s1_annual_specs[, .(scenario_label)]
s1_plot_specs[, colour := JAMA_COLOURS[seq_len(.N)]]
s2_plot_specs <- s2_annual_specs[, .(scenario_label)]
s2_plot_specs[, colour := JAMA_COLOURS[seq_len(.N)]]

prepare_figure_data <- function(path, specs) {
  out <- fread(path)
  out[, scenario_label := factor(
    scenario_label,
    levels = specs$scenario_label
  )]
  out[, gap_metric_label := factor(
    fcase(
      gap_metric == "rd_per_1000", "Rate difference",
      gap_metric == "rr", "Rate ratio"
    ),
    levels = c("Rate difference", "Rate ratio")
  )]
  setorder(out, race_eth, gap_metric_label, scenario_label, year)
  out[]
}

build_race_panel <- function(
    dt, specs, race_value, tag, show_x_title, show_legend, legend_rows) {
  palette <- setNames(specs$colour, specs$scenario_label)

  out <- ggplot(
    dt[race_eth == race_value],
    aes(
      x = year,
      y = value,
      colour = scenario_label,
      group = scenario_label
    )
  ) +
    geom_line(linewidth = 0.72, lineend = "round") +
    geom_point(size = 1.65, stroke = 0) +
    facet_wrap(~gap_metric_label, nrow = 1, scales = "free_y") +
    common_x_scale() +
    scale_y_continuous(
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_colour_manual(
      values = palette,
      breaks = specs$scenario_label,
      name = NULL
    ) +
    labs(
      title = panel_title(tag, RACE_LABELS[[race_value]]),
      x = if (show_x_title) "Year" else NULL,
      y = "Estimate"
    ) +
    theme_parent_loss(
      title_size = 10.5,
      strip_text_size = 8.8,
      axis_text_size = 7.8,
      panel_spacing_x = 0.8,
      panel_spacing_y = 0.5,
      legend_key_width_cm = 0.8
    ) +
    theme(
      plot.margin = margin(
        t = if (identical(tag, "A")) 8 else 0,
        r = 0,
        b = 0,
        l = 0
      )
    )

  if (show_legend) {
    out +
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size = 8.4),
        legend.key.width = grid::unit(0.8, "cm"),
        legend.key.height = grid::unit(0.4, "cm"),
        legend.margin = margin(5, 0, 0, 0)
      ) +
      guides(colour = guide_legend(nrow = legend_rows, byrow = TRUE))
  } else {
    out + theme(legend.position = "none")
  }
}

build_figure <- function(dt, specs, legend_rows) {
  (
    build_race_panel(dt, specs, "white", "A", FALSE, FALSE, legend_rows) /
      build_race_panel(dt, specs, "black", "B", FALSE, FALSE, legend_rows) /
      build_race_panel(dt, specs, "hispanic", "C", TRUE, TRUE, legend_rows)
  ) +
    plot_layout(heights = c(1, 1, 1.12))
}

save_sensitivity_figure <- function(plot_obj, output_path) {
  warmup_path <- tempfile(fileext = ".png")
  on.exit(unlink(warmup_path), add = TRUE)
  save_figure(plot_obj, warmup_path, width = 10.5, height = 10.8)
  save_figure(plot_obj, output_path, width = 10.5, height = 10.8)
}

s1_figure_data <- prepare_figure_data(s1_annual_path, s1_plot_specs)
s2_figure_data <- prepare_figure_data(s2_annual_path, s2_plot_specs)

save_sensitivity_figure(
  build_figure(s1_figure_data, s1_plot_specs, legend_rows = 1L),
  s1_figure_path
)
save_sensitivity_figure(
  build_figure(s2_figure_data, s2_plot_specs, legend_rows = 2L),
  s2_figure_path
)
