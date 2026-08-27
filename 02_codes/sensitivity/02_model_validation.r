# ==============================================================================
# 02_model_validation.r
#
# Compare parental-death prevalence estimates with Schlüter et al. (2024) and
# the National Survey of Children's Health (NSCH), 2016-2020, and generate
# Appendix Figure S8.
# ==============================================================================

# ---- Packages ----

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(ragg)
})

here::i_am("02_codes/sensitivity/02_model_validation.r")

source(here::here("02_codes", "00_functions.r"))

# ---- Paths ----

exposure_path <- here(
  "01_data", "02_processed", "edu_parent_loss",
  "01_education_parent_alive_exposure_pooled.rds"
)
population_path <- here(
  "01_data", "02_processed", "population", "01_pop_main_panel.rds"
)
validation_path <- here(
  "01_data", "01_raw", "external",
  "nsch_validation_race.csv"
)

figure_dir <- here(
  "03_outputs", "03_appendix", "figures"
)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

figure_png_path <- file.path(
  figure_dir, "Fig_S8_model_result_validation.png"
)

# ---- Inputs ----

validation_years <- 2016:2020
exposure <- as.data.table(readRDS(exposure_path))
population <- as.data.table(readRDS(population_path))
published <- fread(validation_path)

# ---- Study estimates ----

parent_survival <- exposure[
  year %in% validation_years &
    race_eth %chin% RACES_R3 &
    age_child %in% AGES_MINOR,
  .(survival = sum(n_parents_alive)),
  by = .(year, race_eth, age_child, sex_parent)
]
parent_survival <- dcast(
  parent_survival,
  year + race_eth + age_child ~ sex_parent,
  value.var = "survival"
)

child_population <- population[
  source == "acs" &
    year %in% validation_years &
    race_eth %chin% RACES_R3 &
    age_years %in% AGES_MINOR,
  .(n_children = sum(n_persons_weighted)),
  by = .(year, race_eth, age_child = age_years)
]

model_by_age <- merge(
  parent_survival,
  child_population,
  by = c("year", "race_eth", "age_child"),
  all = FALSE,
  sort = FALSE
)
model_by_age[, prevalence := 1 - female * male]

ours_by_race <- model_by_age[
  , .(estimate_pct = 100 * weighted.mean(prevalence, n_children)),
  by = .(year, race_eth)
]
ours_total <- model_by_age[
  , .(estimate_pct = 100 * weighted.mean(prevalence, n_children)),
  by = year
]
ours_total[, race_eth := "total"]

ours <- rbindlist(list(ours_by_race, ours_total), use.names = TRUE)
ours[, `:=`(
  source = "Ours",
  lower95_pct = NA_real_,
  upper95_pct = NA_real_
)]

# ---- Published estimates ----

race_lookup <- c(
  "Non-Hispanic White" = "white",
  "Non-Hispanic Black" = "black",
  "Hispanic" = "hispanic"
)
source_lookup <- c(
  "Matrix Projection Model" = "Schlüter et al., 2024",
  "National Survey of Children's Health" = "NSCH"
)

published_r3 <- published[
  race_cat %chin% names(race_lookup),
  .(
    source = unname(source_lookup[source]),
    race_eth = unname(race_lookup[race_cat]),
    year = as.integer(year),
    estimate_pct = as.numeric(median),
    lower95_pct = as.numeric(lower95),
    upper95_pct = as.numeric(upper95)
  )
]
race_weights <- child_population[
  , .(n_children = sum(n_children)),
  by = .(year, race_eth)
]
race_weights[, weight := n_children / sum(n_children), by = year]

published_weighted <- merge(
  published_r3,
  race_weights[, .(year, race_eth, weight)],
  by = c("year", "race_eth"),
  all = FALSE,
  sort = FALSE
)
published_weighted[, se_pct := (upper95_pct - lower95_pct) / (2 * 1.96)]

published_total <- published_weighted[
  , .(
    estimate_pct = sum(weight * estimate_pct),
    se_pct = sqrt(sum((weight * se_pct)^2))
  ),
  by = .(source, year)
]
published_total[, `:=`(
  race_eth = "total",
  lower95_pct = pmax(0, estimate_pct - 1.96 * se_pct),
  upper95_pct = estimate_pct + 1.96 * se_pct
)]
published_total[, se_pct := NULL]

validation_results <- rbindlist(
  list(ours, published_r3, published_total),
  use.names = TRUE
)
validation_results[, race_eth := factor(
  race_eth,
  levels = c("total", "white", "black", "hispanic")
)]
setorder(validation_results, race_eth, source, year)

figure_data_dir <- here::here("03_outputs", "04_figure_data")
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(
  validation_results,
  file.path(figure_data_dir, "Fig_S8_model_result_validation_data.csv")
)

# ---- Figure ----

panel_labels <- c(
  total = "Total",
  white = "Non-Hispanic White",
  black = "Non-Hispanic Black",
  hispanic = "Hispanic"
)
source_levels <- c("Ours", "Schlüter et al., 2024", "NSCH")
source_colours <- c(
  "Ours" = "#009CC6",
  "Schlüter et al., 2024" = "#E38B3A",
  "NSCH" = "#5F6368"
)

plot_data <- copy(validation_results)
plot_data[, panel := factor(
  panel_labels[as.character(race_eth)],
  levels = unname(panel_labels)
)]

p <- ggplot(
  plot_data,
  aes(x = year, y = estimate_pct, colour = source, group = source)
) +
  geom_ribbon(
    data = plot_data[source != "Ours"],
    aes(ymin = lower95_pct, ymax = upper95_pct, fill = source),
    alpha = 0.18,
    colour = NA
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = source), size = 2) +
  facet_wrap(~panel, ncol = 2) +
  scale_colour_manual(
    values = source_colours,
    breaks = source_levels,
    name = NULL
  ) +
  scale_fill_manual(
    values = source_colours,
    breaks = source_levels,
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Ours" = 16,
      "Schlüter et al., 2024" = 15,
      "NSCH" = 16
    ),
    breaks = source_levels,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(2016, 2018, 2020),
    limits = range(validation_years),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(0, 12, 2.5),
    limits = c(0, 12),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Year",
    y = "Children experienced at least 1 parental death (%)"
  ) +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 10,
      colour = "#333333",
      margin = margin(b = 7)
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#D9D9D9", linewidth = 0.4),
    axis.line.x = element_line(colour = "#333333", linewidth = 0.55),
    axis.ticks.x = element_line(colour = "#333333", linewidth = 0.4),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 10, colour = "#222222"),
    axis.text = element_text(size = 9, colour = "#555555"),
    panel.spacing.x = grid::unit(26, "pt"),
    panel.spacing.y = grid::unit(16, "pt"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 9, colour = "#333333"),
    legend.spacing.x = grid::unit(8, "pt"),
    plot.margin = margin(8, 10, 6, 8)
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    shape = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    fill = "none"
  )

ggsave(
  figure_png_path,
  p,
  width = 9,
  height = 7,
  dpi = 300,
  device = ragg::agg_png,
  units = "in"
)
