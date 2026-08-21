# ==============================================================================
# 00_functions.r
#
# Shared study definitions and functions.
# ==============================================================================

# ---- Definitional constants ----

# Analysis and reporting windows.
YEARS_MODEL <- 2000:2024
YEARS_REPORT <- 2010:2024

# Age spans shared across the model and results layers.
AGES_MODEL <- 0:85
AGES_MINOR <- 0:17

# Age spans used in rate construction and positivity checks. Keeping them
# together ensures that both operations refer to the same cells.
AGES_FEMALE_REPRO <- 10:54
AGES_MALE_REPRO <- 15:59
AGES_EDUC_TRANSITION <- 0:97

# Paternal age band used for the age-standardised fertility level.
FATHER_AGE_BAND_WIDTH <- 5L

# Sex and education vocabularies (binary BA split).
SEX_LEVELS <- c("female", "male")
EDUC_LEVELS <- c("No_BA", "BA_plus")

# Race/ethnicity vocabularies. R3 = the three headline groups carried through
# the model and results; R5 = the full harmonised set used upstream.
RACES_R3 <- c("white", "black", "hispanic")
RACES_R5 <- c("white", "black", "aian", "api", "hispanic")

# Sensitivity legs shared by the model, table and QC scripts.
MALE_SENSITIVITY_LEGS <- data.frame(
  leg = c(
    "primary",
    "female_borrow",
    "zero_gradient",
    "complete_case",
    "mnar_all_noba",
    "mnar_all_ba"
  ),
  directory = c(
    "s1_primary",
    "s1_female_borrow",
    "s1_zero_gradient",
    "s2_complete_case",
    "s2_mnar_all_noba",
    "s2_mnar_all_ba"
  ),
  stringsAsFactors = FALSE
)
S1_SENSITIVITY_LEGS <- c("zero_gradient", "primary", "female_borrow")
S2_SENSITIVITY_LEGS <- c(
  "primary", "complete_case", "mnar_all_ba", "mnar_all_noba"
)
S2_MNAR_SCENARIOS <- data.frame(
  leg = c("s2_mnar_all_noba", "s2_mnar_all_ba"),
  missing_ba_odds_ratio = c(0, Inf),
  stringsAsFactors = FALSE
)

# ---- Small shared helpers ----

safe_integer <- function(x) suppressWarnings(as.integer(x))

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

# Trim whitespace and map empty strings to NA so blank cells read as missing.
normalize_missing_character <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x
}

# Divide only where the denominator is finite and positive, else NA.
safe_divide <- function(numerator, denominator) {
  fifelse(is.finite(denominator) & denominator > 0, numerator / denominator, NA_real_)
}

# ---- Education recodes (binary BA split: No_BA / BA_plus) ----

# 1989 revision = completed years of schooling: 0-15 -> No_BA;
# 16-17 (college 4+ years) -> BA_plus; anything else -> Unknown.
recode_educ2_from_1989 <- function(educ_1989_raw) {
  educ_num <- safe_integer(educ_1989_raw)
  out <- rep("Unknown", length(educ_num))
  out[!is.na(educ_num) & educ_num >= 0L & educ_num <= 15L] <- "No_BA"
  out[!is.na(educ_num) & educ_num >= 16L & educ_num <= 17L] <- "BA_plus"
  out
}

# 2003 revision = highest degree: 1-5 (below bachelor's, incl. associate) ->
# No_BA; 6-8 (bachelor's / master's / doctorate or professional) -> BA_plus;
# anything else -> Unknown.
recode_educ2_from_2003 <- function(educ_2003_raw) {
  educ_num <- safe_integer(educ_2003_raw)
  out <- rep("Unknown", length(educ_num))
  out[!is.na(educ_num) & educ_num >= 1L & educ_num <= 5L] <- "No_BA"
  out[!is.na(educ_num) & educ_num >= 6L & educ_num <= 8L] <- "BA_plus"
  out
}

# IPUMS-USA EDUCD [detailed]: everything below a bachelor's degree (incl.
# associate) counts as No_BA; bachelor's / master's / professional / doctoral
# count as BA_plus; codes 0/1 (N/A) and 999 (missing) fall through to Unknown
# and are dropped from rate denominators downstream.
recode_educ2_from_educd <- function(educd_raw) {
  educ_num <- safe_integer(educd_raw)
  out <- rep("Unknown", length(educ_num))
  out[!is.na(educ_num) & educ_num %in% c(2L, 10:64, 65L, 70L, 71L, 80L, 81L, 82L, 83L, 90L, 100L, 110:113)] <- "No_BA"
  out[!is.na(educ_num) & educ_num %in% c(101L, 114L, 115L, 116L)] <- "BA_plus"
  out
}

# ---- Education-panel helpers ----

# Children carry no completed education, so their Unknown is reassigned to No_BA
# (symmetrically in numerator and denominator).
bridge_child_unknown_educ2 <- function(dt, age_col = "age_years", educ_col = "educ2",
                                       max_age = 17L, replacement = "No_BA") {
  out <- data.table::copy(dt)
  out[
    !is.na(get(age_col)) & get(age_col) >= 0L & get(age_col) <= max_age &
      get(educ_col) == "Unknown",
    (educ_col) := replacement
  ]
  out[]
}

# Total of one column per education level, as a named vector over EDUC_LEVELS
# with absent levels counted as zero.
sum_col_by_educ <- function(dt, value_col) {
  counts <- dt[
    ,
    .(value = sum(get(value_col), na.rm = TRUE)),
    by = .(educ2)
  ]
  totals <- setNames(counts$value, counts$educ2)[EDUC_LEVELS]
  totals[is.na(totals)] <- 0
  setNames(as.numeric(totals), EDUC_LEVELS)
}

# ---- Cause classification ----

# Collapse the 15-cause classification to the 6-cause grouping.
classify_cause_6_from_cause_15 <- function(cause_15) {
  cause_15 <- as.character(cause_15)
  out <- rep(NA_character_, length(cause_15))

  out[cause_15 == "drug_poisoning_overdose"] <- "drug_poisoning"
  out[cause_15 %in% c("firearm_homicide", "firearm_suicide", "other_firearm_deaths")] <-
    "firearm_deaths"
  out[cause_15 == "covid19"] <- "covid19"
  out[cause_15 %in% c("heart_disease", "cerebrovascular_disease", "diabetes_mellitus")] <-
    "cardiovascular_disease_and_diabetes"
  out[cause_15 == "cancer"] <- "cancer"
  out[!is.na(cause_15) & is.na(out)] <- "all_other_causes"

  out
}

# ---- Figure and table helpers ----
#
# The figure functions below are called with ggplot2, patchwork, scales and ragg
# attached; nothing here loads a package. That matters because every draw stage
# sources this file on the cluster, where the R environment carries only the
# five packages recorded in batch_info (data.table, here, DemoKin, Matrix,
# tidyselect). Function bodies are never evaluated there, so the figure-only
# dependencies cost nothing, but a top-level call into one of them would break
# every draw.

# Presentation constants. Cause keys are spelled out rather than derived from
# CAUSE_6_LEVELS, which is defined further down this file.

FIGURE_FONT <- "Franklin Gothic Book"
YEAR_BREAKS <- c(2010, 2015, 2020, 2024)

# ggsci::pal_jama()(7) with the alpha byte dropped. Hard-coded because the
# cluster R environment does not carry ggsci and every draw sources this file.
JAMA_COLOURS <- c(
  "#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"
)

RACE_LABELS <- c(white = "Non-Hispanic White", black = "Non-Hispanic Black",
                 hispanic = "Hispanic")
EDUC_LABELS <- c(No_BA = "No BA", BA_plus = "BA or higher")
PARENT_LABELS <- c(mother = "Maternal death", father = "Paternal death")
PARENT_WORDS <- c(mother = "maternal", father = "paternal")
CAUSE_6_LABELS <- c(
  drug_poisoning = "Drug poisoning",
  firearm_deaths = "Firearm-related deaths",
  covid19 = "COVID-19",
  cardiovascular_disease_and_diabetes = "Cardiovascular disease and diabetes",
  cancer = "Cancer",
  all_other_causes = "All other causes"
)
GAP_LABELS <- c(RD = "Absolute difference", RR = "Rate ratio")

EDUC_PALETTE <- c(No_BA = JAMA_COLOURS[2], BA_plus = JAMA_COLOURS[3])
RACE_PALETTE <- c(white = JAMA_COLOURS[1], black = JAMA_COLOURS[4],
                  hispanic = JAMA_COLOURS[5])
CAUSE_6_PALETTE <- c(
  covid19 = JAMA_COLOURS[6],
  drug_poisoning = JAMA_COLOURS[4],
  firearm_deaths = JAMA_COLOURS[3],
  cardiovascular_disease_and_diabetes = JAMA_COLOURS[1],
  cancer = JAMA_COLOURS[5],
  all_other_causes = JAMA_COLOURS[7]
)
# COVID-19 sits at the outer edge of the stack; the rest keep the usual order.
CAUSE_6_STACK_LEVELS <- c(
  "covid19", "drug_poisoning", "firearm_deaths",
  "cardiovascular_disease_and_diabetes", "cancer", "all_other_causes"
)

# Axis titles, kept identical across the incidence and mean-age scripts.
AXIS_NUMBER <- "Number of children"
AXIS_RATE <- "Incidence rate, per 1,000 children"
AXIS_RATE_DIFFERENCE <- "Difference in incidence rate, per 1,000 children"
AXIS_MEAN_AGE <- "Mean age, years"
AXIS_MEAN_AGE_DIFFERENCE <- "Difference in mean age, years"

# Shared x-axis scale over the report window.
common_x_scale <- function(pad = 0.5) {
  scale_x_continuous(
    breaks = YEAR_BREAKS,
    limits = c(min(YEARS_REPORT) - pad, max(YEARS_REPORT) + pad),
    expand = expansion(mult = c(0, 0))
  )
}

theme_parent_loss <- function(
  base_size = 10,
  title_size = 11,
  subtitle_size = 9,
  strip_text_size = 9,
  axis_text_size = 8,
  panel_spacing_x = 0.75,
  panel_spacing_y = 0.6,
  legend_key_width_cm = 1.1
) {
  theme_bw(base_size = base_size, base_family = FIGURE_FONT) +
    theme(
      plot.title = element_text(face = "bold", size = title_size,
                                margin = margin(b = 4), hjust = 0),
      plot.subtitle = element_text(size = subtitle_size, colour = "grey30",
                                   margin = margin(b = 8), hjust = 0),
      plot.caption = element_text(size = subtitle_size - 0.5, colour = "grey40",
                                  hjust = 0, margin = margin(t = 8)),
      plot.title.position = "plot",
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.line.x = element_line(colour = "grey15", linewidth = 0.4),
      axis.line.y = element_blank(),
      axis.ticks.x = element_line(colour = "grey15", linewidth = 0.4),
      axis.ticks.y = element_blank(),
      axis.ticks.length = grid::unit(3, "pt"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text = element_text(size = strip_text_size, colour = "grey10"),
      strip.text.y.right = element_text(angle = 270),
      panel.spacing.x = grid::unit(panel_spacing_x, "lines"),
      panel.spacing.y = grid::unit(panel_spacing_y, "lines"),
      legend.position = "bottom",
      legend.title = element_text(size = base_size * 0.9),
      legend.background = element_rect(fill = "white", colour = "grey25",
                                       linewidth = 0.3),
      legend.margin = margin(5, 8, 5, 8),
      legend.key = element_blank(),
      legend.key.width = grid::unit(legend_key_width_cm, "cm"),
      axis.text = element_text(size = axis_text_size, colour = "grey15"),
      axis.title.y = element_text(margin = margin(r = 8)),
      axis.title.x = element_text(margin = margin(t = 6))
    )
}

# Legend inside the plotting area, boxed. The left inset is much smaller than
# the top one because the panel row is far wider than it is tall, and equal
# fractions would not look equal.
inside_legend <- function(position = c("top", "bottom")) {
  position <- match.arg(position)
  theme(
    legend.position = "inside",
    legend.position.inside = if (identical(position, "top")) {
      c(0.015, 0.95)
    } else {
      c(0.015, 0.04)
    },
    legend.justification.inside = if (identical(position, "top")) {
      c(0, 1)
    } else {
      c(0, 0)
    },
    legend.direction = "vertical",
    legend.key.width = grid::unit(0.8, "cm"),
    legend.key.height = grid::unit(0.42, "cm")
  )
}

# Two guides stacked inside one box, for the decomposition panels.
paired_legend_box <- function() {
  theme(
    legend.background = element_blank(),
    legend.box = "vertical",
    legend.box.background = element_rect(fill = "white", colour = "grey25",
                                         linewidth = 0.3),
    legend.box.margin = margin(5, 8, 5, 8),
    legend.box.just = "left",
    legend.spacing.y = grid::unit(7, "pt"),
    legend.key.width = grid::unit(0.58, "cm"),
    legend.key.height = grid::unit(0.42, "cm"),
    legend.key.spacing.y = grid::unit(1.5, "pt"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8.5),
    legend.margin = margin(1, 2, 1, 2)
  )
}

panel_title <- function(tag, description) sprintf("(%s)  %s", tag, description)

# Maternal and paternal results as tagged panels, one above the other.
compose_ab <- function(panel_a, panel_b) {
  (panel_a / panel_b) +
    patchwork::plot_annotation(
      theme = theme_parent_loss() + theme(plot.margin = margin(6, 6, 4, 6))
    )
}

# Six cause panels in a 3x2 grid under one shared y-axis title.
compose_causes <- function(panels, y_label) {
  y_title <- patchwork::wrap_elements(full = grid::textGrob(
    y_label, rot = 90,
    gp = grid::gpar(fontfamily = FIGURE_FONT, fontsize = 10, col = "black")
  ))
  (y_title | patchwork::wrap_plots(panels, ncol = 2)) +
    patchwork::plot_layout(widths = c(1, 46)) +
    patchwork::plot_annotation(
      theme = theme_parent_loss() + theme(plot.margin = margin(6, 6, 4, 2))
    )
}

make_axis_step <- function(raw_step) {
  magnitude <- 10^floor(log10(raw_step))
  multipliers <- c(1, 2, 2.5, 3, 4, 5, 8, 10)
  multipliers[which(multipliers * magnitude >= raw_step)[1L]] * magnitude
}

make_zero_based_axis_breaks <- function(max_value, target_intervals = 4L) {
  if (!is.finite(max_value) || max_value <= 0) return(c(0, 1))
  step <- make_axis_step(max_value / target_intervals)
  seq(0, ceiling(max_value / step) * step, by = step)
}

# Axis for a stack that may run either side of zero. The step comes from the
# whole span, and the limits are the slackened values themselves rather than
# rounded outwards again, which would buy a whole empty step.
make_stacked_axis <- function(min_value, max_value, target_intervals = 5L,
                              slack = 1.2) {
  lower <- min(min_value * slack, 0)
  upper <- max(max_value * slack, 0)
  span <- max(abs(lower), abs(upper), 1e-9)
  step <- make_axis_step(span / target_intervals)
  breaks <- seq(floor(lower / step) * step, ceiling(upper / step) * step, by = step)
  list(breaks = breaks[breaks >= lower & breaks <= upper], limits = c(lower, upper))
}

save_figure <- function(plot_obj, png_path, width = 10.5, height = 8.4, dpi = 320) {
  normalized <- normalizePath(png_path, mustWork = FALSE)
  write_path <- if (.Platform$OS.type == "windows" && nchar(normalized) >= 255L) {
    tempfile(fileext = ".png")
  } else {
    png_path
  }
  ggsave(filename = write_path, plot = plot_obj, device = ragg::agg_png,
         width = width, height = height, dpi = dpi, bg = "white")
  if (!identical(write_path, png_path)) {
    if (!isTRUE(file.copy(write_path, png_path, overwrite = TRUE))) {
      stop(sprintf("Failed to copy temporary figure to target path: %s", png_path), call. = FALSE)
    }
    unlink(write_path)
  }
  invisible(png_path)
}

# ---- Parent-loss endpoints ----

# 3-cause all-cause partition, cause_6 grouping, parent types, and the education
# levels the endpoint tables carry (the two BA legs plus a pooled Total).
CAUSE_LEVELS <- c("drugs", "firearms", "other")
CAUSE_6_LEVELS <- c(
  "drug_poisoning",
  "firearm_deaths",
  "covid19",
  "cardiovascular_disease_and_diabetes",
  "cancer",
  "all_other_causes"
)
PARENT_TYPE_LEVELS <- c("mother", "father")
EDUC_TABLE_LEVELS <- c(EDUC_LEVELS, "Total")

# All-cause incidence: sum the 3-cause partition; the denominator is the
# cell-level alive-parent exposure (constant across cause), so first() recovers it.
build_total_incidence <- function(annual_dt, educ_levels) {
  dt <- as.data.table(annual_dt)[
    year %in% YEARS_REPORT &
      race_eth %chin% RACES_R3 &
      parent_educ2 %chin% educ_levels &
      cause %chin% CAUSE_LEVELS &
      parent_type %chin% PARENT_TYPE_LEVELS
  ]
  out <- dt[
    , .(
      n_children_losing_parent = sum(n_children_losing_parent, na.rm = TRUE),
      n_children_with_alive_parent_under18 = first(n_children_with_alive_parent_under18)
    ),
    by = .(year, race_eth, parent_educ2, parent_type)
  ]
  out[, parental_loss_per_1000 := 1000 * n_children_losing_parent / n_children_with_alive_parent_under18]
  out[]
}

# Cause_6 incidence: roll the cause_15 panel up to the 6 categories (shared cell denominator).
build_cause6_incidence <- function(cause_6_15_dt, educ_levels) {
  dt <- as.data.table(cause_6_15_dt)[
    year %in% YEARS_REPORT &
      race_eth %chin% RACES_R3 &
      parent_educ2 %chin% educ_levels &
      cause_6 %chin% CAUSE_6_LEVELS &
      parent_type %chin% PARENT_TYPE_LEVELS
  ]
  out <- dt[
    , .(
      n_children_losing_parent = sum(n_children_losing_parent, na.rm = TRUE),
      n_children_with_alive_parent_under18 = first(n_children_with_alive_parent_under18)
    ),
    by = .(year, race_eth, parent_educ2, cause_6, parent_type)
  ]
  out[, parental_loss_per_1000 := 1000 * n_children_losing_parent / n_children_with_alive_parent_under18]
  out[]
}

# All-cause event-weighted mean child age at loss: pool the 3-cause partition to
# child-age loss counts, then take the event-weighted mean child age over 0-17.
build_total_mean_age <- function(child_age_dt, educ_levels) {
  dt <- as.data.table(child_age_dt)[
    year %in% YEARS_REPORT &
      race_eth %chin% RACES_R3 &
      parent_educ2 %chin% educ_levels &
      cause %chin% CAUSE_LEVELS &
      parent_type %chin% PARENT_TYPE_LEVELS &
      age_child %in% AGES_MINOR
  ]
  by_age <- dt[
    , .(loss = sum(as.numeric(n_children_losing_parent_at_child_age), na.rm = TRUE)),
    by = .(year, race_eth, parent_educ2, parent_type, age_child)
  ]
  by_age[
    , {
      total_loss <- sum(loss, na.rm = TRUE)
      .(
        n_children_losing_parent = total_loss,
        mean_age_at_loss = safe_divide(sum(age_child * loss, na.rm = TRUE), total_loss)
      )
    },
    by = .(year, race_eth, parent_educ2, parent_type)
  ]
}

# Cause_6 mean age: same weighting within each cause_6; an undefined mean where a
# cell has no events in that cause-year (e.g. COVID-19 before 2020) is left NA.
build_cause6_mean_age <- function(child_age_dt, educ_levels) {
  dt <- as.data.table(child_age_dt)[
    year %in% YEARS_REPORT &
      race_eth %chin% RACES_R3 &
      parent_educ2 %chin% educ_levels &
      cause_6 %chin% CAUSE_6_LEVELS &
      parent_type %chin% PARENT_TYPE_LEVELS &
      age_child %in% AGES_MINOR
  ]
  by_age <- dt[
    , .(loss = sum(as.numeric(n_children_losing_parent_at_child_age), na.rm = TRUE)),
    by = .(year, race_eth, parent_educ2, cause_6, parent_type, age_child)
  ]
  by_age[
    , {
      total_loss <- sum(loss, na.rm = TRUE)
      .(
        n_children_losing_parent = total_loss,
        mean_age_at_loss = safe_divide(sum(age_child * loss, na.rm = TRUE), total_loss)
      )
    },
    by = .(year, race_eth, parent_educ2, cause_6, parent_type)
  ]
}

# Education gap: No-BA minus / over BA-plus in the all-cause incidence rate.
build_gap_data <- function(total_incidence_dt) {
  wide <- dcast(
    total_incidence_dt[parent_educ2 %chin% c("No_BA", "BA_plus")],
    year + race_eth + parent_type ~ parent_educ2,
    value.var = "parental_loss_per_1000"
  )
  gap <- rbindlist(list(
    wide[, .(year, race_eth, parent_type, gap_metric = "RD", value = No_BA - BA_plus)],
    wide[, .(year, race_eth, parent_type, gap_metric = "RR", value = safe_divide(No_BA, BA_plus))]
  ))
  gap[, `:=`(race_order = match(race_eth, RACES_R3), parent_order = match(parent_type, PARENT_TYPE_LEVELS))]
  setorder(gap, parent_order, race_order, gap_metric, year)
  gap[, c("race_order", "parent_order") := NULL]
  gap[]
}

# The gap meanage(No_BA) - meanage(BA_plus) is a difference of event-weighted
# means, not an additive rate difference. The frame below carries each leg's
# cause shares and cause-specific mean ages.
build_meanage_gap_frame <- function(cause6_annual_dt) {
  dt <- copy(cause6_annual_dt)[
    parent_educ2 %chin% EDUC_LEVELS & cause_6 %chin% CAUSE_6_LEVELS
  ]
  dt[, leg_total := sum(n_children_losing_parent, na.rm = TRUE),
     by = .(year, race_eth, parent_type, parent_educ2)]
  dt[, share := fifelse(leg_total > 0, n_children_losing_parent / leg_total, 0)]
  dt[, leg_total := NULL]

  w <- dcast(
    dt,
    year + race_eth + parent_type + cause_6 ~ parent_educ2,
    value.var = c("share", "mean_age_at_loss")
  )
  # Backfill the unused mean age for causes with zero weight in one education leg.
  for (col in c("share_No_BA", "share_BA_plus")) {
    w[is.na(get(col)), (col) := 0]
  }
  w[, a_N := fifelse(is.finite(mean_age_at_loss_No_BA), mean_age_at_loss_No_BA, mean_age_at_loss_BA_plus)]
  w[, a_B := fifelse(is.finite(mean_age_at_loss_BA_plus), mean_age_at_loss_BA_plus, mean_age_at_loss_No_BA)]
  w[is.na(a_N), a_N := 0]
  w[is.na(a_B), a_B := 0]

  w[, `:=`(
    meanage_No_BA = sum(share_No_BA * a_N),
    meanage_BA_plus = sum(share_BA_plus * a_B)
  ), by = .(year, race_eth, parent_type)]
  w[, mean_age_gap_total := meanage_No_BA - meanage_BA_plus]

  out <- w[, .(
    year, race_eth, parent_type, cause_6,
    share_No_BA, share_BA_plus,
    mean_age_No_BA = a_N, mean_age_BA_plus = a_B,
    mean_age_gap_total
  )]
  out[, `:=`(
    race_order = match(race_eth, RACES_R3),
    parent_order = match(parent_type, PARENT_TYPE_LEVELS),
    cause_order = match(cause_6, CAUSE_6_LEVELS)
  )]
  setorder(out, parent_order, race_order, year, cause_order)
  out[, c("race_order", "parent_order", "cause_order") := NULL]
  out[]
}

# Reference-free cause decomposition of the mean-age education gap. For
# coalition S, causes in S take their No_BA (share, age) and the rest take
# BA_plus; each cause's value averages its marginal effect over all orderings,
# and the values sum to the total gap.
build_meanage_gap_shapley <- function(cause6_annual_dt) {
  n <- length(CAUSE_6_LEVELS)
  masks <- vapply(
    seq_len(n),
    function(j) rep(rep(c(FALSE, TRUE), each = 2^(j - 1)), length.out = 2^n),
    logical(2^n)
  ) # 2^n x n; column j toggles every 2^(j-1) rows (0-based index = sum bit_j 2^(j-1))
  sizes <- rowSums(masks)
  # The grand coalition has no absent-cause weight.
  size_weight <- numeric(length(sizes))
  valid_coalition <- sizes < n
  size_weight[valid_coalition] <- factorial(sizes[valid_coalition]) *
    factorial(n - sizes[valid_coalition] - 1L) / factorial(n)
  # For each cause c: rows where c is absent, and the matching row with c added.
  absent_rows <- lapply(seq_len(n), function(c) which(!masks[, c]))
  with_rows <- lapply(seq_len(n), function(c) absent_rows[[c]] + 2L^(c - 1L))

  shapley_one <- function(sN, aN, sB, aB) {
    sna <- sN * aN
    sba <- sB * aB
    num <- as.numeric(masks %*% sna + (!masks) %*% sba)
    den <- as.numeric(masks %*% sN + (!masks) %*% sB)
    M <- ifelse(den > 0, num / den, NA_real_)
    vapply(seq_len(n), function(c) {
      S <- absent_rows[[c]]
      sum(size_weight[S] * (M[with_rows[[c]]] - M[S]))
    }, numeric(1))
  }

  dt <- cause6_annual_dt[
    parent_educ2 %chin% EDUC_LEVELS & cause_6 %chin% CAUSE_6_LEVELS
  ]
  cells <- unique(dt[, .(year, race_eth, parent_type)])
  res <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    cy <- cells$year[i]; cr <- cells$race_eth[i]; cp <- cells$parent_type[i]
    w <- dcast(
      dt[year == cy & race_eth == cr & parent_type == cp],
      cause_6 ~ parent_educ2,
      value.var = c("n_children_losing_parent", "mean_age_at_loss")
    )
    setkey(w, cause_6)
    w <- w[CAUSE_6_LEVELS]
    nN <- w$n_children_losing_parent_No_BA; nN[is.na(nN)] <- 0
    nB <- w$n_children_losing_parent_BA_plus; nB[is.na(nB)] <- 0
    sN <- if (sum(nN) > 0) nN / sum(nN) else nN
    sB <- if (sum(nB) > 0) nB / sum(nB) else nB
    aN <- w$mean_age_at_loss_No_BA; aN[is.na(aN)] <- 0
    aB <- w$mean_age_at_loss_BA_plus; aB[is.na(aB)] <- 0
    res[[i]] <- data.table(
      year = cy, race_eth = cr, parent_type = cp,
      cause_6 = CAUSE_6_LEVELS,
      shapley_contribution = shapley_one(sN, aN, sB, aB)
    )
  }
  rbindlist(res)
}

build_meanage_decomposition <- function(cause6_annual_dt) {
  out <- merge(
    build_meanage_gap_frame(cause6_annual_dt),
    build_meanage_gap_shapley(cause6_annual_dt),
    by = c("year", "race_eth", "parent_type", "cause_6"), all.x = TRUE, sort = FALSE
  )
  setcolorder(out, c(
    "year", "race_eth", "parent_type", "cause_6",
    "share_No_BA", "share_BA_plus", "mean_age_No_BA", "mean_age_BA_plus",
    "shapley_contribution", "mean_age_gap_total"
  ))
  out[]
}

# Reference row order for the total and cause_6 endpoint tables.
order_total <- function(dt) {
  dt <- copy(dt)
  dt[, `:=`(
    race_order = match(race_eth, RACES_R3),
    educ_order = match(parent_educ2, EDUC_TABLE_LEVELS),
    parent_order = match(parent_type, PARENT_TYPE_LEVELS)
  )]
  setorder(dt, race_order, educ_order, parent_order, year)
  dt[, c("race_order", "educ_order", "parent_order") := NULL]
  dt[]
}

order_cause6 <- function(dt) {
  dt <- copy(dt)
  dt[, `:=`(
    race_order = match(race_eth, RACES_R3),
    educ_order = match(parent_educ2, EDUC_TABLE_LEVELS),
    cause_order = match(cause_6, CAUSE_6_LEVELS),
    parent_order = match(parent_type, PARENT_TYPE_LEVELS)
  )]
  setorder(dt, race_order, educ_order, cause_order, parent_order, year)
  dt[, c("race_order", "educ_order", "cause_order", "parent_order") := NULL]
  dt[]
}

# ---- Uncertainty shards ----

# Check the four endpoint tables and their closure identities.
check_shard_tables <- function(tables) {
  errs <- character(0)
  add <- function(msg) errs <<- c(errs, msg)

  n_cell <- length(YEARS_REPORT) * length(RACES_R3) *
    length(EDUC_TABLE_LEVELS) * length(PARENT_TYPE_LEVELS)
  n_cause <- n_cell * length(CAUSE_6_LEVELS)
  key_total <- c("year", "race_eth", "parent_educ2", "parent_type")
  key_cause <- c(key_total, "cause_6")

  it <- tables$incidence_total
  ic <- tables$incidence_cause6
  mt <- tables$meanage_total
  mc <- tables$meanage_cause6

  # Row grid and key uniqueness.
  if (nrow(it) != n_cell) add(sprintf("incidence_total has %d rows, expected %d.", nrow(it), n_cell))
  if (nrow(mt) != n_cell) add(sprintf("meanage_total has %d rows, expected %d.", nrow(mt), n_cell))
  if (nrow(ic) != n_cause) add(sprintf("incidence_cause6 has %d rows, expected %d.", nrow(ic), n_cause))
  if (nrow(mc) != n_cause) add(sprintf("meanage_cause6 has %d rows, expected %d.", nrow(mc), n_cause))
  if (anyDuplicated(it, by = key_total)) add("incidence_total keys not unique.")
  if (anyDuplicated(mt, by = key_total)) add("meanage_total keys not unique.")
  if (anyDuplicated(ic, by = key_cause)) add("incidence_cause6 keys not unique.")
  if (anyDuplicated(mc, by = key_cause)) add("meanage_cause6 keys not unique.")

  # Incidence numeric domain: finite, non-negative, denominator strictly positive.
  for (nm in c("incidence_total", "incidence_cause6")) {
    d <- tables[[nm]]
    if (!all(is.finite(d$n_children_losing_parent) & d$n_children_losing_parent >= 0)) {
      add(sprintf("%s: n_children_losing_parent not finite and non-negative.", nm))
    }
    if (!all(is.finite(d$n_children_with_alive_parent_under18) & d$n_children_with_alive_parent_under18 > 0)) {
      add(sprintf("%s: denominator not finite and positive.", nm))
    }
    if (!all(is.finite(d$parental_loss_per_1000) & d$parental_loss_per_1000 >= 0)) {
      add(sprintf("%s: parental_loss_per_1000 not finite and non-negative.", nm))
    }
  }

  # Mean age: counts finite and non-negative; mean age NA only where the count is
  # 0, otherwise finite and within [0, oldest minor age].
  for (nm in c("meanage_total", "meanage_cause6")) {
    d <- tables[[nm]]
    if (!all(is.finite(d$n_children_losing_parent) & d$n_children_losing_parent >= 0)) {
      add(sprintf("%s: n_children_losing_parent not finite and non-negative.", nm))
    }
    bad_na <- d[is.na(mean_age_at_loss) & n_children_losing_parent > 0]
    if (nrow(bad_na)) add(sprintf("%s: mean_age_at_loss is NA where events > 0 (%d cells).", nm, nrow(bad_na)))
    fin <- d[is.finite(mean_age_at_loss)]
    if (nrow(fin) && (min(fin$mean_age_at_loss) < 0 || max(fin$mean_age_at_loss) > max(AGES_MINOR))) {
      add(sprintf("%s: finite mean_age_at_loss outside [0, %d].", nm, max(AGES_MINOR)))
    }
  }

  # Closure (1): cause_6 incidence denominator constant across cause within a cell.
  denom_var <- ic[
    , .(u = uniqueN(round(n_children_with_alive_parent_under18, 6))), by = key_total
  ][u != 1L]
  if (nrow(denom_var)) add("cause_6 incidence denominator varies across cause within a cell.")

  # Closure (2): the six cause loss counts sum to the all-cause loss.
  chk2 <- merge(
    ic[, .(loss6 = sum(n_children_losing_parent)), by = key_total],
    it[, c(key_total, "n_children_losing_parent"), with = FALSE],
    by = key_total
  )
  if (nrow(chk2) != n_cell || chk2[, max(abs(loss6 - n_children_losing_parent))] > 1e-6) {
    add("cause_6 loss counts do not sum to the all-cause loss.")
  }

  # Closure (3): incidence and mean-age loss counts agree per cell (two panels).
  chk3t <- merge(
    it[, c(key_total, "n_children_losing_parent"), with = FALSE],
    mt[, c(key_total, "n_children_losing_parent"), with = FALSE],
    by = key_total, suffixes = c("_inc", "_age")
  )
  if (chk3t[, max(abs(n_children_losing_parent_inc - n_children_losing_parent_age))] > 1e-6) {
    add("total incidence and mean-age loss counts disagree.")
  }
  chk3c <- merge(
    ic[, c(key_cause, "n_children_losing_parent"), with = FALSE],
    mc[, c(key_cause, "n_children_losing_parent"), with = FALSE],
    by = key_cause, suffixes = c("_inc", "_age")
  )
  if (chk3c[, max(abs(n_children_losing_parent_inc - n_children_losing_parent_age))] > 1e-6) {
    add("cause_6 incidence and mean-age loss counts disagree.")
  }

  # Closure (4): cause_6 mean ages weighted by loss close to the all-cause mean age.
  wt <- mc[
    n_children_losing_parent > 0 & is.finite(mean_age_at_loss),
    .(w_age = sum(n_children_losing_parent * mean_age_at_loss) / sum(n_children_losing_parent)),
    by = key_total
  ]
  chk4 <- merge(wt, mt[n_children_losing_parent > 0, c(key_total, "mean_age_at_loss"), with = FALSE], by = key_total)
  if (nrow(chk4) && chk4[, max(abs(w_age - mean_age_at_loss))] > 1e-6) {
    add("cause_6 mean ages do not close to the all-cause mean age.")
  }

  errs
}

uncertainty_field <- function(object, name, default = NA) {
  if (!is.list(object)) return(default)
  value <- object[[name]]
  if (is.null(value) || !length(value)) default else value[[1]]
}

# ---- ACS harmonisation ----

# Locate the single IPUMS-USA DDI/XML companion file and its matching data file.
# The caller refuses to guess fixed-width positions, so missing metadata is a
# hard error.
find_acs_files <- function(raw_dir) {
  if (!dir.exists(raw_dir)) {
    stop("Raw ACS directory does not exist: ", raw_dir, call. = FALSE)
  }

  xml_files <- list.files(raw_dir, pattern = "^usa_.*\\.xml$", full.names = TRUE, ignore.case = TRUE)
  data_files <- list.files(raw_dir, pattern = "\\.(dat|dat\\.gz)$", full.names = TRUE, ignore.case = TRUE)

  if (!length(xml_files)) {
    stop(
      "Missing IPUMS-USA companion metadata in ", raw_dir,
      ". Add the matching DDI/XML file (for example, usa_00001.xml) before running.",
      call. = FALSE
    )
  }
  if (length(xml_files) > 1L) {
    stop("Multiple DDI/XML files found in ", raw_dir, ". Keep only one companion metadata file.", call. = FALSE)
  }
  if (!length(data_files)) {
    stop("No IPUMS data file (.dat/.dat.gz) found in ", raw_dir, call. = FALSE)
  }

  ddi_file <- xml_files[[1]]
  ddi_base <- sub("\\.xml$", "", basename(ddi_file), ignore.case = TRUE)

  matching_data <- data_files[
    basename(data_files) %in% c(paste0(ddi_base, ".dat"), paste0(ddi_base, ".dat.gz"))
  ]
  if (!length(matching_data)) matching_data <- data_files

  # Prefer the uncompressed .dat when both forms are present.
  preferred_data <- matching_data[grepl("\\.dat$", matching_data, ignore.case = TRUE)]
  if (!length(preferred_data)) preferred_data <- matching_data
  if (length(preferred_data) > 1L) {
    stop("Multiple IPUMS data files found in ", raw_dir, ". Keep only one .dat/.dat.gz file.", call. = FALSE)
  }

  list(ddi_file = ddi_file, data_file = preferred_data[[1]])
}

# IPUMS-USA SEX: 1 = male, 2 = female (no missing code on person records).
harmonize_acs_sex <- function(sex_raw) {
  sex_num <- safe_integer(sex_raw)
  out <- rep(NA_character_, length(sex_num))
  out[!is.na(sex_num) & sex_num == 1L] <- "male"
  out[!is.na(sex_num) & sex_num == 2L] <- "female"
  out
}

# IPUMS AGE is single-year and already top-coded by IPUMS; treat clearly
# out-of-range values as missing.
harmonize_acs_age <- function(age_raw) {
  age_num <- safe_integer(age_raw)
  age_num[!is.na(age_num) & (age_num < 0L | age_num > 120L)] <- NA_integer_
  age_num
}

# IPUMS-USA HISPAN [general version]: 0 = not Hispanic, 1-4 = Hispanic-origin
# groups (Mexican / Puerto Rican / Cuban / other), 9 = not reported. Unreported
# origin stays NA, matching the single-race death and birth numerators.
acs_hispanic_flag <- function(hispanic_raw) {
  hisp_num <- safe_integer(hispanic_raw)
  out <- rep(NA, length(hisp_num))
  out[!is.na(hisp_num) & hisp_num == 0L] <- FALSE
  out[!is.na(hisp_num) & hisp_num >= 1L & hisp_num <= 4L] <- TRUE
  out
}

# IPUMS-USA RACE [general version] single-race codes: 1 white, 2 black,
# 3 American Indian/Alaska Native, 4-6 Asian or Pacific Islander
# (Chinese / Japanese / other Asian or PI). Codes 7-9 (other race n.e.c., two or
# three or more major races) have no single-race equivalent and become NA, so
# multiracial respondents are dropped symmetrically with the single-race
# death/birth numerators (Asian + Pacific Islander merged per Villaveces et al., 2025).
acs_race_group <- function(race_raw) {
  race_num <- safe_integer(race_raw)
  out <- rep(NA_character_, length(race_num))
  out[!is.na(race_num) & race_num == 1L] <- "white"
  out[!is.na(race_num) & race_num == 2L] <- "black"
  out[!is.na(race_num) & race_num == 3L] <- "aian"
  out[!is.na(race_num) & race_num >= 4L & race_num <= 6L] <- "api"
  out
}

harmonize_acs_race_eth <- function(race_raw, hispanic_raw) {
  race_group <- acs_race_group(race_raw)
  hispanic_flag <- acs_hispanic_flag(hispanic_raw)

  out <- race_group
  out[hispanic_flag %in% TRUE] <- "hispanic"
  out
}

# Mortality denominator: cap age at 85 to match the 85+ top-code in the death
# numerator, then count by analysis cell.
derive_death_denominator <- function(panel_df) {
  if (!nrow(panel_df)) {
    return(data.frame(
      year = integer(0),
      sex = character(0),
      age_years = integer(0),
      race_eth = character(0),
      educ2 = character(0),
      n_pop_unweighted = integer(0),
      n_pop_weighted = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  dt <- data.table::as.data.table(panel_df)
  dt[!is.na(age_years) & age_years > max(AGES_MODEL), age_years := max(AGES_MODEL)]
  out <- dt[, .(
    n_pop_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
    n_pop_weighted = sum(n_persons_weighted, na.rm = TRUE)
  ), by = .(year, sex, age_years, race_eth, educ2)]
  data.table::setorderv(out, c("year", "sex", "age_years", "race_eth", "educ2"))
  as.data.frame(out)
}

# Fertility denominator: women only, counted by analysis cell. Age is not capped
# because the fertility rate layer works within the childbearing age window.
derive_birth_denominator <- function(panel_df) {
  if (!nrow(panel_df)) {
    return(data.frame(
      year = integer(0),
      age_years = integer(0),
      race_eth = character(0),
      educ2 = character(0),
      n_women_unweighted = integer(0),
      n_women_weighted = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  dt <- data.table::as.data.table(panel_df)
  dt <- dt[sex %in% "female"]
  out <- dt[, .(
    n_women_unweighted = sum(n_persons_unweighted, na.rm = TRUE),
    n_women_weighted = sum(n_persons_weighted, na.rm = TRUE)
  ), by = .(year, age_years, race_eth, educ2)]
  data.table::setorderv(out, c("year", "age_years", "race_eth", "educ2"))
  as.data.frame(out)
}

# ---- ACS replicate weights ----

# The ACS variance design is successive-difference replication: 80 replicate
# weights and the 4/80 scale (Census Bureau, ACS Design and Methodology ch. 12).
# The uncertainty draw perturbs the population with a discrete Rademacher
# multiplier on the replicate deviations. Replicate weights are published from
# 2005 on, so earlier ACS years carry the
# point weight unperturbed and contribute no population uncertainty.
ACS_REPLICATE_COUNT <- 80L
ACS_SDR_SCALE <- 4 / 80
ACS_REPLICATE_MIN_YEAR <- 2005L

# A replicate count of zero would send the log below to -Inf. Half a weighted
# person is far under any real cell and keeps the deviation finite.
ACS_REPLICATE_ZERO_FLOOR <- 0.5

# Column names carrying the replicate-weighted counts.
acs_replicate_columns <- function(prefix = "n_persons_weighted_rep") {
  sprintf("%s%02d", prefix, seq_len(ACS_REPLICATE_COUNT))
}

# A population is a positive count, so the multiplier acts on the log scale and
# keeps the draw inside the support; the linear form leaves it in cells the
# replicate design cannot resolve. To first order the two agree, so the interval
# changes only where the linear form was producing negative populations.
acs_log_deviations <- function(point, replicate_matrix) {
  log(pmax(replicate_matrix, ACS_REPLICATE_ZERO_FLOOR)) -
    log(pmax(point, ACS_REPLICATE_ZERO_FLOOR))
}

# ---- Uncertainty run definitions ----

# L'Ecuyer-CMRG seed used to construct the draw ledger.
DRAW_MASTER_SEED <- 20260725L
DRAW_RNG_KIND <- "L'Ecuyer-CMRG"

# Split one draw stream into fixed components so that changes in one
# perturbation do not alter the others.
draw_component_streams <- function(draw_stream) {
  death <- parallel::nextRNGSubStream(draw_stream)
  birth <- parallel::nextRNGSubStream(death)
  acs <- parallel::nextRNGSubStream(birth)
  list(death = death, birth = birth, acs = acs)
}

current_commit <- function() {
  from_environment <- Sys.getenv(
    "EDU_PARENT_LOSS_CODE_VERSION", unset = ""
  )
  if (nzchar(from_environment)) return(from_environment)

  from_git <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)[1],
    error = function(e) NA_character_
  )
  if (length(from_git) && !is.na(from_git) && nzchar(from_git)) return(from_git)
  NA_character_
}
