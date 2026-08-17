#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
OUT_DIR <- file.path(PROJECT_DIR, "outputs", "manuscript_story_analysis", "aerobic_only_physiology_figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

theme_phys <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
      strip.background = element_rect(fill = "grey94", color = "grey75"),
      strip.text = element_text(face = "bold", size = 10),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 11),
      legend.title = element_text(face = "bold")
    )
}

summary_se <- function(df, value_col) {
  df %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = ifelse(n > 1, sd(.data[[value_col]], na.rm = TRUE), 0),
      se = ifelse(n > 1, sd / sqrt(n), 0),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
  )
}

ppm_per_mM <- 12

day_from_timepoint <- function(timepoint) {
  case_when(
    timepoint == "T0" ~ 0,
    timepoint == "T1" ~ 2,
    timepoint == "T2" ~ 6,
    timepoint == "T3" ~ 11,
    timepoint == "T4" ~ 13,
    timepoint == "T5" ~ 27,
    TRUE ~ NA_real_
  )
}

summarise_toc_workbook <- function(path, sheet = "sample_only_mM") {
  df <- read_excel(path, sheet = sheet) %>%
    mutate(sample_id = as.character(sample_id))

  value_col <- if ("toc_mM_actual_nonneg" %in% names(df)) {
    "toc_mM_actual_nonneg"
  } else if ("toc_mM_actual" %in% names(df)) {
    "toc_mM_actual"
  } else if ("toc_mM_nonneg" %in% names(df)) {
    "toc_mM_nonneg"
  } else {
    "toc_mM"
  }

  df %>%
    extract(
      sample_id,
      into = c("source_type", "strain", "oxygen", "replicate", "timepoint"),
      regex = "^(Fresh|Recyc)_(Arthro|Pseudo)_7(Aer|Ana)(?:_(\\d+))?_(T\\d+)$",
      remove = FALSE
    ) %>%
    filter(
      !is.na(source_type),
      !is.na(strain),
      oxygen == "Aer",
      timepoint %in% c("T0", "T1", "T2", "T3", "T4", "T5")
    ) %>%
    mutate(
      day = day_from_timepoint(timepoint),
      source_type = factor(source_type, levels = c("Fresh", "Recyc")),
      source_carbon = factor(
        if_else(strain == "Arthro", "Arthrobacter", "Pseudomonas"),
        levels = c("Arthrobacter", "Pseudomonas")
      ),
      ppm_c = .data[[value_col]] * ppm_per_mM
    ) %>%
    group_by(source_type, source_carbon, timepoint, day) %>%
    summarise(
      n = sum(!is.na(ppm_c)),
      mean = mean(ppm_c, na.rm = TRUE),
      sd = ifelse(n > 1, sd(ppm_c, na.rm = TRUE), 0),
      se = ifelse(n > 1, sd / sqrt(n), 0),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
    ) %>%
    arrange(source_type, source_carbon, day)
}

toc_legacy_summary <- function(path) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(
      pH == 7,
      oxygen == "Aer",
      state %in% c("Fresh", "Recyc"),
      strain %in% c("Arthro", "Pseudo")
    ) %>%
    mutate(
      source_type = factor(state, levels = c("Fresh", "Recyc")),
      source_carbon = factor(
        if_else(strain == "Arthro", "Arthrobacter", "Pseudomonas"),
        levels = c("Arthrobacter", "Pseudomonas")
      ),
      day = day_from_timepoint(timepoint),
      mean = mean_ppmC,
      sd = sd_ppmC,
      se = se_ppmC,
      ymin = mean_ppmC - se_ppmC,
      ymax = mean_ppmC + se_ppmC
    ) %>%
    select(source_type, source_carbon, timepoint, day, mean, sd, se, ymin, ymax) %>%
    arrange(source_type, source_carbon, day)
}

fcm_path <- file.path(PROJECT_DIR, "outputs", "manuscript_story_analysis", "physiology_fcm_gctcd_ph7", "fcm_clean_all_samples.csv")
co2_path <- file.path(PROJECT_DIR, "outputs", "manuscript_story_analysis", "physiology_fcm_gctcd_ph7", "gctcd_clean_all_samples_with_co2.csv")
toc_path <- file.path(PROJECT_DIR, "outputs", "toc_plate_calibration", "calibrated_summary_by_group.csv")
toc_clean_path <- file.path("/Users", "mingfeichen", "Recyc_Necromass_TOC_clean_070326.xlsx")

fcm_raw <- read_csv(fcm_path, show_col_types = FALSE)
co2_raw <- read_csv(co2_path, show_col_types = FALSE)
toc_raw <- read_csv(toc_path, show_col_types = FALSE)
toc_t0_new <- summarise_toc_workbook(toc_clean_path, sheet = "Sheet2") %>%
  filter(timepoint == "T0")
toc_new_all <- summarise_toc_workbook(toc_clean_path, sheet = "Sheet2")

normalize_toc_to_t0 <- function(df) {
  baseline <- df %>%
    filter(timepoint == "T0") %>%
    select(source_type, source_carbon, baseline_mean = mean, baseline_se = se)

  df %>%
    left_join(baseline, by = c("source_type", "source_carbon")) %>%
    mutate(
      mean = mean - baseline_mean,
      se = case_when(
        timepoint == "T0" ~ 0,
        is.na(baseline_se) ~ se,
        TRUE ~ sqrt(se^2 + baseline_se^2)
      ),
      ymin = mean - se,
      ymax = mean + se
    ) %>%
    select(source_type, source_carbon, timepoint, day, mean, sd, se, ymin, ymax)
}

toc_delta_all <- normalize_toc_to_t0(toc_new_all)

fcm <- fcm_raw %>%
  mutate(
    timepoint = as.character(timepoint),
    oxygen = as.character(oxygen),
    source_type = as.character(source_type),
    source_carbon = as.character(source_carbon),
    day = as.numeric(day)
  ) %>%
  filter(
    pH == 7,
    oxygen == "Aerobic",
    source_type %in% c("Fresh", "Recyc"),
    source_carbon %in% c("Arthrobacter", "Pseudomonas"),
    timepoint %in% c("T1", "T2", "T3", "T4", "T5")
  )

co2 <- co2_raw %>%
  mutate(
    timepoint = as.character(timepoint),
    oxygen = as.character(oxygen),
    source_type = as.character(source_type),
    source_carbon = as.character(source_carbon),
    day = as.numeric(day)
  ) %>%
  filter(
    pH == 7,
    oxygen == "Aerobic",
    source_type %in% c("Fresh", "Recyc"),
    source_carbon %in% c("Arthrobacter", "Pseudomonas"),
    timepoint %in% c("T1", "T2", "T3", "T4", "T5")
  )

toc <- toc_raw %>%
  mutate(
    state = as.character(state),
    strain = as.character(strain),
    pH = as.numeric(pH),
    oxygen = as.character(oxygen),
    timepoint = as.character(timepoint),
    day = day_from_timepoint(timepoint)
  ) %>%
  filter(
    pH == 7,
    oxygen == "Aer",
    state %in% c("Fresh", "Recyc"),
    strain %in% c("Arthro", "Pseudo"),
    timepoint %in% c("T0", "T1", "T2", "T3", "T4", "T5")
  ) %>%
  mutate(
    source_type = factor(state, levels = c("Fresh", "Recyc")),
    source_carbon = factor(
      if_else(strain == "Arthro", "Arthrobacter", "Pseudomonas"),
      levels = c("Arthrobacter", "Pseudomonas")
    )
  )

fcm_t0 <- fcm_raw %>%
  filter(timepoint == "T0") %>%
  summarise(
    n = sum(!is.na(log10_cells_ml)),
    mean = mean(log10_cells_ml, na.rm = TRUE),
    sd = ifelse(n > 1, sd(log10_cells_ml, na.rm = TRUE), 0),
    se = ifelse(n > 1, sd / sqrt(n), 0),
    ymin = mean - se,
    ymax = mean + se,
    .groups = "drop"
  ) %>%
  crossing(source_type = c("Fresh", "Recyc")) %>%
  mutate(day = 0)

fcm_t0_plot <- fcm_t0 %>%
  slice(1)

fcm_sum <- fcm %>%
  group_by(source_type, source_carbon, timepoint, day) %>%
  summary_se("log10_cells_ml") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Arthrobacter", "Pseudomonas"))
  ) %>%
  arrange(source_type, source_carbon, day)

co2_sum <- co2 %>%
  group_by(source_type, source_carbon, timepoint, day) %>%
  summary_se("co2_c_ug_per_ml_liquid_cumulative") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Arthrobacter", "Pseudomonas"))
  ) %>%
  arrange(source_type, source_carbon, day)

toc_sum <- toc %>%
  group_by(source_type, source_carbon, timepoint, day) %>%
  summary_se("mean_ppmC") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Arthrobacter", "Pseudomonas"))
  ) %>%
  arrange(source_type, source_carbon, day)

toc_main <- bind_rows(toc_sum, toc_t0_new) %>%
  arrange(source_type, source_carbon, day)

carbon_cols <- c("Arthrobacter" = "#3F6DA8", "Pseudomonas" = "#D95F02")

metric_limits <- function(df, extra = 0.05) {
  lim <- range(c(df$ymin, df$ymax), na.rm = TRUE)
  pad <- diff(lim) * extra
  c(lim[1] - pad, lim[2] + pad)
}

cell_ylim <- metric_limits(bind_rows(fcm_sum, fcm_t0))
co2_ylim <- metric_limits(co2_sum)
toc_ylim_main <- metric_limits(toc_main)
toc_ylim_new <- metric_limits(toc_new_all)

make_metric_panel <- function(df, metric_title, ylab, y_limits, x_breaks, x_limits, source_type_label, show_t0 = FALSE) {
  p <- ggplot(df, aes(day, mean, color = source_carbon, group = source_carbon)) +
    geom_line(linewidth = 0.85) +
    geom_point(size = 2.8) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.25, linewidth = 0.85, color = "grey20", alpha = 1) +
    scale_color_manual(values = carbon_cols, name = "Source carbon") +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.02, 0.04))) +
    labs(
      title = metric_title,
      subtitle = source_type_label,
      x = NULL,
      y = ylab
    ) +
    theme_phys(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 10, hjust = 0),
      legend.position = "none"
    )

  if (show_t0) {
    p <- p +
      geom_point(
        data = fcm_t0_plot,
        inherit.aes = FALSE,
        aes(x = day, y = mean),
        shape = 21,
        fill = "black",
        color = "black",
        stroke = 0.3,
        size = 4.6
      )
  }

  p
}

cell_fresh <- fcm_sum %>% filter(source_type == "Fresh")
cell_recyc <- fcm_sum %>% filter(source_type == "Recyc")
co2_fresh <- co2_sum %>% filter(source_type == "Fresh")
co2_recyc <- co2_sum %>% filter(source_type == "Recyc")
toc_main_fresh <- toc_main %>% filter(source_type == "Fresh")
toc_main_recyc <- toc_main %>% filter(source_type == "Recyc")
toc_new_fresh <- toc_new_all %>% filter(source_type == "Fresh")
toc_new_recyc <- toc_new_all %>% filter(source_type == "Recyc")
toc_delta_fresh <- toc_delta_all %>% filter(source_type == "Fresh")
toc_delta_recyc <- toc_delta_all %>% filter(source_type == "Recyc")

p_cell_fresh <- make_metric_panel(
  cell_fresh,
  "Cell count",
  expression(log[10]~cells~mL^{-1}),
  cell_ylim,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Fresh",
  show_t0 = TRUE
)
p_cell_recyc <- make_metric_panel(
  cell_recyc,
  "Cell count",
  expression(log[10]~cells~mL^{-1}),
  cell_ylim,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Recyc",
  show_t0 = TRUE
)
p_co2_fresh <- make_metric_panel(
  co2_fresh,
  "CO2",
  expression(mu*g~CO[2]*"-C per mL liquid"),
  co2_ylim,
  c(2, 6, 11, 13, 27),
  c(1.5, 27.5),
  "Fresh"
)
p_co2_recyc <- make_metric_panel(
  co2_recyc,
  "CO2",
  expression(mu*g~CO[2]*"-C per mL liquid"),
  co2_ylim,
  c(2, 6, 11, 13, 27),
  c(1.5, 27.5),
  "Recyc"
)
p_toc_main_fresh <- make_metric_panel(
  toc_main_fresh,
  "TOC",
  "TOC (ppm C)",
  toc_ylim_main,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Fresh"
)
p_toc_main_recyc <- make_metric_panel(
  toc_main_recyc,
  "TOC",
  "TOC (ppm C)",
  toc_ylim_main,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Recyc"
)

combined <- (
  p_co2_fresh | p_co2_recyc
) / (
  p_toc_main_fresh | p_toc_main_recyc
) / (
  p_cell_fresh | p_cell_recyc
) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7.png"), combined, width = 11.5, height = 11.5, dpi = 300)
ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7.pdf"), combined, width = 11.5, height = 11.5)

p_toc_new_fresh <- make_metric_panel(
  toc_new_fresh,
  "TOC",
  "TOC (ppm C)",
  toc_ylim_new,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Fresh"
)
p_toc_new_recyc <- make_metric_panel(
  toc_new_recyc,
  "TOC",
  "TOC (ppm C)",
  toc_ylim_new,
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Recyc"
)

combined_new_toc <- (
  p_co2_fresh | p_co2_recyc
) / (
  p_toc_new_fresh | p_toc_new_recyc
) / (
  p_cell_fresh | p_cell_recyc
) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7_new_toc.png"), combined_new_toc, width = 11.5, height = 11.5, dpi = 300)
ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7_new_toc.pdf"), combined_new_toc, width = 11.5, height = 11.5)

p_toc_delta_fresh <- make_metric_panel(
  toc_delta_fresh,
  "TOC change from T0",
  expression(Delta~TOC~"(ppm C)"),
  metric_limits(toc_delta_all),
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Fresh"
)
p_toc_delta_recyc <- make_metric_panel(
  toc_delta_recyc,
  "TOC change from T0",
  expression(Delta~TOC~"(ppm C)"),
  metric_limits(toc_delta_all),
  c(0, 2, 6, 11, 13, 27),
  c(-0.5, 27.5),
  "Recyc"
)

combined_delta_toc <- (
  p_co2_fresh | p_co2_recyc
) / (
  p_toc_delta_fresh | p_toc_delta_recyc
) / (
  p_cell_fresh | p_cell_recyc
) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7_delta_toc.png"), combined_delta_toc, width = 11.5, height = 11.5, dpi = 300)
ggsave(file.path(OUT_DIR, "03_combined_aerobic_only_ph7_delta_toc.pdf"), combined_delta_toc, width = 11.5, height = 11.5)

message("Wrote aerobic-only physiology figures to: ", OUT_DIR)
