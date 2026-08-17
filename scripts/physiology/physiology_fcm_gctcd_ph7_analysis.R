suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(broom)
})

out_dir <- "outputs/manuscript_story_analysis/physiology_fcm_gctcd_ph7"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fcm_path <- "/Users/mingfeichen/Recyc_necromass_FCM.xlsx"
gc_path <- "/Users/mingfeichen/Recyc_necromass_GC-TCD.xlsx"

# Incubation/sampling geometry.
bottle_ml <- 50
initial_liquid_ml <- 15
initial_headspace_ml <- bottle_ml - initial_liquid_ml
gas_sample_ml <- 1
liquid_sample_ml <- 3
atmospheric_co2_ppm <- 420
temperature_c <- 25
molar_volume_l_mol <- 0.082057 * (273.15 + temperature_c)
ug_c_per_mol_co2 <- 12.0107e6

time_map <- c(T0 = 0, T1 = 2, T2 = 6, T3 = 11, T4 = 13, T5 = 27)
interval_map <- c(T1 = 2, T2 = 4, T3 = 5, T4 = 2, T5 = 14)

sampling_geometry <- tibble(
  timepoint = names(interval_map),
  day = unname(time_map[timepoint]),
  interval_days = unname(interval_map[timepoint]),
  sample_index = seq_along(timepoint),
  headspace_ml_before_sampling = initial_headspace_ml + (sample_index - 1) * liquid_sample_ml,
  liquid_ml_during_interval = initial_liquid_ml - (sample_index - 1) * liquid_sample_ml
)

parse_sample_id <- function(df) {
  df %>%
    mutate(
      source_type = case_when(
        str_detect(Sample_ID, "^Fresh") ~ "Fresh",
        str_detect(Sample_ID, "^Recyc") ~ "Recyc",
        str_detect(Sample_ID, "^cntrl") ~ "Control",
        Sample_ID == "inoculum" ~ "Inoculum",
        TRUE ~ NA_character_
      ),
      necromass_state = recode(source_type, Fresh = "Fresh", Recyc = "Recyc"),
      source_carbon = case_when(
        str_detect(Sample_ID, "Pseudo") ~ "Pseudomonas",
        str_detect(Sample_ID, "Arthro") ~ "Arthrobacter",
        str_detect(Sample_ID, "pyruva") ~ "Pyruvate control",
        Sample_ID == "inoculum" ~ "Inoculum",
        TRUE ~ NA_character_
      ),
      organism_panel = recode(source_carbon, Arthrobacter = "Arthro", Pseudomonas = "Pseudo"),
      pH = as.numeric(str_match(Sample_ID, "_(6|7)(Aer|Ana)_")[, 2]),
      oxygen = case_when(
        str_detect(Sample_ID, "Aer_") ~ "Aerobic",
        str_detect(Sample_ID, "Ana_") ~ "Anaerobic",
        TRUE ~ NA_character_
      ),
      replicate_label = as.integer(str_match(Sample_ID, "_([0-9]+)$")[, 2])
    )
}

mean_se <- function(df, value_col) {
  df %>%
    summarise(
      n = sum(!is.na(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = sd / sqrt(n),
      ymin = mean - se,
      ymax = mean + se,
      .groups = "drop"
    )
}

fcm_raw <- read_excel(fcm_path, sheet = "Sheet1") %>%
  rename(timepoint = `Sampling point`, cell_count_ml = `Cell Counts`) %>%
  parse_sample_id() %>%
  mutate(
    day = unname(time_map[timepoint]),
    log10_cells_ml = log10(cell_count_ml)
  )

fcm_ph7 <- fcm_raw %>%
  filter(source_type %in% c("Fresh", "Recyc"), pH == 7, timepoint %in% names(time_map))

fcm_t0 <- fcm_raw %>%
  filter(Sample_ID == "inoculum", timepoint == "T0") %>%
  summarise(cell_count_ml = mean(cell_count_ml, na.rm = TRUE), log10_cells_ml = log10(cell_count_ml)) %>%
  mutate(timepoint = "T0", day = 0)

fcm_plot_data <- fcm_ph7 %>%
  select(Sample_ID, source_type, necromass_state, source_carbon, organism_panel,
         oxygen, timepoint, day, cell_count_ml, log10_cells_ml) %>%
  bind_rows(
    expand_grid(
      source_type = c("Fresh", "Recyc"),
      source_carbon = c("Arthrobacter", "Pseudomonas"),
      oxygen = c("Aerobic", "Anaerobic")
    ) %>%
      mutate(
        necromass_state = source_type,
        organism_panel = recode(source_carbon, Arthrobacter = "Arthro", Pseudomonas = "Pseudo"),
        Sample_ID = "inoculum",
        timepoint = fcm_t0$timepoint,
        day = fcm_t0$day,
        cell_count_ml = fcm_t0$cell_count_ml,
        log10_cells_ml = fcm_t0$log10_cells_ml
      )
  )

fcm_summary <- fcm_plot_data %>%
  group_by(source_carbon, organism_panel, source_type, necromass_state, oxygen, timepoint, day) %>%
  mean_se("log10_cells_ml") %>%
  arrange(source_carbon, source_type, oxygen, day)

write_csv(fcm_raw, file.path(out_dir, "fcm_clean_all_samples.csv"))
write_csv(fcm_ph7, file.path(out_dir, "fcm_clean_ph7_necromass_samples.csv"))
write_csv(fcm_summary, file.path(out_dir, "fcm_log10_cell_count_summary_ph7.csv"))

fcm_lm <- fcm_ph7 %>%
  mutate(timepoint = factor(timepoint, levels = names(time_map))) %>%
  filter(timepoint != "T0") %>%
  lm(log10_cells_ml ~ day + source_type + oxygen + source_carbon, data = .) %>%
  anova() %>%
  tidy()
write_csv(fcm_lm, file.path(out_dir, "fcm_log10_cell_count_lm_anova_ph7.csv"))

std <- read_excel(gc_path, sheet = "STD")
std_fit <- lm(PPM ~ 0 + Area, data = std)
std_fit_with_intercept <- lm(PPM ~ Area, data = std)

std_summary <- bind_rows(
  glance(std_fit) %>% mutate(model = "through_zero"),
  glance(std_fit_with_intercept) %>% mutate(model = "with_intercept")
) %>%
  select(model, everything())
std_coefficients <- bind_rows(
  tidy(std_fit) %>% mutate(model = "through_zero"),
  tidy(std_fit_with_intercept) %>% mutate(model = "with_intercept")
) %>%
  select(model, everything())
write_csv(std_summary, file.path(out_dir, "gctcd_std_model_summary.csv"))
write_csv(std_coefficients, file.path(out_dir, "gctcd_std_model_coefficients.csv"))

gc_raw <- read_excel(gc_path, sheet = "GC-TCD") %>%
  rename(timepoint = `Sampling point`, area = Area) %>%
  parse_sample_id() %>%
  left_join(sampling_geometry, by = "timepoint") %>%
  mutate(
    co2_ppm = as.numeric(predict(std_fit, newdata = tibble(Area = area))),
    replacement_co2_ppm = if_else(oxygen == "Aerobic", atmospheric_co2_ppm, 0),
    initial_headspace_co2_ppm = replacement_co2_ppm,
    co2_ppm_ml_before_sampling = co2_ppm * headspace_ml_before_sampling,
    initial_co2_ppm_ml = initial_headspace_co2_ppm * initial_headspace_ml,
    co2_umol_headspace_before_sampling =
      (co2_ppm * 1e-6) * (headspace_ml_before_sampling / 1000) / molar_volume_l_mol * 1e6
  )

gc_raw <- gc_raw %>%
  arrange(Sample_ID, day) %>%
  group_by(Sample_ID) %>%
  mutate(
    previous_ppm = lag(co2_ppm),
    previous_headspace_ml = lag(headspace_ml_before_sampling),
    previous_replacement_co2_ppm = lag(replacement_co2_ppm),
    previous_post_sampling_co2_ppm_ml = previous_ppm * (previous_headspace_ml - gas_sample_ml) +
      previous_replacement_co2_ppm * (gas_sample_ml + liquid_sample_ml),
    expected_start_interval_co2_ppm_ml = if_else(
      row_number() == 1,
      initial_co2_ppm_ml,
      previous_post_sampling_co2_ppm_ml
    ),
    interval_produced_co2_ppm_ml = co2_ppm_ml_before_sampling - expected_start_interval_co2_ppm_ml,
    cumulative_produced_co2_ppm_ml = cumsum(interval_produced_co2_ppm_ml),
    interval_produced_co2_umol =
      (interval_produced_co2_ppm_ml * 1e-6 / 1000) / molar_volume_l_mol * 1e6,
    cumulative_produced_co2_umol =
      (cumulative_produced_co2_ppm_ml * 1e-6 / 1000) / molar_volume_l_mol * 1e6,
    interval_co2_c_ug = interval_produced_co2_umol * 12.0107,
    cumulative_co2_c_ug = cumulative_produced_co2_umol * 12.0107,
    co2_c_ug_per_ml_liquid_interval = interval_co2_c_ug / liquid_ml_during_interval,
    co2_c_ug_ml_day_interval = co2_c_ug_per_ml_liquid_interval / interval_days,
    co2_c_ug_ml_day_since_start = (cumulative_co2_c_ug / initial_liquid_ml) / day,
    co2_c_ug_per_ml_liquid_cumulative = cumulative_co2_c_ug / initial_liquid_ml
  ) %>%
  ungroup()

gc_ph7 <- gc_raw %>%
  filter(source_type %in% c("Fresh", "Recyc"), pH == 7, timepoint %in% names(interval_map))

gc_summary <- gc_ph7 %>%
  group_by(source_carbon, organism_panel, source_type, necromass_state, oxygen, timepoint, day) %>%
  summarise(
    n = n(),
    area_mean = mean(area, na.rm = TRUE),
    area_se = sd(area, na.rm = TRUE) / sqrt(n),
    co2_ppm_mean = mean(co2_ppm, na.rm = TRUE),
    co2_ppm_se = sd(co2_ppm, na.rm = TRUE) / sqrt(n),
    co2_cumulative_ug_ml_mean = mean(co2_c_ug_per_ml_liquid_cumulative, na.rm = TRUE),
    co2_cumulative_ug_ml_se = sd(co2_c_ug_per_ml_liquid_cumulative, na.rm = TRUE) / sqrt(n),
    co2_rate_since_start_mean = mean(co2_c_ug_ml_day_since_start, na.rm = TRUE),
    co2_rate_since_start_se = sd(co2_c_ug_ml_day_since_start, na.rm = TRUE) / sqrt(n),
    co2_rate_interval_mean = mean(co2_c_ug_ml_day_interval, na.rm = TRUE),
    co2_rate_interval_se = sd(co2_c_ug_ml_day_interval, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    co2_rate_since_start_ymin = co2_rate_since_start_mean - co2_rate_since_start_se,
    co2_rate_since_start_ymax = co2_rate_since_start_mean + co2_rate_since_start_se,
    co2_rate_interval_ymin = co2_rate_interval_mean - co2_rate_interval_se,
    co2_rate_interval_ymax = co2_rate_interval_mean + co2_rate_interval_se
  ) %>%
  arrange(source_carbon, source_type, oxygen, day)

write_csv(gc_raw, file.path(out_dir, "gctcd_clean_all_samples_with_co2.csv"))
write_csv(gc_ph7, file.path(out_dir, "gctcd_clean_ph7_necromass_samples_with_co2.csv"))
write_csv(gc_summary, file.path(out_dir, "gctcd_co2_summary_ph7.csv"))

gc_lm_since_start <- gc_ph7 %>%
  lm(co2_c_ug_ml_day_since_start ~ day + source_type + oxygen + source_carbon, data = .) %>%
  anova() %>%
  tidy()
gc_lm_interval <- gc_ph7 %>%
  lm(co2_c_ug_ml_day_interval ~ day + source_type + oxygen + source_carbon, data = .) %>%
  anova() %>%
  tidy()
write_csv(gc_lm_since_start, file.path(out_dir, "gctcd_co2_rate_since_start_lm_anova_ph7.csv"))
write_csv(gc_lm_interval, file.path(out_dir, "gctcd_co2_interval_rate_lm_anova_ph7.csv"))

condition_endpoint_summary <- fcm_summary %>%
  filter(timepoint == "T5") %>%
  select(source_carbon, source_type, oxygen, endpoint_log10_cells_ml = mean, endpoint_log10_cells_se = se) %>%
  left_join(
    gc_summary %>%
      filter(timepoint == "T5") %>%
      select(source_carbon, source_type, oxygen,
             endpoint_co2_rate_since_start = co2_rate_since_start_mean,
             endpoint_co2_rate_since_start_se = co2_rate_since_start_se,
             endpoint_co2_interval_rate = co2_rate_interval_mean,
             endpoint_co2_interval_rate_se = co2_rate_interval_se),
    by = c("source_carbon", "source_type", "oxygen")
  )
write_csv(condition_endpoint_summary, file.path(out_dir, "endpoint_ph7_cell_count_co2_summary.csv"))

plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.45),
    panel.grid.minor = element_line(color = "grey93", linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(size = 24, face = "plain", hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

color_vals <- c(Aerobic = "#1b9e77", Anaerobic = "#d95f02")
shape_vals <- c(Fresh = 16, Recyc = 17)
linetype_vals <- c(Fresh = "solid", Recyc = "dashed")

p_cells <- ggplot(
  fcm_summary,
  aes(x = day, y = mean, color = oxygen, shape = necromass_state,
      linetype = necromass_state, group = interaction(oxygen, necromass_state))
) +
  geom_line(linewidth = 0.85) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0, linewidth = 0.65, color = "black") +
  geom_point(size = 2.8) +
  facet_grid(. ~ organism_panel) +
  scale_color_manual(values = color_vals, name = "oxygen status") +
  scale_shape_manual(values = shape_vals, name = "necromass_state") +
  scale_linetype_manual(values = linetype_vals, name = "necromass_state") +
  scale_x_continuous(breaks = c(0, 10, 20), minor_breaks = seq(0, 30, 5), limits = c(-1.5, 28.5)) +
  labs(title = "Cell growth", x = "Day (numeric)", y = expression(log[10]~cells~mL^{-1})) +
  guides(linetype = "none") +
  plot_theme

p_co2 <- ggplot(
  gc_summary,
  aes(x = day, y = co2_rate_since_start_mean, color = oxygen, shape = necromass_state,
      linetype = necromass_state, group = interaction(oxygen, necromass_state))
) +
  geom_line(linewidth = 0.85) +
  geom_errorbar(aes(ymin = co2_rate_since_start_ymin, ymax = co2_rate_since_start_ymax),
                width = 0, linewidth = 0.65, color = "black") +
  geom_point(size = 2.8) +
  facet_grid(. ~ organism_panel) +
  scale_color_manual(values = color_vals, name = "oxygen status") +
  scale_shape_manual(values = shape_vals, name = "necromass_state") +
  scale_linetype_manual(values = linetype_vals, name = "necromass_state") +
  scale_x_continuous(breaks = c(0, 10, 20), minor_breaks = seq(0, 30, 5), limits = c(-1.5, 28.5)) +
  labs(
    title = "CO2 respired",
    x = "Day (numeric)",
    y = expression(mu*g~CO[2]*"-C per mL liquid per day")
  ) +
  guides(linetype = "none") +
  plot_theme

p_co2_interval <- ggplot(
  gc_summary,
  aes(x = day, y = co2_rate_interval_mean, color = oxygen, shape = necromass_state,
      linetype = necromass_state, group = interaction(oxygen, necromass_state))
) +
  geom_line(linewidth = 0.85) +
  geom_errorbar(
    aes(ymin = co2_rate_interval_ymin, ymax = co2_rate_interval_ymax),
    width = 0, linewidth = 0.65, color = "black"
  ) +
  geom_point(size = 2.8) +
  facet_grid(. ~ organism_panel) +
  scale_color_manual(values = color_vals, name = "oxygen status") +
  scale_shape_manual(values = shape_vals, name = "necromass_state") +
  scale_linetype_manual(values = linetype_vals, name = "necromass_state") +
  scale_x_continuous(breaks = c(0, 10, 20), minor_breaks = seq(0, 30, 5), limits = c(-1.5, 28.5)) +
  labs(
    title = "CO2 interval rate",
    x = "Day (numeric)",
    y = expression(mu*g~CO[2]*"-C per mL liquid per interval day")
  ) +
  guides(linetype = "none") +
  plot_theme

ggsave(file.path(out_dir, "01_fcm_cell_growth_ph7.png"), p_cells, width = 7.2, height = 4.3, dpi = 300)
ggsave(file.path(out_dir, "01_fcm_cell_growth_ph7.pdf"), p_cells, width = 7.2, height = 4.3)
ggsave(file.path(out_dir, "02_gctcd_co2_respired_ph7.png"), p_co2, width = 7.2, height = 4.3, dpi = 300)
ggsave(file.path(out_dir, "02_gctcd_co2_respired_ph7.pdf"), p_co2, width = 7.2, height = 4.3)
ggsave(file.path(out_dir, "03_gctcd_co2_interval_rate_ph7.png"), p_co2_interval, width = 7.2, height = 4.3, dpi = 300)
ggsave(file.path(out_dir, "03_gctcd_co2_interval_rate_ph7.pdf"), p_co2_interval, width = 7.2, height = 4.3)

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- (p_cells + theme(legend.position = "none")) + p_co2 + patchwork::plot_layout(widths = c(1, 1))
  ggsave(file.path(out_dir, "04_fcm_and_co2_ph7_combined.png"), combined, width = 13.5, height = 4.8, dpi = 300)
  ggsave(file.path(out_dir, "04_fcm_and_co2_ph7_combined.pdf"), combined, width = 13.5, height = 4.8)
}

conversion_notes <- tibble(
  parameter = c("standard_curve_model", "standard_curve_slope_ppm_per_area", "headspace_ml", "liquid_ml", "temperature_c", "molar_volume_l_mol", "main_co2_rate"),
  value = c(
    "PPM ~ 0 + Area",
    as.character(unname(coef(std_fit)[["Area"]])),
    paste0("starts at ", initial_headspace_ml, " mL; increases by ", liquid_sample_ml, " mL after each liquid sample"),
    paste0("starts at ", initial_liquid_ml, " mL; decreases by ", liquid_sample_ml, " mL after each liquid sample"),
    as.character(temperature_c),
    as.character(molar_volume_l_mol),
    "cumulative CO2-C corrected for repeated gas sampling and replacement, divided by initial liquid volume and day since inoculation"
  )
)
write_csv(conversion_notes, file.path(out_dir, "co2_conversion_notes.csv"))

message("Wrote physiology pH7 outputs to: ", out_dir)
