suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(broom)
})

input_path <- "/Users/mingfeichen/Downloads/Recyc_Necromass_IC.xlsx"
out_dir <- "outputs/manuscript_story_analysis/ic_nitrate_nitrite"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

time_lookup <- tibble(
  timepoint = c("T1", "T2", "T3", "T4", "T5"),
  day = c(2, 6, 11, 13, 27)
)

ic <- read_excel(input_path, sheet = "Sheet1") %>%
  rename(timepoint = `Sampling point`) %>%
  mutate(
    source_type = str_extract(Sample_ID, "^(Fresh|Recyc|cntrl)"),
    source_type = recode(source_type, cntrl = "Control"),
    source_carbon = case_when(
      str_detect(Sample_ID, "Pseudo") ~ "Pseudomonas",
      str_detect(Sample_ID, "Arthro") ~ "Arthrobacter",
      str_detect(Sample_ID, "pyruva") ~ "Pyruvate control",
      TRUE ~ NA_character_
    ),
    pH = as.integer(str_extract(Sample_ID, "(?<=_)[67](?=Aer|Ana)")),
    oxygen = case_when(
      str_detect(Sample_ID, "Aer") ~ "Aerobic",
      str_detect(Sample_ID, "Ana") ~ "Anaerobic",
      TRUE ~ NA_character_
    ),
    replicate = as.integer(str_extract(Sample_ID, "[0-9]+$")),
    nitrate_mM = Nitrate_ppm / 62.0049,
    nitrite_mM = Nitrite_ppm / 46.0055,
    nox_mM = nitrate_mM + nitrite_mM,
    nitrite_fraction_of_nox = if_else(nox_mM > 0, nitrite_mM / nox_mM, 0)
  ) %>%
  left_join(time_lookup, by = "timepoint") %>%
  mutate(
    amendment_phase = case_when(
      oxygen == "Anaerobic" & timepoint %in% c("T1", "T2", "T3") ~ "Pre-T4 amendment",
      oxygen == "Anaerobic" & timepoint %in% c("T4", "T5") ~ "Post-T4 amendment",
      oxygen == "Aerobic" ~ "Aerobic, no nitrate amendment",
      TRUE ~ "Other"
    ),
    condition_label = paste(source_type, source_carbon, paste0("pH", pH), oxygen, sep = " | ")
  )

write_csv(ic, file.path(out_dir, "ic_clean_parsed.csv"))

sample_id_issues <- ic %>%
  count(Sample_ID, timepoint, name = "rows_per_sample_time") %>%
  filter(rows_per_sample_time > 1) %>%
  arrange(Sample_ID, timepoint)
write_csv(sample_id_issues, file.path(out_dir, "ic_sample_id_duplicates_to_review.csv"))

summary_by_condition <- ic %>%
  group_by(source_type, source_carbon, pH, oxygen, timepoint, day, amendment_phase) %>%
  summarise(
    n = n(),
    nitrate_ppm_mean = mean(Nitrate_ppm, na.rm = TRUE),
    nitrate_ppm_sd = sd(Nitrate_ppm, na.rm = TRUE),
    nitrate_mM_mean = mean(nitrate_mM, na.rm = TRUE),
    nitrate_mM_sd = sd(nitrate_mM, na.rm = TRUE),
    nitrite_ppm_mean = mean(Nitrite_ppm, na.rm = TRUE),
    nitrite_ppm_sd = sd(Nitrite_ppm, na.rm = TRUE),
    nitrite_mM_mean = mean(nitrite_mM, na.rm = TRUE),
    nitrite_mM_sd = sd(nitrite_mM, na.rm = TRUE),
    nox_mM_mean = mean(nox_mM, na.rm = TRUE),
    nox_mM_sd = sd(nox_mM, na.rm = TRUE),
    nitrite_fraction_mean = mean(nitrite_fraction_of_nox, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source_type, source_carbon, pH, oxygen, day)

write_csv(summary_by_condition, file.path(out_dir, "ic_summary_by_condition_timepoint.csv"))

summary_by_oxygen_time <- ic %>%
  group_by(oxygen, timepoint, day, amendment_phase) %>%
  summarise(
    n = n(),
    nitrate_ppm_mean = mean(Nitrate_ppm, na.rm = TRUE),
    nitrate_ppm_sd = sd(Nitrate_ppm, na.rm = TRUE),
    nitrate_mM_mean = mean(nitrate_mM, na.rm = TRUE),
    nitrate_mM_sd = sd(nitrate_mM, na.rm = TRUE),
    nitrite_ppm_mean = mean(Nitrite_ppm, na.rm = TRUE),
    nitrite_ppm_sd = sd(Nitrite_ppm, na.rm = TRUE),
    nitrite_mM_mean = mean(nitrite_mM, na.rm = TRUE),
    nitrite_mM_sd = sd(nitrite_mM, na.rm = TRUE),
    nox_mM_mean = mean(nox_mM, na.rm = TRUE),
    nox_mM_sd = sd(nox_mM, na.rm = TRUE),
    nitrite_fraction_mean = mean(nitrite_fraction_of_nox, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(oxygen, day)

write_csv(summary_by_oxygen_time, file.path(out_dir, "ic_summary_by_oxygen_timepoint.csv"))

endpoint_changes <- summary_by_condition %>%
  filter(oxygen == "Anaerobic", timepoint %in% c("T1", "T3", "T4", "T5")) %>%
  select(source_type, source_carbon, pH, timepoint, nitrate_mM_mean, nitrite_mM_mean, nox_mM_mean, nitrite_fraction_mean) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(nitrate_mM_mean, nitrite_mM_mean, nox_mM_mean, nitrite_fraction_mean)
  ) %>%
  mutate(
    nitrate_delta_T3_minus_T1 = nitrate_mM_mean_T3 - nitrate_mM_mean_T1,
    nitrite_delta_T3_minus_T1 = nitrite_mM_mean_T3 - nitrite_mM_mean_T1,
    nitrate_delta_T5_minus_T4 = nitrate_mM_mean_T5 - nitrate_mM_mean_T4,
    nitrite_delta_T5_minus_T4 = nitrite_mM_mean_T5 - nitrite_mM_mean_T4,
    nox_delta_T5_minus_T4 = nox_mM_mean_T5 - nox_mM_mean_T4
  ) %>%
  arrange(source_carbon, source_type, pH)

write_csv(endpoint_changes, file.path(out_dir, "ic_anaerobic_pre_post_amendment_changes.csv"))

safe_lm <- function(data, formula) {
  fit <- lm(formula, data = data)
  tidy(anova(fit)) %>%
    mutate(term = rownames(anova(fit))) %>%
    select(term, df, sumsq, meansq, statistic, p.value)
}

pre_amend_stats <- ic %>%
  filter(oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3"), source_type != "Control") %>%
  group_by(analyte = "nitrate_mM") %>%
  group_modify(~safe_lm(.x, nitrate_mM ~ day + source_type + source_carbon + factor(pH))) %>%
  ungroup() %>%
  bind_rows(
    ic %>%
      filter(oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3"), source_type != "Control") %>%
      group_by(analyte = "nitrite_mM") %>%
      group_modify(~safe_lm(.x, nitrite_mM ~ day + source_type + source_carbon + factor(pH))) %>%
      ungroup()
  )
write_csv(pre_amend_stats, file.path(out_dir, "ic_pre_T4_anaerobic_lm_stats.csv"))

post_amend_stats <- ic %>%
  filter(oxygen == "Anaerobic", timepoint %in% c("T4", "T5"), source_type != "Control") %>%
  group_by(analyte = "nitrate_mM") %>%
  group_modify(~safe_lm(.x, nitrate_mM ~ day + source_type + source_carbon + factor(pH))) %>%
  ungroup() %>%
  bind_rows(
    ic %>%
      filter(oxygen == "Anaerobic", timepoint %in% c("T4", "T5"), source_type != "Control") %>%
      group_by(analyte = "nitrite_mM") %>%
      group_modify(~safe_lm(.x, nitrite_mM ~ day + source_type + source_carbon + factor(pH))) %>%
      ungroup()
  )
write_csv(post_amend_stats, file.path(out_dir, "ic_post_T4_anaerobic_lm_stats.csv"))

fit_p7_t0 <- function(df, keys) {
  fit <- lm(nitrate_mM ~ day, data = df)
  pred <- predict(fit, newdata = data.frame(day = 0), se.fit = TRUE, interval = "confidence")
  gl <- broom::glance(fit)
  coefs <- summary(fit)$coefficients

  tibble(
    n = nrow(df),
    day_min = min(df$day, na.rm = TRUE),
    day_max = max(df$day, na.rm = TRUE),
    observed_mean_t1 = mean(df$nitrate_mM[df$day == 2], na.rm = TRUE),
    observed_sd_t1 = sd(df$nitrate_mM[df$day == 2], na.rm = TRUE),
    observed_mean_t2 = mean(df$nitrate_mM[df$day == 6], na.rm = TRUE),
    observed_sd_t2 = sd(df$nitrate_mM[df$day == 6], na.rm = TRUE),
    observed_mean_t3 = mean(df$nitrate_mM[df$day == 11], na.rm = TRUE),
    observed_sd_t3 = sd(df$nitrate_mM[df$day == 11], na.rm = TRUE),
    slope_mM_per_day = unname(coef(fit)[["day"]]),
    intercept_mM = unname(coef(fit)[["(Intercept)"]]),
    t0_mM = as.numeric(pred$fit[1]),
    t0_ci_low = as.numeric(pred$fit[2]),
    t0_ci_high = as.numeric(pred$fit[3]),
    t0_se = as.numeric(pred$se.fit[1]),
    r_squared = gl$r.squared,
    adj_r_squared = gl$adj.r.squared,
    p_day = coefs["day", "Pr(>|t|)"]
  )
}

fit_p7_t0_nox <- function(df, keys) {
  fit <- lm(nox_mM ~ day, data = df)
  pred <- predict(fit, newdata = data.frame(day = 0), se.fit = TRUE, interval = "confidence")
  gl <- broom::glance(fit)
  coefs <- summary(fit)$coefficients

  tibble(
    n = nrow(df),
    day_min = min(df$day, na.rm = TRUE),
    day_max = max(df$day, na.rm = TRUE),
    observed_mean_t1 = mean(df$nox_mM[df$day == 2], na.rm = TRUE),
    observed_sd_t1 = sd(df$nox_mM[df$day == 2], na.rm = TRUE),
    observed_mean_t2 = mean(df$nox_mM[df$day == 6], na.rm = TRUE),
    observed_sd_t2 = sd(df$nox_mM[df$day == 6], na.rm = TRUE),
    observed_mean_t3 = mean(df$nox_mM[df$day == 11], na.rm = TRUE),
    observed_sd_t3 = sd(df$nox_mM[df$day == 11], na.rm = TRUE),
    slope_mM_per_day = unname(coef(fit)[["day"]]),
    intercept_mM = unname(coef(fit)[["(Intercept)"]]),
    t0_mM = as.numeric(pred$fit[1]),
    t0_ci_low = as.numeric(pred$fit[2]),
    t0_ci_high = as.numeric(pred$fit[3]),
    t0_se = as.numeric(pred$se.fit[1]),
    r_squared = gl$r.squared,
    adj_r_squared = gl$adj.r.squared,
    p_day = coefs["day", "Pr(>|t|)"]
  )
}

fit_p7_t0_nitrite <- function(df, keys) {
  fit <- lm(nitrite_mM ~ day, data = df)
  pred <- predict(fit, newdata = data.frame(day = 0), se.fit = TRUE, interval = "confidence")
  gl <- broom::glance(fit)
  coefs <- summary(fit)$coefficients

  tibble(
    n = nrow(df),
    day_min = min(df$day, na.rm = TRUE),
    day_max = max(df$day, na.rm = TRUE),
    observed_mean_t1 = mean(df$nitrite_mM[df$day == 2], na.rm = TRUE),
    observed_sd_t1 = sd(df$nitrite_mM[df$day == 2], na.rm = TRUE),
    observed_mean_t2 = mean(df$nitrite_mM[df$day == 6], na.rm = TRUE),
    observed_sd_t2 = sd(df$nitrite_mM[df$day == 6], na.rm = TRUE),
    observed_mean_t3 = mean(df$nitrite_mM[df$day == 11], na.rm = TRUE),
    observed_sd_t3 = sd(df$nitrite_mM[df$day == 11], na.rm = TRUE),
    slope_mM_per_day = unname(coef(fit)[["day"]]),
    intercept_mM = unname(coef(fit)[["(Intercept)"]]),
    t0_mM = as.numeric(pred$fit[1]),
    t0_ci_low = as.numeric(pred$fit[2]),
    t0_ci_high = as.numeric(pred$fit[3]),
    t0_se = as.numeric(pred$se.fit[1]),
    r_squared = gl$r.squared,
    adj_r_squared = gl$adj.r.squared,
    p_day = coefs["day", "Pr(>|t|)"]
  )
}

p7_anaerobic_pre <- ic %>%
  filter(pH == 7, oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3"), source_type != "Control") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Pseudomonas", "Arthrobacter"))
  )

p7_t0_estimates <- p7_anaerobic_pre %>%
  group_by(source_type, source_carbon) %>%
  group_modify(~fit_p7_t0(.x, .y)) %>%
  ungroup() %>%
  arrange(source_carbon, source_type)

write_csv(p7_t0_estimates, file.path(out_dir, "ic_pH7_anaerobic_nitrate_t0_estimates.csv"))

p7_nox_t0_estimates <- p7_anaerobic_pre %>%
  group_by(source_type, source_carbon) %>%
  group_modify(~fit_p7_t0_nox(.x, .y)) %>%
  ungroup() %>%
  arrange(source_carbon, source_type)

write_csv(p7_nox_t0_estimates, file.path(out_dir, "ic_pH7_anaerobic_nox_t0_estimates.csv"))

p7_nitrite_t0_estimates <- p7_anaerobic_pre %>%
  group_by(source_type, source_carbon) %>%
  group_modify(~fit_p7_t0_nitrite(.x, .y)) %>%
  ungroup() %>%
  arrange(source_carbon, source_type)

write_csv(p7_nitrite_t0_estimates, file.path(out_dir, "ic_pH7_anaerobic_nitrite_t0_estimates.csv"))

p7_t0_plot_obs <- summary_by_condition %>%
  filter(pH == 7, oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3"), source_type != "Control") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Pseudomonas", "Arthrobacter"))
  )

p7_t0_alltime_obs <- summary_by_condition %>%
  filter(pH == 7, oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3", "T4", "T5"), source_type != "Control") %>%
  select(source_type, source_carbon, day, nitrate_mM_mean, nitrate_mM_sd) %>%
  bind_rows(
    p7_nox_t0_estimates %>%
      transmute(
        source_type,
        source_carbon,
        day = 0,
        nitrate_mM_mean = t0_mM,
        nitrate_mM_sd = NA_real_,
        t0_ci_low,
        t0_ci_high
      )
  ) %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Pseudomonas", "Arthrobacter"))
  ) %>%
  arrange(source_carbon, source_type, day)

p7_t0_plot <- ggplot(p7_t0_alltime_obs, aes(day, nitrate_mM_mean, color = source_type, group = interaction(source_type, source_carbon))) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  geom_errorbar(
    data = subset(p7_t0_alltime_obs, day != 0),
    aes(ymin = pmax(0, nitrate_mM_mean - nitrate_mM_sd), ymax = nitrate_mM_mean + nitrate_mM_sd),
    width = 0.35,
    alpha = 0.5
  ) +
  geom_linerange(
    data = p7_nox_t0_estimates,
    aes(x = 0, ymin = t0_ci_low, ymax = t0_ci_high, color = source_type),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_point(
    data = p7_nox_t0_estimates,
    aes(x = 0, y = t0_mM, color = source_type),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 3.2,
    stroke = 0.9
  ) +
  facet_wrap(~source_carbon, nrow = 1) +
  scale_color_manual(values = c("Fresh" = "#3F6DA8", "Recyc" = "#B65A44")) +
  scale_x_continuous(breaks = c(0, 2, 6, 11, 13, 27), limits = c(0, 27)) +
  labs(
    x = "Day",
    y = "Nitrate (mM, calculated from ppm NO3-)",
    color = "Necromass"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "04_ic_pH7_anaerobic_nitrate_timecourse_with_t0.png"), p7_t0_plot, width = 10.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "04_ic_pH7_anaerobic_nitrate_timecourse_with_t0.pdf"), p7_t0_plot, width = 10.5, height = 4.8)

p7_nitrite_alltime_obs <- summary_by_condition %>%
  filter(pH == 7, oxygen == "Anaerobic", timepoint %in% c("T1", "T2", "T3", "T4", "T5"), source_type != "Control") %>%
  select(source_type, source_carbon, day, nitrite_mM_mean, nitrite_mM_sd) %>%
  bind_rows(
    p7_nitrite_t0_estimates %>%
      transmute(
        source_type,
        source_carbon,
        day = 0,
        nitrite_mM_mean = 0,
        nitrite_mM_sd = 0,
        t0_ci_low = 0,
        t0_ci_high = 0
      )
  ) %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Pseudomonas", "Arthrobacter"))
  ) %>%
  arrange(source_carbon, source_type, day)

nitrate_peak <- max(c(p7_t0_alltime_obs$nitrate_mM_mean, p7_nox_t0_estimates$t0_mM), na.rm = TRUE)
nitrite_peak <- max(p7_nitrite_alltime_obs$nitrite_mM_mean[p7_nitrite_alltime_obs$day != 0], na.rm = TRUE)
nitrite_scale <- nitrate_peak / nitrite_peak

p7_dual_alltime_obs <- bind_rows(
  p7_t0_alltime_obs %>%
    transmute(
      source_type,
      source_carbon,
      day,
      analyte = "nitrate",
      value_mM = nitrate_mM_mean,
      sd_mM = nitrate_mM_sd,
      ci_low = t0_ci_low,
      ci_high = t0_ci_high
    ),
  p7_nitrite_alltime_obs %>%
    transmute(
      source_type,
      source_carbon,
      day,
      analyte = "nitrite",
      value_mM = nitrite_mM_mean,
      sd_mM = nitrite_mM_sd,
      ci_low = t0_ci_low,
      ci_high = t0_ci_high
    )
) %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Pseudomonas", "Arthrobacter")),
    analyte = factor(analyte, levels = c("nitrate", "nitrite")),
    value_scaled = if_else(analyte == "nitrite", value_mM * nitrite_scale, value_mM),
    sd_scaled = if_else(analyte == "nitrite", sd_mM * nitrite_scale, sd_mM),
    ci_low_scaled = if_else(analyte == "nitrite", ci_low * nitrite_scale, ci_low),
    ci_high_scaled = if_else(analyte == "nitrite", ci_high * nitrite_scale, ci_high)
  ) %>%
  arrange(source_carbon, source_type, analyte, day)

p7_dual_plot <- ggplot(p7_dual_alltime_obs, aes(day, value_scaled, color = source_type, linetype = analyte, group = interaction(source_type, source_carbon, analyte))) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = analyte), size = 2.5, fill = "white") +
  geom_errorbar(
    data = subset(p7_dual_alltime_obs, day != 0),
    aes(ymin = pmax(0, value_scaled - sd_scaled), ymax = value_scaled + sd_scaled),
    width = 0.35,
    alpha = 0.45
  ) +
  geom_linerange(
    data = subset(p7_dual_alltime_obs, day == 0),
    aes(x = day, ymin = ci_low_scaled, ymax = ci_high_scaled, color = source_type),
    inherit.aes = FALSE,
    linewidth = 0.7
  ) +
  geom_point(
    data = subset(p7_dual_alltime_obs, day == 0),
    aes(x = day, y = value_scaled, color = source_type, shape = analyte),
    inherit.aes = FALSE,
    size = 3.1,
    fill = "white",
    stroke = 0.9
  ) +
  facet_wrap(~source_carbon, nrow = 1) +
  scale_color_manual(values = c("Fresh" = "#3F6DA8", "Recyc" = "#B65A44")) +
  scale_linetype_manual(values = c("nitrate" = "solid", "nitrite" = "dashed")) +
  scale_shape_manual(values = c("nitrate" = 21, "nitrite" = 24)) +
  scale_x_continuous(breaks = c(0, 2, 6, 11, 13, 27), limits = c(0, 27)) +
  scale_y_continuous(
    name = "Nitrate (mM)",
    sec.axis = sec_axis(~ . / nitrite_scale, name = "Nitrite (mM)")
  ) +
  labs(
    x = "Day",
    color = "Necromass",
    linetype = "Analyte",
    shape = "Analyte"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "04_ic_pH7_anaerobic_nitrate_nitrite_timecourse_with_t0.png"), p7_dual_plot, width = 11.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "04_ic_pH7_anaerobic_nitrate_nitrite_timecourse_with_t0.pdf"), p7_dual_plot, width = 11.5, height = 4.8)

plot_df <- summary_by_condition %>%
  filter(source_type != "Control") %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic")),
    pH_label = paste0("pH ", pH),
    condition = paste(source_type, source_carbon, sep = " | ")
  )

p_nitrate <- ggplot(plot_df, aes(day, nitrate_mM_mean, color = source_type, linetype = oxygen, group = interaction(source_type, oxygen))) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, nitrate_mM_mean - nitrate_mM_sd), ymax = nitrate_mM_mean + nitrate_mM_sd),
                width = 0.45, alpha = 0.45) +
  facet_grid(pH_label ~ source_carbon, scales = "free_y") +
  scale_color_manual(values = c("Fresh" = "#3F6DA8", "Recyc" = "#B65A44")) +
  labs(
    x = "Day",
    y = "Nitrate (mM, calculated from ppm NO3-)",
    color = "Necromass",
    linetype = "Oxygen"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "01_ic_nitrate_timecourse.png"), p_nitrate, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "01_ic_nitrate_timecourse.pdf"), p_nitrate, width = 10, height = 7)

p_nitrite <- ggplot(plot_df, aes(day, nitrite_mM_mean, color = source_type, linetype = oxygen, group = interaction(source_type, oxygen))) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(0, nitrite_mM_mean - nitrite_mM_sd), ymax = nitrite_mM_mean + nitrite_mM_sd),
                width = 0.45, alpha = 0.45) +
  facet_grid(pH_label ~ source_carbon, scales = "free_y") +
  scale_color_manual(values = c("Fresh" = "#3F6DA8", "Recyc" = "#B65A44")) +
  labs(
    x = "Day",
    y = "Nitrite (mM, calculated from ppm NO2-)",
    color = "Necromass",
    linetype = "Oxygen"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "02_ic_nitrite_timecourse.png"), p_nitrite, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "02_ic_nitrite_timecourse.pdf"), p_nitrite, width = 10, height = 7)

p_nox <- summary_by_oxygen_time %>%
  filter(!is.na(oxygen)) %>%
  ggplot(aes(day, nox_mM_mean, color = oxygen, group = oxygen)) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = pmax(0, nox_mM_mean - nox_mM_sd), ymax = nox_mM_mean + nox_mM_sd),
                width = 0.45, alpha = 0.45) +
  scale_color_manual(values = c("Aerobic" = "#3F6DA8", "Anaerobic" = "#B65A44")) +
  labs(
    x = "Day",
    y = "NOx (mM)",
    color = "Oxygen"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "03_ic_nox_by_oxygen.png"), p_nox, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "03_ic_nox_by_oxygen.pdf"), p_nox, width = 7, height = 5)

message("Wrote IC nitrate/nitrite outputs to: ", out_dir)
