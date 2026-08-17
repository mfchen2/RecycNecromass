#!/usr/bin/env Rscript

# Aerobic-only alpha diversity plots with T0 inoculum overlaid as a reference.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- Sys.getenv("R_PLOT_DIR", file.path(QIIME_DIR, "R-plots"))

suppressPackageStartupMessages({
  library(tidyverse)
})

theme_recyc <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey94", color = "grey70"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggsave(file.path(OUT_DIR, filename), plot, width = width, height = height, units = "in")
}

alpha_df <- read_csv(file.path(OUT_DIR, "alpha_diversity_rarefied.csv"), show_col_types = FALSE) %>%
  mutate(
    source_type = factor(source_type, levels = c("Inoculum", "Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas")),
    time_point = factor(time_point, levels = paste0("T", 0:5))
  )

alpha_plot_df <- bind_rows(
  alpha_df %>%
    filter(source_type %in% c("Fresh", "Recyc"), oxygen == "Aerobic", time_point %in% paste0("T", 1:5)) %>%
    mutate(source_type = as.character(source_type)),
  alpha_df %>%
    filter(source_type == "Inoculum", time_point == "T0", source_carbon == "Inoculum") %>%
    select(-source_type) %>%
    crossing(source_type = c("Fresh", "Recyc"))
) %>%
  mutate(
    source_type = factor(source_type, levels = c("Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas"))
  )

write_csv(alpha_plot_df, file.path(OUT_DIR, "alpha_diversity_plot_data_aerobic_only_with_t0_overlay.csv"))

carbon_cols <- c("Inoculum" = "#666666", "Arthrobacter" = "#3F6DA8", "Pseudomonas" = "#D95F02")

alpha_shannon <- ggplot(alpha_plot_df, aes(x = time_point, y = Shannon, color = source_carbon)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.85) +
  stat_summary(
    data = alpha_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = source_carbon),
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  stat_summary(
    data = alpha_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = source_carbon),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    alpha = 0.7
  ) +
  facet_wrap(~ source_type, nrow = 1) +
  scale_color_manual(values = carbon_cols) +
  labs(x = "Time point", y = "Shannon diversity", color = "Source carbon") +
  theme_recyc()
save_plot(alpha_shannon, "02_alpha_shannon_aerobic_only_with_t0_overlay.pdf", width = 9, height = 5)
save_plot(alpha_shannon, "02_alpha_shannon_aerobic_only_with_t0_overlay.png", width = 9, height = 5)

alpha_plot_long <- alpha_plot_df %>%
  pivot_longer(c(Observed, Shannon, InvSimpson), names_to = "metric", values_to = "value")

alpha_all <- ggplot(alpha_plot_long, aes(x = time_point, y = value, color = source_carbon)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.8, alpha = 0.8) +
  stat_summary(
    data = alpha_plot_long %>% filter(source_carbon != "Inoculum"),
    aes(group = source_carbon),
    fun = mean,
    geom = "line",
    linewidth = 0.7
  ) +
  facet_grid(metric ~ source_type, scales = "free_y") +
  scale_color_manual(values = carbon_cols) +
  labs(x = "Time point", y = NULL, color = "Source carbon") +
  theme_recyc()
save_plot(alpha_all, "02b_alpha_all_metrics_aerobic_only_with_t0_overlay.pdf", width = 10, height = 8)
save_plot(alpha_all, "02b_alpha_all_metrics_aerobic_only_with_t0_overlay.png", width = 10, height = 8)

message("Wrote aerobic-only alpha diversity plots with T0 overlay to: ", OUT_DIR)
