#!/usr/bin/env Rscript

# Regenerate beta diversity plots in a single, non-faceted ordination space.
# Encodings:
#   fill = source type (Fresh/Recyc/Inoculum)
#   shape = source carbon (Arthrobacter/Pseudomonas/Inoculum)
#   outline color = oxygen condition
#   light ellipses = source type x source carbon x oxygen groups

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
      legend.position = "right"
    )
}

save_plot <- function(plot, filename, width = 9, height = 7) {
  ggsave(file.path(OUT_DIR, filename), plot, width = width, height = height, units = "in")
  ggsave(file.path(OUT_DIR, sub("\\.pdf$", ".png", filename)), plot, width = width, height = height, units = "in", dpi = 300)
}

prepare_scores <- function(score_file) {
  read_csv(file.path(OUT_DIR, score_file), show_col_types = FALSE) %>%
    filter(
      (
        source_type %in% c("Fresh", "Recyc") &
          source_carbon %in% c("Arthrobacter", "Pseudomonas") &
          oxygen %in% c("Aerobic", "Anaerobic") &
          time_point %in% paste0("T", 1:5)
      ) |
        (
          source_type == "Inoculum" &
            source_carbon == "Inoculum" &
            time_point == "T0"
        )
    ) %>%
    mutate(
      source_type = factor(source_type, levels = c("Inoculum", "Fresh", "Recyc")),
      source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas")),
      oxygen_plot = if_else(time_point == "T0" & source_carbon == "Inoculum", "Inoculum", as.character(oxygen)),
      oxygen_plot = factor(oxygen_plot, levels = c("Inoculum", "Aerobic", "Anaerobic")),
      time_point = factor(time_point, levels = paste0("T", 0:5)),
      group_label = interaction(source_type, source_carbon, oxygen_plot, drop = TRUE, sep = " | ")
    )
}

plot_combined_categories <- function(score_file, x_col, y_col, output_file, title, with_ellipse = TRUE) {
  score_path <- file.path(OUT_DIR, score_file)
  if (!file.exists(score_path)) return(invisible(NULL))

  plot_df <- prepare_scores(score_file)
  write_csv(plot_df, file.path(OUT_DIR, sub("_scores\\.csv$", "_combined_categories_plot_data.csv", score_file)))

  x_sym <- rlang::sym(x_col)
  y_sym <- rlang::sym(y_col)

  p <- ggplot(plot_df, aes(x = !!x_sym, y = !!y_sym))

  if (with_ellipse) {
    p <- p +
      stat_ellipse(
        data = plot_df %>% filter(source_type != "Inoculum"),
        aes(group = group_label, color = oxygen_plot, linetype = source_type),
        linewidth = 0.45,
        alpha = 0.55
      )
  }

  p <- p +
    geom_point(
      aes(fill = source_type, shape = source_carbon, color = oxygen_plot),
      size = 3.2,
      stroke = 0.9,
      alpha = 0.92
    ) +
    scale_fill_manual(
      values = c("Inoculum" = "#BDBDBD", "Fresh" = "#67A9CF", "Recyc" = "#F4A582"),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c("Inoculum" = "#555555", "Aerobic" = "#1F78B4", "Anaerobic" = "#E31A1C"),
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("Inoculum" = 21, "Arthrobacter" = 22, "Pseudomonas" = 24),
      drop = FALSE
    ) +
    labs(
      title = title,
      x = x_col,
      y = y_col,
      fill = "Source type",
      color = "Oxygen",
      shape = "Source carbon"
    ) +
    theme_recyc()

  if (with_ellipse) {
    p <- p +
      scale_linetype_manual(values = c("Fresh" = "solid", "Recyc" = "dashed")) +
      labs(linetype = "Source type ellipse")
  }

  save_plot(p, output_file)
}

plot_combined_timepoints <- function(score_file, x_col, y_col, output_file, title) {
  score_path <- file.path(OUT_DIR, score_file)
  if (!file.exists(score_path)) return(invisible(NULL))

  plot_df <- prepare_scores(score_file)

  x_sym <- rlang::sym(x_col)
  y_sym <- rlang::sym(y_col)

  p <- ggplot(plot_df, aes(x = !!x_sym, y = !!y_sym)) +
    geom_path(
      data = plot_df %>% filter(source_type != "Inoculum"),
      aes(group = interaction(source_type, source_carbon, oxygen_plot, replicate), color = oxygen_plot, linetype = source_type),
      linewidth = 0.35,
      alpha = 0.35
    ) +
    geom_point(
      aes(shape = time_point, color = oxygen_plot, fill = source_type),
      size = 3,
      stroke = 0.8,
      alpha = 0.92
    ) +
    scale_fill_manual(
      values = c("Inoculum" = "#BDBDBD", "Fresh" = "#67A9CF", "Recyc" = "#F4A582"),
      drop = FALSE
    ) +
    scale_color_manual(
      values = c("Inoculum" = "#555555", "Aerobic" = "#1F78B4", "Anaerobic" = "#E31A1C"),
      drop = FALSE
    ) +
    scale_shape_manual(values = c("T0" = 21, "T1" = 16, "T2" = 17, "T3" = 15, "T4" = 18, "T5" = 8), drop = FALSE) +
    scale_linetype_manual(values = c("Fresh" = "solid", "Recyc" = "dashed")) +
    labs(
      title = title,
      x = x_col,
      y = y_col,
      fill = "Source type",
      color = "Oxygen",
      shape = "Time point",
      linetype = "Source type trajectory"
    ) +
    theme_recyc()

  save_plot(p, output_file)
}

plot_combined_categories("beta_bray_pcoa_scores.csv", "PCoA1", "PCoA2", "03j_beta_bray_pcoa_combined_categories.pdf", "Bray-Curtis PCoA, combined categorical groups")
plot_combined_categories("beta_bray_nmds_scores.csv", "NMDS1", "NMDS2", "03k_beta_bray_nmds_combined_categories.pdf", "Bray-Curtis NMDS, combined categorical groups")
plot_combined_categories("beta_weighted_unifrac_pcoa_scores.csv", "PCoA1", "PCoA2", "03l_beta_weighted_unifrac_pcoa_combined_categories.pdf", "Weighted UniFrac PCoA, combined categorical groups", with_ellipse = FALSE)

plot_combined_timepoints("beta_bray_pcoa_scores.csv", "PCoA1", "PCoA2", "03m_beta_bray_pcoa_combined_timepoints.pdf", "Bray-Curtis PCoA, combined time trajectories")
plot_combined_timepoints("beta_bray_nmds_scores.csv", "NMDS1", "NMDS2", "03n_beta_bray_nmds_combined_timepoints.pdf", "Bray-Curtis NMDS, combined time trajectories")
plot_combined_timepoints("beta_weighted_unifrac_pcoa_scores.csv", "PCoA1", "PCoA2", "03o_beta_weighted_unifrac_pcoa_combined_timepoints.pdf", "Weighted UniFrac PCoA, combined time trajectories")

message("Wrote combined beta diversity plots to: ", OUT_DIR)
