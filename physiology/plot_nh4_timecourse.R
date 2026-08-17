library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(svglite)
library(ragg)
library(grid)

theme_nature_contract <- function(base_size = 7, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.2),
      plot.caption = element_text(size = base_size - 1, colour = "#4D4D4D"),
      panel.grid = element_blank(),
      legend.position = "top"
    )
}

save_pub_r <- function(plot, filename, width_mm = 183, height_mm = 140, dpi = 600) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4
  h <- height_mm / 25.4

  svglite::svglite(paste0(filename, ".svg"), width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(paste0(filename, ".pdf"), width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_tiff(paste0(filename, ".tiff"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()

  ragg::agg_png(paste0(filename, ".png"), width = w, height = h, units = "in", res = dpi)
  print(plot)
  dev.off()
}

input_xlsx <- "/Users/mingfeichen/Recyc_necromass/Recyc_Necromass_NH4_070226_individual_ugN_mL.xlsx"
out_base <- "/Users/mingfeichen/Recyc_necromass/nh4_lineplots/nh4_timecourse_tocstyle"

df <- readxl::read_xlsx(input_xlsx, sheet = "sample_only_ugN_mL") %>%
  filter(kind == "sample", !is.na(state), !is.na(strain), !is.na(oxygen)) %>%
  mutate(
    timepoint = factor(timepoint, levels = c("T0", "T1", "T2", "T3", "T4", "T5")),
    state = factor(state, levels = c("Fresh", "Recyc")),
    strain = factor(strain, levels = c("Pseudo", "Arthro")),
    oxygen = factor(oxygen, levels = c("Aer", "Ana")),
    replicate = factor(replicate, levels = c(1, 2, 3)),
    replicate_key = ifelse(is.na(replicate), sample_id, as.character(replicate))
  )

# Collapse technical duplicates within each biological replicate, then summarize
# across the three biological replicates for plotting.
replicate_df <- df %>%
  group_by(timepoint, state, strain, oxygen, replicate_key) %>%
  summarise(
    conc = mean(conc_ugN_mL_actual_nonneg, na.rm = TRUE),
    .groups = "drop"
  )

summary_df <- replicate_df %>%
  group_by(timepoint, state, strain, oxygen) %>%
  summarise(
    mean = mean(conc, na.rm = TRUE),
    sd = ifelse(n() > 1, sd(conc, na.rm = TRUE), 0),
    se = sd / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

make_panel <- function(state_level, strain_level, show_x = TRUE, show_y = TRUE, show_legend = FALSE) {
  sdf <- summary_df %>% filter(state == state_level, strain == strain_level)
  rdf <- replicate_df %>% filter(state == state_level, strain == strain_level)

  p <- ggplot(sdf, aes(x = timepoint, y = mean, colour = oxygen, group = oxygen)) +
    geom_point(
      data = rdf,
      aes(x = timepoint, y = conc, colour = oxygen),
      alpha = 0.25, size = 1.6,
      position = position_jitter(width = 0.06, height = 0),
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_errorbar(
      aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
      width = 0.10, linewidth = 0.35, alpha = 0.95
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.8, stroke = 0.2) +
    scale_colour_manual(
      values = c(Aer = "#3182BD", Ana = "#E28E2C"),
      breaks = c("Aer", "Ana"),
      labels = c(Aer = "Aerobic", Ana = "Anaerobic"),
      name = "Oxygen"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme_nature_contract() +
    theme(
      legend.position = if (show_legend) "top" else "none",
      panel.spacing = unit(0.6, "lines"),
      plot.margin = margin(2, 2, 2, 2),
      axis.title.y = element_blank(),
      axis.text.x = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.line.x = if (show_x) element_line() else element_blank(),
      axis.text.y = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.line.y = if (show_y) element_line() else element_blank(),
      strip.background = element_blank(),
      strip.text = element_blank()
    )

  if (!show_x) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
  }
  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())
  }
  p
}

library(cowplot)

label_box <- function(label, rotate = 0, fill = "#F2F2F2") {
  ggdraw() +
    draw_grob(rectGrob(gp = gpar(fill = fill, col = "#D0D0D0", lwd = 0.6)), x = 0, y = 0, width = 1, height = 1) +
    draw_label(label, fontface = "bold", size = 11, angle = rotate)
}

blank_box <- ggdraw() +
  draw_grob(rectGrob(gp = gpar(fill = "white", col = NA)), x = 0, y = 0, width = 1, height = 1)
top_arthro <- label_box("Arthro", rotate = 0)
top_pseudo <- label_box("Pseudo", rotate = 0)
row_fresh <- label_box("Fresh", rotate = 90)
row_recyc <- label_box("Recyc", rotate = 90)

p_legend_src <- make_panel("Fresh", "Arthro", show_x = FALSE, show_y = TRUE, show_legend = TRUE)
legend_grob <- cowplot::get_legend(p_legend_src)
legend_row <- ggdraw() +
  draw_grob(legend_grob, x = 0.02, y = 0.05, width = 0.55, height = 0.90)

p_fresh_arthro <- make_panel("Fresh", "Arthro", show_x = FALSE, show_y = TRUE, show_legend = FALSE)
p_fresh_pseudo <- make_panel("Fresh", "Pseudo", show_x = FALSE, show_y = FALSE, show_legend = FALSE)
p_recyc_arthro <- make_panel("Recyc", "Arthro", show_x = TRUE, show_y = TRUE, show_legend = FALSE)
p_recyc_pseudo <- make_panel("Recyc", "Pseudo", show_x = TRUE, show_y = FALSE, show_legend = FALSE)

header_row <- plot_grid(
  blank_box, top_arthro, top_pseudo, blank_box,
  ncol = 4,
  rel_widths = c(0.14, 1, 1, 0.11),
  align = "h"
)

panel_row_fresh <- plot_grid(
  blank_box, p_fresh_arthro, p_fresh_pseudo, row_fresh,
  ncol = 4,
  rel_widths = c(0.14, 1, 1, 0.11),
  align = "h"
)

panel_row_recyc <- plot_grid(
  blank_box, p_recyc_arthro, p_recyc_pseudo, row_recyc,
  ncol = 4,
  rel_widths = c(0.14, 1, 1, 0.11),
  align = "h"
)

x_label_row <- ggdraw() +
  draw_label("Timepoint", x = 0.50, y = 0.50, size = 12)

fig_core <- plot_grid(
  legend_row,
  header_row,
  panel_row_fresh,
  panel_row_recyc,
  x_label_row,
  ncol = 1,
  rel_heights = c(0.18, 0.10, 1, 1, 0.08),
  align = "v"
)

fig <- ggdraw(fig_core) +
  draw_label("NH4-N (ug/mL)", x = 0.03, y = 0.50, angle = 90, size = 18)

save_pub_r(fig, out_base, width_mm = 183, height_mm = 138, dpi = 600)

message("Saved figure to: ", out_base, ".[svg|pdf|tiff|png]")
message("Replicate summary rows: ", nrow(summary_df))
