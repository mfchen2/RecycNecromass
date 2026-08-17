#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(grid)
  library(svglite)
  library(ragg)
})

PROJECT_DIR <- "/Users/mingfeichen/Recyc_necromass"
OUT_DIR <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_t0_heatmap")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

safe_row_mean <- function(df, pattern) {
  cols <- names(df)[str_detect(names(df), fixed(pattern))]
  if (length(cols) == 0) {
    stop("No columns matched pattern: ", pattern)
  }
  values <- as.matrix(df[, cols, drop = FALSE])
  means <- rowMeans(values, na.rm = TRUE)
  means[is.nan(means)] <- NA_real_
  means
}

normalize_label <- function(x) {
  dplyr::if_else(
    str_to_lower(x) == "malic acid peak 1",
    "Malic acid",
    x
  )
}

read_targeted_matrix <- function(path, organism_prefix) {
  df <- read_csv(path, show_col_types = FALSE)
  meta <- df %>%
    mutate(label = normalize_label(label)) %>%
    select(label, formula, polarity)

  raw <- df %>%
    select(-label, -formula, -polarity)

  fresh_pattern <- paste0("TxCtrl-100C-", organism_prefix, "-A-fresh-d0-non")
  recyc_pattern <- paste0("TxCtrl-100C-", organism_prefix, "-A-recyc-d0-non")

  tibble(
    label = meta$label,
    formula = meta$formula,
    polarity = meta$polarity,
    fresh_d0 = safe_row_mean(raw, fresh_pattern),
    recyc_d0 = safe_row_mean(raw, recyc_pattern)
  )
}

pseudo <- read_targeted_matrix(
  file.path(PROJECT_DIR, "outputs", "targeted_metabolites_pseudomonas", "pseudomonas_targeted_filtered_working.csv"),
  "Pseu"
)

arthro <- read_targeted_matrix(
  file.path(PROJECT_DIR, "outputs", "targeted_metabolites_arthrobacter", "arthrobacter_targeted_filtered_working.csv"),
  "Arth"
)

combined <- full_join(
  pseudo %>% rename(
    formula_pseudo = formula,
    polarity_pseudo = polarity,
    Pseudomonas_Fresh = fresh_d0,
    Pseudomonas_Recycled = recyc_d0
  ),
  arthro %>% rename(
    formula_arthro = formula,
    polarity_arthro = polarity,
    Arthrobacter_Fresh = fresh_d0,
    Arthrobacter_Recycled = recyc_d0
  ),
  by = "label"
) %>%
  mutate(
    formula = coalesce(formula_pseudo, formula_arthro),
    polarity = coalesce(polarity_pseudo, polarity_arthro)
  ) %>%
  select(label, formula, polarity, Pseudomonas_Fresh, Pseudomonas_Recycled, Arthrobacter_Fresh, Arthrobacter_Recycled)

heat_raw <- combined %>%
  select(-label, -formula, -polarity) %>%
  as.matrix()

heat_raw[is.na(heat_raw)] <- 0
colnames(heat_raw) <- c(
  "Pseudomonas\nFresh",
  "Pseudomonas\nRecycled",
  "Arthrobacter\nFresh",
  "Arthrobacter\nRecycled"
)

heat_log <- log2(heat_raw + 1)
row_var <- apply(heat_log, 1, stats::var)
row_var[is.na(row_var)] <- -Inf

keep_idx <- order(row_var, decreasing = TRUE)
top_n <- min(50, length(keep_idx))
keep_idx <- keep_idx[seq_len(top_n)]

plot_mat <- heat_log[keep_idx, , drop = FALSE]
plot_labels <- combined$label[keep_idx]
plot_meta <- combined %>%
  slice(keep_idx) %>%
  select(label, formula, polarity)

row_z <- t(scale(t(plot_mat)))
row_z[is.na(row_z)] <- 0
rownames(row_z) <- plot_labels

annotation_colors <- list(
  Necromass = c(Fresh = "#8DBCE3", Recycled = "#F0B37E")
)

group_split <- factor(
  c("Pseudomonas", "Pseudomonas", "Arthrobacter", "Arthrobacter"),
  levels = c("Pseudomonas", "Arthrobacter")
)

top_anno <- HeatmapAnnotation(
  Necromass = factor(c("Fresh", "Recycled", "Fresh", "Recycled"), levels = c("Fresh", "Recycled")),
  col = annotation_colors,
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
  simple_anno_size = unit(4, "mm")
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B2182B")
)

ht <- Heatmap(
  row_z,
  name = "row z-score",
  col = col_fun,
  top_annotation = top_anno,
  column_split = group_split,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 5.4),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 0,
  column_gap = unit(3, "mm"),
  heatmap_legend_param = list(
    title = "Row z-score",
    at = c(-2, 0, 2),
    labels = c("-2", "0", "2")
  ),
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  show_heatmap_legend = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.25),
  border = TRUE,
  column_title = "t0 necromass metabolite profiles",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

source_data <- plot_meta %>%
  mutate(
    Pseudomonas_Fresh = plot_mat[, "Pseudomonas\nFresh"],
    Pseudomonas_Recycled = plot_mat[, "Pseudomonas\nRecycled"],
    Arthrobacter_Fresh = plot_mat[, "Arthrobacter\nFresh"],
    Arthrobacter_Recycled = plot_mat[, "Arthrobacter\nRecycled"],
    row_variance = row_var[keep_idx]
  )

write_csv(source_data, file.path(OUT_DIR, "t0_fresh_recycled_heatmap_source_data.csv"))

save_heatmap <- function(filename, device = c("png", "pdf", "svg"), width = 9.2, height = 10.5) {
  device <- match.arg(device)
  if (device == "png") {
    ragg::agg_png(filename, width = width, height = height, units = "in", res = 320, background = "white")
  } else if (device == "pdf") {
    grDevices::cairo_pdf(filename, width = width, height = height, family = "Arial")
  } else {
    svglite::svglite(filename, width = width, height = height)
  }
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

base <- file.path(OUT_DIR, "t0_fresh_recycled_arthro_pseudo_heatmap")
save_heatmap(paste0(base, ".png"), device = "png")
save_heatmap(paste0(base, ".pdf"), device = "pdf")
save_heatmap(paste0(base, ".svg"), device = "svg")

read_enrichment_results <- function(path, organism) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(label = normalize_label(label)) %>%
    mutate(
      organism = organism,
      source = str_match(comparison, "^(fresh|recyc)-(A|Ana)-d(\\d+)-")[, 2],
      oxygen = str_match(comparison, "^(fresh|recyc)-(A|Ana)-d(\\d+)-")[, 3],
      day = paste0("d", str_match(comparison, "^(fresh|recyc)-(A|Ana)-d(\\d+)-")[, 4]),
      comparison_label = paste(
        organism,
        if_else(source == "fresh", "Fresh", "Recycled"),
        if_else(oxygen == "A", "Aerobic", "Anaerobic"),
        day
      )
    )
}

pseudo_enrich <- read_enrichment_results(
  file.path(PROJECT_DIR, "outputs", "targeted_metabolites_pseudomonas", "inoc_vs_d0_results.csv"),
  "Pseudomonas"
)

arthro_enrich <- read_enrichment_results(
  file.path(PROJECT_DIR, "outputs", "targeted_metabolites_arthrobacter", "inoc_vs_d0_results_arthro.csv"),
  "Arthrobacter"
)

comparison_levels <- c(
  paste("Pseudomonas Fresh Aerobic", c("d2", "d6", "d11", "d27")),
  paste("Pseudomonas Fresh Anaerobic", c("d2", "d6", "d11", "d27")),
  paste("Pseudomonas Recycled Aerobic", c("d2", "d6", "d11", "d27")),
  paste("Pseudomonas Recycled Anaerobic", c("d2", "d6", "d11", "d27")),
  paste("Arthrobacter Fresh Aerobic", c("d2", "d6", "d11", "d27")),
  paste("Arthrobacter Fresh Anaerobic", c("d2", "d6", "d11", "d27")),
  paste("Arthrobacter Recycled Aerobic", c("d2", "d6", "d11", "d27")),
  paste("Arthrobacter Recycled Anaerobic", c("d2", "d6", "d11", "d27"))
)

enrichment <- bind_rows(pseudo_enrich, arthro_enrich) %>%
  select(label, formula, polarity, organism, source, oxygen, day, comparison_label, log2FC) %>%
  group_by(label, formula, polarity, organism, source, oxygen, day, comparison_label) %>%
  summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(comparison_label = factor(comparison_label, levels = comparison_levels))

enrich_wide <- enrichment %>%
  select(label, formula, polarity, comparison_label, log2FC) %>%
  pivot_wider(
    names_from = comparison_label,
    values_from = log2FC,
    values_fn = mean,
    values_fill = 0
  )

enrich_mat <- enrich_wide %>%
  select(-label, -formula, -polarity) %>%
  as.matrix()
enrich_mat[is.na(enrich_mat)] <- 0

enrich_var <- apply(enrich_mat, 1, function(x) max(abs(x), na.rm = TRUE))
enrich_keep <- order(enrich_var, decreasing = TRUE)
enrich_top_n <- min(50, length(enrich_keep))
enrich_keep <- enrich_keep[seq_len(enrich_top_n)]

enrich_plot <- enrich_mat[enrich_keep, , drop = FALSE]
rownames(enrich_plot) <- enrich_wide$label[enrich_keep]
colnames(enrich_plot) <- rep(c("d2", "d6", "d11", "d27"), times = 8)
enrich_plot <- pmax(pmin(enrich_plot, 4), -4)

enrich_split <- factor(
  c(
    rep("Pseudomonas | Fresh | Aerobic", 4),
    rep("Pseudomonas | Fresh | Anaerobic", 4),
    rep("Pseudomonas | Recycled | Aerobic", 4),
    rep("Pseudomonas | Recycled | Anaerobic", 4),
    rep("Arthrobacter | Fresh | Aerobic", 4),
    rep("Arthrobacter | Fresh | Anaerobic", 4),
    rep("Arthrobacter | Recycled | Aerobic", 4),
    rep("Arthrobacter | Recycled | Anaerobic", 4)
  ),
  levels = c(
    "Pseudomonas | Fresh | Aerobic",
    "Pseudomonas | Fresh | Anaerobic",
    "Pseudomonas | Recycled | Aerobic",
    "Pseudomonas | Recycled | Anaerobic",
    "Arthrobacter | Fresh | Aerobic",
    "Arthrobacter | Fresh | Anaerobic",
    "Arthrobacter | Recycled | Aerobic",
    "Arthrobacter | Recycled | Anaerobic"
  )
)

day_colors <- c(d2 = "#D0E1F2", d6 = "#9EC5E3", d11 = "#6FA8D6", d27 = "#3182BD")

enrich_slice_info <- tibble(
  Organism = rep(c("Pseudomonas", "Pseudomonas", "Pseudomonas", "Pseudomonas",
                   "Arthrobacter", "Arthrobacter", "Arthrobacter", "Arthrobacter"), each = 4),
  Source = rep(c("Fresh", "Fresh", "Recycled", "Recycled",
                 "Fresh", "Fresh", "Recycled", "Recycled"), each = 4),
  Oxygen = rep(c("Aerobic", "Anaerobic", "Aerobic", "Anaerobic",
                 "Aerobic", "Anaerobic", "Aerobic", "Anaerobic"), each = 4),
  Day = rep(c("d2", "d6", "d11", "d27"), times = 8)
)

enrich_factor_colors <- list(
  Organism = c(Pseudomonas = "#4C78A8", Arthrobacter = "#F58518"),
  Source = c(Fresh = "#8DBCE3", Recycled = "#F0B37E"),
  Oxygen = c(Aerobic = "#54A24B", Anaerobic = "#B279A2"),
  Day = day_colors
)

enrich_anno <- HeatmapAnnotation(
  df = enrich_slice_info,
  col = enrich_factor_colors,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(4, "mm")
)

enrich_ht <- Heatmap(
  enrich_plot,
  name = "log2FC",
  col = colorRamp2(c(-4, 0, 4), c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = enrich_anno,
  column_split = enrich_split,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 7.4),
  row_names_max_width = unit(6.8, "cm"),
  column_names_gp = gpar(fontsize = 6.2),
  column_names_rot = 0,
  column_gap = unit(3, "mm"),
  heatmap_legend_param = list(
    title = "log2FC",
    at = c(-4, -2, 0, 2, 4),
    labels = c("-4", "-2", "0", "2", "4")
  ),
  rect_gp = gpar(col = "white", lwd = 0.25),
  border = TRUE,
  column_title = "Incubation-driven metabolite enrichment relative to d0",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

write_csv(
  enrich_wide %>%
    mutate(row_variance = enrich_var[seq_len(nrow(enrich_wide))]),
  file.path(OUT_DIR, "incubation_enrichment_heatmap_source_data.csv")
)

enrich_base <- file.path(OUT_DIR, "incubation_enrichment_shift_arthro_pseudo_heatmap")
save_heatmap_enrich <- function(filename, device = c("png", "pdf", "svg"), width = 13.5, height = 9.5) {
  device <- match.arg(device)
  if (device == "png") {
    ragg::agg_png(filename, width = width, height = height, units = "in", res = 320, background = "white")
  } else if (device == "pdf") {
    grDevices::cairo_pdf(filename, width = width, height = height, family = "Arial")
  } else {
    svglite::svglite(filename, width = width, height = height)
  }
  draw(enrich_ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

save_heatmap_enrich(paste0(enrich_base, ".png"), device = "png", width = 15.5, height = 10.2)
save_heatmap_enrich(paste0(enrich_base, ".pdf"), device = "pdf", width = 15.5, height = 10.2)
save_heatmap_enrich(paste0(enrich_base, ".svg"), device = "svg", width = 15.5, height = 10.2)

message("Wrote heatmaps to: ", base, " and ", enrich_base, ".[png|pdf|svg]")
