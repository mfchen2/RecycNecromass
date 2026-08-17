#!/usr/bin/env Rscript

# Aerobic-only distance-based redundancy analysis (dbRDA) on Bray-Curtis distances.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- Sys.getenv("R_PLOT_DIR", file.path(QIIME_DIR, "R-plots"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(scales)
  library(grid)
})

theme_recyc <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10),
      plot.margin = margin(8, 26, 8, 8),
      axis.title = element_text(size = 12)
    )
}

save_plot <- function(plot, filename, width = 10, height = 7) {
  ggsave(file.path(OUT_DIR, filename), plot, width = width, height = height, units = "in")
  ggsave(file.path(OUT_DIR, sub("\\.pdf$", ".png", filename)), plot, width = width, height = height, units = "in", dpi = 300)
}

build_physeq_diversity <- function() {
  features_qza <- file.path(QIIME_DIR, "artifacts/table.qza")
  metadata_tsv <- file.path(PROJECT_DIR, "metadata/sample-metadata.tsv")
  stopifnot(file.exists(features_qza), file.exists(metadata_tsv))

  ps <- qza_to_phyloseq(features = features_qza, metadata = metadata_tsv)
  metadata <- data.frame(as(sample_data(ps), "data.frame"), check.names = FALSE)
  if ("oxic.level" %in% names(metadata) && !"oxic level" %in% names(metadata)) {
    names(metadata)[names(metadata) == "oxic.level"] <- "oxic level"
  }
  metadata$sample_id <- rownames(metadata)
  metadata <- metadata %>%
    mutate(
      necromass_type = as.character(necromass_type),
      necromass_source = as.character(necromass_source),
      `oxic level` = as.character(`oxic level`),
      time_point = as.character(time_point),
      source_type = na_if(necromass_type, "NA"),
      source_carbon = na_if(necromass_source, "NA"),
      oxygen = na_if(`oxic level`, "NA"),
      source_type = if_else(is.na(source_type) & source_carbon %in% c("Inoculum", "SynCom"), source_carbon, source_type),
      time_point = factor(time_point, levels = paste0("T", 0:5)),
      source_type = factor(source_type, levels = c("Inoculum", "SynCom", "Predegrading community", "Fresh", "Recyc", "Control")),
      source_carbon = factor(source_carbon, levels = c("Inoculum", "SynCom", "Pseudomonas", "Arthrobacter", "Pyruvate")),
      oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic")),
      replicate = factor(as.character(replicate)),
      sample_label = paste(time_point, oxygen, replicate, sep = "_")
    )

  metadata_for_physeq <- metadata
  rownames(metadata_for_physeq) <- metadata_for_physeq$sample_id
  metadata_for_physeq$sample_id <- NULL
  sample_data(ps) <- sample_data(metadata_for_physeq)
  ps <- prune_samples(sample_sums(ps) > 0, ps)

  default_depth <- min(sample_sums(ps))
  rarefy_depth <- as.integer(Sys.getenv("RAREFY_DEPTH", default_depth))
  if (is.na(rarefy_depth) || rarefy_depth <= 0) rarefy_depth <- default_depth

  set.seed(111)
  ps_rarefied <- rarefy_even_depth(
    ps,
    sample.size = rarefy_depth,
    rngseed = 111,
    replace = FALSE,
    verbose = FALSE
  )

  subset_samples(
    ps_rarefied,
    source_type %in% c("Fresh", "Recyc") &
      source_carbon %in% c("Pseudomonas", "Arthrobacter") &
      oxygen == "Aerobic" &
      time_point %in% paste0("T", 1:5)
  )
}

if (exists("physeq_diversity", inherits = TRUE) && inherits(get("physeq_diversity", inherits = TRUE), "phyloseq")) {
  physeq_dbrda <- get("physeq_diversity", inherits = TRUE)
} else {
  physeq_dbrda <- build_physeq_diversity()
}

physeq_dbrda <- prune_samples(sample_sums(physeq_dbrda) > 0, physeq_dbrda)
physeq_dbrda <- prune_taxa(taxa_sums(physeq_dbrda) > 0, physeq_dbrda)

metadata_dbrda <- data.frame(as(sample_data(physeq_dbrda), "data.frame"), check.names = FALSE)
metadata_dbrda$sample_id <- rownames(metadata_dbrda)
metadata_dbrda <- metadata_dbrda %>%
  mutate(
    Source = case_when(
      source_type == "Fresh" ~ "Fresh",
      source_type == "Recyc" ~ "Recycled",
      TRUE ~ as.character(source_type)
    ),
    CarbonType = case_when(
      source_carbon == "Pseudomonas" ~ "Pseudomonas",
      source_carbon == "Arthrobacter" ~ "Arthrobacter",
      TRUE ~ as.character(source_carbon)
    ),
    TimeNum = as.numeric(gsub("^T", "", as.character(time_point))),
    SourceRecycled = as.numeric(Source == "Recycled"),
    CarbonTypePseudo = as.numeric(CarbonType == "Pseudomonas"),
    Source = factor(Source, levels = c("Fresh", "Recycled")),
    CarbonType = factor(CarbonType, levels = c("Arthrobacter", "Pseudomonas")),
    time_point = factor(time_point, levels = paste0("T", 0:5))
  )

otu <- as(otu_table(physeq_dbrda), "matrix")
if (taxa_are_rows(physeq_dbrda)) otu <- t(otu)
otu <- otu[metadata_dbrda$sample_id, , drop = FALSE]
otu <- otu[, colSums(otu) > 0, drop = FALSE]

dbrda_fit <- capscale(
  otu ~ SourceRecycled + CarbonTypePseudo + TimeNum,
  data = metadata_dbrda,
  distance = "bray",
  add = TRUE
)

bray_dist <- vegdist(otu, method = "bray")
permanova_terms <- adonis2(
  bray_dist ~ SourceRecycled + CarbonTypePseudo + TimeNum,
  data = metadata_dbrda,
  permutations = 999,
  by = "margin"
)
dbrda_overall <- anova.cca(dbrda_fit, permutations = 999)
dbrda_terms <- anova.cca(dbrda_fit, by = "terms", permutations = 999)
dbrda_axes <- anova.cca(dbrda_fit, by = "axis", permutations = 999)

write_csv(as.data.frame(permanova_terms) %>% rownames_to_column("term"), file.path(OUT_DIR, "dbrda_bray_permanova_margin_aerobic_only.csv"))
write_csv(as.data.frame(dbrda_overall) %>% rownames_to_column("term"), file.path(OUT_DIR, "dbrda_bray_anova_overall_aerobic_only.csv"))
write_csv(as.data.frame(dbrda_terms) %>% rownames_to_column("term"), file.path(OUT_DIR, "dbrda_bray_anova_terms_aerobic_only.csv"))
write_csv(as.data.frame(dbrda_axes) %>% rownames_to_column("axis"), file.path(OUT_DIR, "dbrda_bray_anova_axes_aerobic_only.csv"))

site_scores <- as.data.frame(scores(dbrda_fit, display = "sites", choices = 1:2, scaling = 1))
names(site_scores)[1:2] <- c("CAP1", "CAP2")
site_scores <- site_scores %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata_dbrda, by = "sample_id")
write_csv(site_scores, file.path(OUT_DIR, "dbrda_bray_site_scores_aerobic_only.csv"))

biplot_scores <- as.data.frame(scores(dbrda_fit, display = "bp", choices = 1:2, scaling = 1))
names(biplot_scores)[1:2] <- c("CAP1", "CAP2")
biplot_scores <- biplot_scores %>%
  rownames_to_column("term") %>%
  mutate(
    label = recode(
      term,
      "SourceRecycled" = "SourceTypeRecycled",
      "CarbonTypePseudo" = "CarbonTypePseudo",
      "TimeNum" = "TimeNum"
    )
  )

arrow_scale <- {
  point_radius <- max(sqrt(site_scores$CAP1^2 + site_scores$CAP2^2), na.rm = TRUE)
  arrow_radius <- max(sqrt(biplot_scores$CAP1^2 + biplot_scores$CAP2^2), na.rm = TRUE)
  if (is.finite(point_radius) && is.finite(arrow_radius) && arrow_radius > 0) point_radius * 0.82 / arrow_radius else 1
}

arrow_df <- biplot_scores %>%
  mutate(
    CAP1_end = CAP1 * arrow_scale,
    CAP2_end = CAP2 * arrow_scale,
    label_x = CAP1_end * 1.08,
    label_y = CAP2_end * 1.08,
    label_hjust = if_else(CAP1_end >= 0, 0, 1),
    label_vjust = if_else(CAP2_end >= 0, 0, 1)
  )
write_csv(arrow_df, file.path(OUT_DIR, "dbrda_bray_biplot_arrows_aerobic_only.csv"))

constrained_eig <- dbrda_fit$CCA$eig
cap1_pct <- if (length(constrained_eig) >= 1) 100 * constrained_eig[1] / sum(constrained_eig) else NA_real_
cap2_pct <- if (length(constrained_eig) >= 2) 100 * constrained_eig[2] / sum(constrained_eig) else NA_real_

permanova_p <- as.data.frame(permanova_terms) %>%
  rownames_to_column("term") %>%
  filter(term %in% c("SourceRecycled", "CarbonTypePseudo", "TimeNum")) %>%
  transmute(
    term = recode(
      term,
      "SourceRecycled" = "Source",
      "CarbonTypePseudo" = "Carbon",
      "TimeNum" = "Time"
    ),
    p = `Pr(>F)`
  )
subtitle_text <- permanova_p %>%
  mutate(p_label = if_else(is.na(p), "NA", pvalue(p, accuracy = 0.001, add_p = FALSE))) %>%
  transmute(label = paste0(term, " p=", p_label)) %>%
  pull(label) %>%
  paste(collapse = "; ")

dbrda_plot <- ggplot(site_scores, aes(CAP1, CAP2)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.5) +
  geom_point(aes(color = Source, shape = CarbonType), size = 3.1, alpha = 0.92) +
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = CAP1_end, yend = CAP2_end),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.22, "cm")),
    linewidth = 0.55,
    color = "grey20"
  ) +
  geom_text(
    data = arrow_df,
    aes(x = label_x, y = label_y, label = label, hjust = label_hjust, vjust = label_vjust),
    inherit.aes = FALSE,
    size = 3.8,
    color = "grey10"
  ) +
  scale_color_manual(values = c("Fresh" = "#2166AC", "Recycled" = "#E31A1C"), drop = FALSE) +
  scale_shape_manual(values = c("Arthrobacter" = 16, "Pseudomonas" = 17), drop = FALSE) +
  labs(
    title = "dbRDA analysis with permutation tests (aerobic only)",
    subtitle = subtitle_text,
    x = paste0("CAP1 (", round(cap1_pct, 1), "% constrained variance)"),
    y = paste0("CAP2 (", round(cap2_pct, 1), "% constrained variance)"),
    color = "Source",
    shape = "Carbon"
  ) +
  expand_limits(
    x = range(c(site_scores$CAP1, arrow_df$label_x), na.rm = TRUE),
    y = range(c(site_scores$CAP2, arrow_df$label_y), na.rm = TRUE)
  ) +
  coord_equal(clip = "off") +
  theme_recyc(base_size = 12)

save_plot(dbrda_plot, "03p_dbrda_bray_permanova_aerobic_only.pdf", width = 11, height = 7.5)

message("Wrote aerobic-only dbRDA outputs to: ", OUT_DIR)
