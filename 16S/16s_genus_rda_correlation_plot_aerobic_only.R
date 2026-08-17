#!/usr/bin/env Rscript

# Aerobic-only genus-level RDA correlation plot for 16S data.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- file.path(QIIME_DIR, "R-plots", "16s_genus_rda")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(scales)
  library(ggrepel)
  library(grid)
})

theme_rda <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.35),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 13),
      plot.margin = margin(8, 24, 8, 8),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 12)
    )
}

clean_rank <- function(x, fallback) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  missing <- is.na(x) | x == "" | x == "NA"
  x[missing] <- fallback[missing]
  x
}

make_genus_phyloseq <- function() {
  features_qza <- file.path(QIIME_DIR, "artifacts", "table.qza")
  taxonomy_qza <- file.path(QIIME_DIR, "artifacts", "taxonomy.qza")
  metadata_tsv <- file.path(PROJECT_DIR, "metadata", "sample-metadata.tsv")
  stopifnot(file.exists(features_qza), file.exists(taxonomy_qza), file.exists(metadata_tsv))

  ps <- qza_to_phyloseq(
    features = features_qza,
    taxonomy = taxonomy_qza,
    metadata = metadata_tsv
  )

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
      source_type = if_else(
        is.na(source_type) & source_carbon %in% c("Inoculum", "SynCom"),
        source_carbon,
        source_type
      ),
      time_point = factor(time_point, levels = paste0("T", 0:5)),
      source_type = factor(source_type, levels = c("Inoculum", "Fresh", "Recyc", "SynCom", "Predegrading community", "Control")),
      source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas", "SynCom", "Pyruvate")),
      oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic")),
      replicate = factor(as.character(replicate))
    )

  rownames(metadata) <- metadata$sample_id
  metadata$sample_id <- NULL
  sample_data(ps) <- sample_data(metadata)

  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  tax <- data.frame(as(tax_table(ps_genus), "matrix"), check.names = FALSE)
  family_clean <- clean_rank(tax$Family, "Unclassified family")
  genus_clean <- clean_rank(tax$Genus, paste("Unclassified", family_clean))
  genus_clean[is.na(genus_clean) | genus_clean == "Unclassified NA"] <- "Unclassified"
  taxa_names(ps_genus) <- make.unique(genus_clean)

  ps_genus
}

ps_genus <- make_genus_phyloseq()

meta <- data.frame(as(sample_data(ps_genus), "data.frame"), check.names = FALSE) %>%
  rownames_to_column("sample_id") %>%
  mutate(
    source_type = as.character(source_type),
    source_carbon = as.character(source_carbon),
    oxygen = as.character(oxygen),
    time_point = as.character(time_point),
    SourceType = case_when(
      source_type == "Fresh" ~ "Fresh",
      source_type == "Recyc" ~ "Recycled",
      TRUE ~ NA_character_
    ),
    CarbonType = case_when(
      source_carbon == "Pseudomonas" ~ "Pseudo",
      source_carbon == "Arthrobacter" ~ "Arthro",
      TRUE ~ NA_character_
    ),
    Oxygen = as.character(oxygen),
    TimeNum = as.numeric(gsub("^T", "", time_point)),
    SourceRecycled = as.numeric(SourceType == "Recycled"),
    CarbonTypePseudo = as.numeric(CarbonType == "Pseudo")
  ) %>%
  filter(
    SourceType %in% c("Fresh", "Recycled"),
    CarbonType %in% c("Pseudo", "Arthro"),
    Oxygen == "Aerobic",
    time_point %in% paste0("T", 1:5)
  )

ps_genus <- prune_samples(meta$sample_id, ps_genus)
ps_genus <- prune_samples(sample_sums(ps_genus) > 0, ps_genus)
meta <- meta %>% filter(sample_id %in% sample_names(ps_genus))
meta <- meta[match(sample_names(ps_genus), meta$sample_id), ]
rownames(meta) <- meta$sample_id
sample_data(ps_genus) <- sample_data(meta)

otu <- as(otu_table(ps_genus), "matrix")
if (taxa_are_rows(ps_genus)) otu <- t(otu)
otu <- otu[meta$sample_id, , drop = FALSE]
otu <- otu[, colSums(otu) > 0, drop = FALSE]

hellinger <- decostand(otu, method = "hellinger")

rda_fit <- rda(
  hellinger ~ SourceRecycled + CarbonTypePseudo + TimeNum,
  data = meta
)

anova_overall <- anova.cca(rda_fit, permutations = 999)
anova_terms <- anova.cca(rda_fit, by = "terms", permutations = 999)
anova_axes <- anova.cca(rda_fit, by = "axis", permutations = 999)

write_csv(as.data.frame(anova_overall) %>% rownames_to_column("term"), file.path(OUT_DIR, "16s_genus_rda_anova_overall_aerobic_only.csv"))
write_csv(as.data.frame(anova_terms) %>% rownames_to_column("term"), file.path(OUT_DIR, "16s_genus_rda_anova_terms_aerobic_only.csv"))
write_csv(as.data.frame(anova_axes) %>% rownames_to_column("axis"), file.path(OUT_DIR, "16s_genus_rda_anova_axes_aerobic_only.csv"))

site_scores <- as.data.frame(scores(rda_fit, display = "sites", choices = 1:2, scaling = 2)) %>%
  rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")
names(site_scores)[names(site_scores) == "RDA1"] <- "Axis1"
names(site_scores)[names(site_scores) == "RDA2"] <- "Axis2"
write_csv(site_scores, file.path(OUT_DIR, "16s_genus_rda_site_scores_aerobic_only.csv"))

species_scores <- as.data.frame(scores(rda_fit, display = "species", choices = 1:2, scaling = 2)) %>%
  rownames_to_column("genus") %>%
  rename(RDA1 = 2, RDA2 = 3)
species_scores <- species_scores %>%
  mutate(
    genus = as.character(genus),
    mean_abundance = colMeans(otu[, genus, drop = FALSE], na.rm = TRUE),
    score_mag = sqrt(RDA1^2 + RDA2^2)
  )
write_csv(species_scores, file.path(OUT_DIR, "16s_genus_rda_species_scores_aerobic_only.csv"))

biplot_scores <- as.data.frame(scores(rda_fit, display = "bp", choices = 1:2, scaling = 2)) %>%
  rownames_to_column("term") %>%
  rename(RDA1 = 2, RDA2 = 3) %>%
  mutate(
    label = recode(
      term,
      "SourceRecycled" = "SourceType\nRecycled",
      "CarbonTypePseudo" = "CarbonType\nPseudo",
      "TimeNum" = "Time\npoint 5"
    ),
    score_mag = sqrt(RDA1^2 + RDA2^2)
  )

arrow_scale <- {
  point_radius <- max(sqrt(species_scores$RDA1^2 + species_scores$RDA2^2), na.rm = TRUE)
  arrow_radius <- max(sqrt(biplot_scores$RDA1^2 + biplot_scores$RDA2^2), na.rm = TRUE)
  if (is.finite(point_radius) && is.finite(arrow_radius) && arrow_radius > 0) point_radius * 0.8 / arrow_radius else 1
}

arrow_df <- biplot_scores %>%
  mutate(
    RDA1_end = RDA1 * arrow_scale,
    RDA2_end = RDA2 * arrow_scale,
    label_x = RDA1_end * 1.08,
    label_y = RDA2_end * 1.08,
    label_hjust = if_else(RDA1_end >= 0, 0, 1),
    label_vjust = if_else(RDA2_end >= 0, 0, 1)
  )
write_csv(arrow_df, file.path(OUT_DIR, "16s_genus_rda_arrows_aerobic_only.csv"))

constrained_eig <- rda_fit$CCA$eig
axis1_pct <- if (length(constrained_eig) >= 1) 100 * constrained_eig[1] / sum(constrained_eig) else NA_real_
axis2_pct <- if (length(constrained_eig) >= 2) 100 * constrained_eig[2] / sum(constrained_eig) else NA_real_

selected_genera <- species_scores %>%
  filter(!str_detect(genus, "^Unclassified")) %>%
  mutate(
    mean_abundance_pct = 100 * mean_abundance,
    rank_score = score_mag * 1000 + mean_abundance_pct
  ) %>%
  arrange(desc(rank_score)) %>%
  slice_head(n = 25) %>%
  pull(genus)

label_df <- species_scores %>%
  filter(genus %in% selected_genera) %>%
  mutate(
    label_x = RDA1 * 1.03,
    label_y = RDA2 * 1.03
  )

taxa_df <- species_scores %>%
  mutate(
    label_x = RDA1,
    label_y = RDA2
  )

plot_x_limits <- range(c(taxa_df$label_x, arrow_df$label_x), na.rm = TRUE)
plot_y_limits <- range(c(taxa_df$label_y, arrow_df$label_y), na.rm = TRUE)

subtitle_text <- paste0(
  "Aerobic samples only; Hellinger-transformed genus abundances; ",
  "permutation p overall = ", pvalue(anova_overall$`Pr(>F)`[1], accuracy = 0.001, add_p = FALSE)
)

p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey55", linewidth = 0.45) +
  geom_point(
    data = taxa_df,
    aes(RDA1, RDA2),
    color = "#A63D3D",
    size = 1.9,
    alpha = 0.7
  ) +
  geom_text_repel(
    data = label_df,
    aes(label_x, label_y, label = genus),
    color = "#A63D3D",
    size = 7 / ggplot2::.pt,
    box.padding = 0.3,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.color = NA,
    seed = 123
  ) +
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = RDA1_end, yend = RDA2_end),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.24, "cm")),
    linewidth = 0.6,
    color = "grey22"
  ) +
  geom_text(
    data = arrow_df,
    aes(x = label_x, y = label_y, label = label, hjust = label_hjust, vjust = label_vjust),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "grey22"
  ) +
  scale_x_continuous(breaks = pretty_breaks(5)) +
  scale_y_continuous(breaks = pretty_breaks(5)) +
  labs(
    subtitle = subtitle_text,
    x = paste0("RDA1 (", round(axis1_pct, 1), "% constrained variance)"),
    y = paste0("RDA2 (", round(axis2_pct, 1), "% constrained variance)")
  ) +
  coord_equal(
    xlim = c(plot_x_limits[1] * 1.05, plot_x_limits[2] * 1.05),
    ylim = c(plot_y_limits[1] * 1.05, plot_y_limits[2] * 1.05),
    clip = "off"
  ) +
  theme_rda(base_size = 13)

ggsave(file.path(OUT_DIR, "16s_genus_rda_correlation_plot_aerobic_only.pdf"), p, width = 11, height = 8, units = "in")
ggsave(file.path(OUT_DIR, "16s_genus_rda_correlation_plot_aerobic_only.png"), p, width = 11, height = 8, units = "in", dpi = 300)

message("Wrote aerobic-only genus RDA correlation outputs to: ", OUT_DIR)
