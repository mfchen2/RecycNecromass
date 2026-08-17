#!/usr/bin/env Rscript

# Recyc necromass QIIME2 -> phyloseq analysis
# Outputs are written to qiime2-results-local/R-plots by default.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- Sys.getenv("R_PLOT_DIR", file.path(QIIME_DIR, "R-plots"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

INSTALL_MISSING <- tolower(Sys.getenv("INSTALL_MISSING_R_PACKAGES", "true")) %in% c("true", "1", "yes")

install_if_missing <- function(pkg, source = c("cran", "bioc", "github"), repo = NULL) {
  source <- match.arg(source)
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
  if (!INSTALL_MISSING) {
    stop("Missing R package: ", pkg, ". Install it or set INSTALL_MISSING_R_PACKAGES=true.")
  }
  if (source == "cran") {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else if (source == "bioc") {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else if (source == "github") {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org")
    }
    remotes::install_github(repo, upgrade = "never")
  }
}

install_if_missing("qiime2R", "github", "jbisanz/qiime2R")
install_if_missing("phyloseq", "bioc")
install_if_missing("tidyverse", "cran")
install_if_missing("vegan", "cran")
install_if_missing("ape", "cran")
install_if_missing("scales", "cran")

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(vegan)
  library(ape)
  library(scales)
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

clean_rank <- function(x, fallback) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  x[is.na(x) | x == "" | x == "NA"] <- fallback[is.na(x) | x == "" | x == "NA"]
  x
}

make_sample_label <- function(df) {
  paste(df$time_point, df$oxygen, df$replicate, sep = "_")
}

features_qza <- file.path(QIIME_DIR, "artifacts/table.qza")
taxonomy_qza <- file.path(QIIME_DIR, "artifacts/taxonomy.qza")
tree_qza <- file.path(QIIME_DIR, "artifacts/rooted-tree.qza")
metadata_tsv <- file.path(PROJECT_DIR, "metadata/sample-metadata.tsv")

stopifnot(file.exists(features_qza), file.exists(taxonomy_qza), file.exists(tree_qza), file.exists(metadata_tsv))

physeq <- qza_to_phyloseq(
  features = features_qza,
  taxonomy = taxonomy_qza,
  tree = tree_qza,
  metadata = metadata_tsv
)

metadata <- data.frame(as(sample_data(physeq), "data.frame"), check.names = FALSE)
if ("oxic.level" %in% names(metadata) && !"oxic level" %in% names(metadata)) {
  names(metadata)[names(metadata) == "oxic.level"] <- "oxic level"
}
metadata$sample_id <- rownames(metadata)

# Keep original columns, and add analysis-friendly aliases requested here.
metadata <- metadata %>%
  mutate(
    necromass_type = as.character(necromass_type),
    necromass_source = as.character(necromass_source),
    `oxic level` = as.character(`oxic level`),
    time_point = as.character(time_point),
    source_type = na_if(necromass_type, "NA"),
    source_carbon = na_if(necromass_source, "NA"),
    oxygen = na_if(`oxic level`, "NA"),
    time_point = factor(time_point, levels = paste0("T", 0:5)),
    replicate = as.factor(as.character(replicate)),
    source_type = if_else(is.na(source_type) & source_carbon %in% c("Inoculum", "SynCom"), source_carbon, source_type),
    oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic")),
    source_carbon = factor(source_carbon, levels = c("Inoculum", "SynCom", "Pseudomonas", "Arthrobacter", "Pyruvate")),
    source_type = factor(source_type, levels = c("Inoculum", "SynCom", "Predegrading community", "Fresh", "Recyc", "Control")),
    sample_label = paste(time_point, oxygen, replicate, sep = "_")
  )

metadata_for_physeq <- metadata
rownames(metadata_for_physeq) <- metadata_for_physeq$sample_id
metadata_for_physeq$sample_id <- NULL
sample_data(physeq) <- sample_data(metadata_for_physeq)

# Remove samples with no reads, if any.
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

depths <- tibble(
  sample_id = sample_names(physeq),
  sequencing_depth = sample_sums(physeq)
) %>%
  left_join(metadata, by = "sample_id")

write_csv(depths, file.path(OUT_DIR, "sequencing_depths.csv"))

read_depth_plot <- ggplot(depths, aes(x = reorder(sample_id, sequencing_depth), y = sequencing_depth, fill = source_type)) +
  geom_col(width = 0.85) +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(x = NULL, y = "DADA2 read depth", fill = "Source type") +
  theme_recyc(base_size = 8) +
  theme(axis.text.y = element_text(size = 5))
save_plot(read_depth_plot, "00_read_depths.pdf", width = 8, height = 14)

default_depth <- min(sample_sums(physeq))
rarefy_depth <- as.integer(Sys.getenv("RAREFY_DEPTH", default_depth))
if (is.na(rarefy_depth) || rarefy_depth <= 0) rarefy_depth <- default_depth

set.seed(111)
physeq_rarefied <- rarefy_even_depth(
  physeq,
  sample.size = rarefy_depth,
  rngseed = 111,
  replace = FALSE,
  verbose = TRUE
)

write_lines(
  c(
    paste("Rarefaction depth:", rarefy_depth),
    paste("Samples before rarefaction:", nsamples(physeq)),
    paste("Samples after rarefaction:", nsamples(physeq_rarefied))
  ),
  file.path(OUT_DIR, "rarefaction_depth.txt")
)

physeq_experiment <- subset_samples(
  physeq_rarefied,
  source_type %in% c("Fresh", "Recyc") &
    source_carbon %in% c("Pseudomonas", "Arthrobacter") &
    !is.na(oxygen) &
    !is.na(time_point)
)

physeq_diversity <- subset_samples(
  physeq_rarefied,
  (
    source_type %in% c("Fresh", "Recyc") &
      source_carbon %in% c("Pseudomonas", "Arthrobacter") &
      !is.na(oxygen) &
      !is.na(time_point)
  ) |
    (
      source_type == "Inoculum" &
        source_carbon == "Inoculum" &
        time_point == "T0"
    )
)

physeq_bar <- subset_samples(
  physeq_rarefied,
  time_point != "T0" | source_carbon == "Inoculum"
)
physeq_bar <- prune_taxa(taxa_sums(physeq_bar) > 0, physeq_bar)

physeq_bar_genus <- tax_glom(physeq_bar, taxrank = "Genus", NArm = FALSE)
physeq_bar_rel <- transform_sample_counts(physeq_bar_genus, function(x) x / sum(x) * 100)
bar_df <- psmelt(physeq_bar_rel) %>%
  mutate(
    Family_clean = clean_rank(Family, "Unclassified family"),
    Genus_clean = clean_rank(Genus, paste("Unclassified", Family_clean)),
    Genus_clean = if_else(is.na(Genus_clean) | Genus_clean == "Unclassified NA", "Unclassified", Genus_clean)
  )

genus_keep <- bar_df %>%
  group_by(Genus_clean) %>%
  summarise(max_abundance = max(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(max_abundance >= 1) %>%
  pull(Genus_clean)

bar_df <- bar_df %>%
  mutate(Genus_plot = if_else(Genus_clean %in% genus_keep, Genus_clean, "Other <1%")) %>%
  group_by(
    Sample, sample_label, source_type, source_carbon, oxygen, time_point, replicate, Genus_plot
  ) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  mutate(
    Genus_plot = fct_relevel(Genus_plot, "Other <1%", after = Inf),
    sample_label = factor(sample_label, levels = unique(sample_label[order(time_point, oxygen, replicate)]))
  )

write_csv(bar_df, file.path(OUT_DIR, "genus_relative_abundance_barplot_data.csv"))

genus_bar_sample <- ggplot(bar_df, aes(x = sample_label, y = Abundance, fill = Genus_plot)) +
  geom_col(width = 0.9) +
  facet_grid(source_type ~ source_carbon, scales = "free_x", space = "free_x", drop = TRUE) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Genus") +
  theme_recyc() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
save_plot(genus_bar_sample, "01_genus_barplot_by_sample_source_type_carbon.pdf", width = 15, height = 9)

bar_mean_df <- bar_df %>%
  filter(!is.na(source_type), !is.na(source_carbon)) %>%
  group_by(source_type, source_carbon, oxygen, time_point, Genus_plot) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

genus_bar_mean <- ggplot(bar_mean_df, aes(x = time_point, y = mean_abundance, fill = Genus_plot)) +
  geom_col(width = 0.85) +
  facet_grid(source_type + oxygen ~ source_carbon, scales = "free_x", drop = TRUE) +
  labs(x = NULL, y = "Mean relative abundance (%)", fill = "Genus") +
  theme_recyc()
save_plot(genus_bar_mean, "01b_genus_barplot_mean_by_time_oxygen.pdf", width = 13, height = 10)

# Alpha diversity
alpha_df <- estimate_richness(physeq_diversity, measures = c("Observed", "Shannon", "InvSimpson")) %>%
  rownames_to_column("sample_id") %>%
  left_join(data.frame(as(sample_data(physeq_diversity), "data.frame"), check.names = FALSE) %>% rownames_to_column("sample_id"), by = "sample_id") %>%
  mutate(
    source_type = factor(source_type, levels = c("Inoculum", "Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas")),
    oxygen_plot = if_else(time_point == "T0" & source_carbon == "Inoculum", "Inoculum", as.character(oxygen)),
    oxygen_plot = factor(oxygen_plot, levels = c("Inoculum", "Aerobic", "Anaerobic"))
  )

write_csv(alpha_df, file.path(OUT_DIR, "alpha_diversity_rarefied.csv"))

make_t0_overlay_by_source_type <- function(df) {
  bind_rows(
    df %>%
      filter(source_type %in% c("Fresh", "Recyc")) %>%
      mutate(source_type = as.character(source_type)),
    df %>%
      filter(source_type == "Inoculum", time_point == "T0", source_carbon == "Inoculum") %>%
      select(-source_type) %>%
      crossing(source_type = c("Fresh", "Recyc"))
  ) %>%
    mutate(
      source_type = factor(source_type, levels = c("Fresh", "Recyc")),
      source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas")),
      oxygen_plot = factor(oxygen_plot, levels = c("Inoculum", "Aerobic", "Anaerobic"))
    )
}

alpha_plot_df <- make_t0_overlay_by_source_type(alpha_df)
write_csv(alpha_plot_df, file.path(OUT_DIR, "alpha_diversity_plot_data_with_t0_overlay.csv"))

alpha_shannon <- ggplot(alpha_plot_df, aes(x = time_point, y = Shannon, color = source_carbon, linetype = oxygen_plot)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2, alpha = 0.85) +
  stat_summary(
    data = alpha_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = interaction(source_carbon, oxygen_plot)),
    fun = mean,
    geom = "line",
    linewidth = 0.8
  ) +
  stat_summary(
    data = alpha_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = interaction(source_carbon, oxygen_plot)),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    alpha = 0.7
  ) +
  facet_wrap(~ source_type, nrow = 1) +
  scale_color_manual(values = c("Inoculum" = "#666666", "Arthrobacter" = "#1b9e77", "Pseudomonas" = "#7570b3")) +
  scale_linetype_manual(values = c("Inoculum" = "dotted", "Aerobic" = "solid", "Anaerobic" = "dashed")) +
  labs(x = "Time point", y = "Shannon diversity", color = "Source carbon", linetype = "Oxygen") +
  theme_recyc()
save_plot(alpha_shannon, "02_alpha_shannon_fresh_vs_recyc.pdf", width = 9, height = 5)

alpha_long <- alpha_df %>%
  pivot_longer(c(Observed, Shannon, InvSimpson), names_to = "metric", values_to = "value")

alpha_plot_long <- alpha_plot_df %>%
  pivot_longer(c(Observed, Shannon, InvSimpson), names_to = "metric", values_to = "value")

alpha_all <- ggplot(alpha_plot_long, aes(x = time_point, y = value, color = source_carbon, linetype = oxygen_plot)) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.8, alpha = 0.8) +
  stat_summary(
    data = alpha_plot_long %>% filter(source_carbon != "Inoculum"),
    aes(group = interaction(source_carbon, oxygen_plot)),
    fun = mean,
    geom = "line",
    linewidth = 0.7
  ) +
  facet_grid(metric ~ source_type, scales = "free_y") +
  scale_color_manual(values = c("Inoculum" = "#666666", "Arthrobacter" = "#1b9e77", "Pseudomonas" = "#7570b3")) +
  scale_linetype_manual(values = c("Inoculum" = "dotted", "Aerobic" = "solid", "Anaerobic" = "dashed")) +
  labs(x = "Time point", y = NULL, color = "Source carbon", linetype = "Oxygen") +
  theme_recyc()
save_plot(alpha_all, "02b_alpha_all_metrics_fresh_vs_recyc.pdf", width = 10, height = 8)

alpha_tests <- alpha_long %>%
  filter(source_type %in% c("Fresh", "Recyc")) %>%
  group_by(metric, source_type, time_point, oxygen) %>%
  summarise(
    p_arthro_vs_pseudo = if (n_distinct(source_carbon) == 2) {
      tryCatch(wilcox.test(value ~ source_carbon)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_adj_bh = p.adjust(p_arthro_vs_pseudo, method = "BH"))
write_csv(alpha_tests, file.path(OUT_DIR, "alpha_wilcoxon_arthrobacter_vs_pseudomonas.csv"))

# Beta diversity: Bray PCoA, NMDS, PERMANOVA, ANOSIM, dispersion.
physeq_beta <- prune_taxa(taxa_sums(physeq_diversity) > 0, physeq_diversity)
metadata_beta <- data.frame(as(sample_data(physeq_beta), "data.frame"), check.names = FALSE)
metadata_beta$sample_id <- rownames(metadata_beta)
metadata_beta <- metadata_beta %>%
  mutate(
    source_type = factor(source_type, levels = c("Inoculum", "Fresh", "Recyc")),
    source_carbon = factor(source_carbon, levels = c("Inoculum", "Arthrobacter", "Pseudomonas")),
    oxygen_plot = if_else(time_point == "T0" & source_carbon == "Inoculum", "Inoculum", as.character(oxygen)),
    oxygen_plot = factor(oxygen_plot, levels = c("Inoculum", "Aerobic", "Anaerobic"))
  )

bray_dist <- phyloseq::distance(physeq_beta, method = "bray")

pcoa <- ordinate(physeq_beta, method = "PCoA", distance = bray_dist)
pcoa_var <- round(100 * pcoa$values$Relative_eig[1:2], 1)
pcoa_df <- as.data.frame(pcoa$vectors[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  rename(PCoA1 = Axis.1, PCoA2 = Axis.2) %>%
  left_join(metadata_beta, by = "sample_id")

write_csv(pcoa_df, file.path(OUT_DIR, "beta_bray_pcoa_scores.csv"))

pcoa_plot_df <- make_t0_overlay_by_source_type(pcoa_df)
write_csv(pcoa_plot_df, file.path(OUT_DIR, "beta_bray_pcoa_plot_data_with_t0_overlay.csv"))

pcoa_plot <- ggplot(pcoa_plot_df, aes(PCoA1, PCoA2, color = source_carbon, shape = oxygen_plot)) +
  stat_ellipse(
    data = pcoa_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = interaction(source_type, source_carbon, oxygen), linetype = oxygen),
    linewidth = 0.45,
    alpha = 0.6
  ) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_wrap(~ source_type) +
  scale_color_manual(values = c("Inoculum" = "#666666", "Arthrobacter" = "#1b9e77", "Pseudomonas" = "#7570b3")) +
  scale_shape_manual(values = c("Inoculum" = 4, "Aerobic" = 16, "Anaerobic" = 17)) +
  scale_linetype_manual(values = c("Aerobic" = "solid", "Anaerobic" = "dashed")) +
  labs(
    x = paste0("PCoA1 (", pcoa_var[1], "%)"),
    y = paste0("PCoA2 (", pcoa_var[2], "%)"),
    color = "Source carbon",
    shape = "Oxygen",
    linetype = "Oxygen"
  ) +
  theme_recyc()
save_plot(pcoa_plot, "03_beta_bray_pcoa.pdf", width = 11, height = 5.5)

set.seed(111)
nmds <- metaMDS(bray_dist, k = 2, trymax = 200, trace = FALSE)
nmds_df <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("sample_id") %>%
  left_join(metadata_beta, by = "sample_id")
write_csv(nmds_df, file.path(OUT_DIR, "beta_bray_nmds_scores.csv"))
write_lines(paste("NMDS stress:", round(nmds$stress, 4)), file.path(OUT_DIR, "beta_bray_nmds_stress.txt"))

nmds_plot_df <- make_t0_overlay_by_source_type(nmds_df)
write_csv(nmds_plot_df, file.path(OUT_DIR, "beta_bray_nmds_plot_data_with_t0_overlay.csv"))

nmds_plot <- ggplot(nmds_plot_df, aes(NMDS1, NMDS2, color = source_carbon, shape = oxygen_plot)) +
  stat_ellipse(
    data = nmds_plot_df %>% filter(source_carbon != "Inoculum"),
    aes(group = interaction(source_type, source_carbon, oxygen), linetype = oxygen),
    linewidth = 0.45,
    alpha = 0.6
  ) +
  geom_point(size = 2.8, alpha = 0.9) +
  facet_wrap(~ source_type) +
  scale_color_manual(values = c("Inoculum" = "#666666", "Arthrobacter" = "#1b9e77", "Pseudomonas" = "#7570b3")) +
  scale_shape_manual(values = c("Inoculum" = 4, "Aerobic" = 16, "Anaerobic" = 17)) +
  scale_linetype_manual(values = c("Aerobic" = "solid", "Anaerobic" = "dashed")) +
  labs(color = "Source carbon", shape = "Oxygen", linetype = "Oxygen") +
  theme_recyc()
save_plot(nmds_plot, "03b_beta_bray_nmds.pdf", width = 11, height = 5.5)

physeq_beta_stats <- prune_taxa(taxa_sums(physeq_experiment) > 0, physeq_experiment)
metadata_beta_stats <- data.frame(as(sample_data(physeq_beta_stats), "data.frame"), check.names = FALSE)
metadata_beta_stats$sample_id <- rownames(metadata_beta_stats)
bray_dist_stats <- phyloseq::distance(physeq_beta_stats, method = "bray")

permanova_formula <- bray_dist_stats ~ source_type * source_carbon * oxygen + time_point
permanova <- adonis2(permanova_formula, data = metadata_beta_stats, permutations = 999, by = "margin")
write_csv(as.data.frame(permanova) %>% rownames_to_column("term"), file.path(OUT_DIR, "beta_permanova_bray_margin.csv"))

permanova_terms <- adonis2(bray_dist_stats ~ source_type + source_carbon + oxygen + time_point, data = metadata_beta_stats, permutations = 999, by = "terms")
write_csv(as.data.frame(permanova_terms) %>% rownames_to_column("term"), file.path(OUT_DIR, "beta_permanova_bray_terms.csv"))

group_combo <- interaction(metadata_beta_stats$source_type, metadata_beta_stats$source_carbon, metadata_beta_stats$oxygen, drop = TRUE)
anosim_out <- anosim(bray_dist_stats, group_combo, permutations = 999)
write_lines(
  c(
    paste("ANOSIM grouping: source_type:source_carbon:oxygen"),
    paste("R:", round(anosim_out$statistic, 4)),
    paste("p:", anosim_out$signif)
  ),
  file.path(OUT_DIR, "beta_anosim_bray_source_type_carbon_oxygen.txt")
)

dispersion <- betadisper(bray_dist_stats, group_combo)
dispersion_anova <- anova(dispersion)
write_csv(as.data.frame(dispersion_anova) %>% rownames_to_column("term"), file.path(OUT_DIR, "beta_betadisper_anova.csv"))

pairwise_permanova <- function(distance_object, metadata_df, group_col, strata_cols = NULL) {
  out <- list()
  metadata_df <- as.data.frame(metadata_df)
  if ("sample_id" %in% names(metadata_df)) {
    rownames(metadata_df) <- metadata_df$sample_id
  }
  groups <- na.omit(unique(metadata_df[[group_col]]))
  if (length(groups) < 2) return(tibble())
  pairs <- combn(as.character(groups), 2, simplify = FALSE)
  for (pair in pairs) {
    keep <- metadata_df[[group_col]] %in% pair
    meta_sub <- droplevels(metadata_df[keep, , drop = FALSE])
    dist_sub <- as.dist(as.matrix(distance_object)[rownames(meta_sub), rownames(meta_sub)])
    fit <- adonis2(dist_sub ~ meta_sub[[group_col]], permutations = 999)
    out[[paste(pair, collapse = "_vs_")]] <- tibble(
      comparison = paste(pair, collapse = " vs "),
      r2 = fit$R2[1],
      p = fit$`Pr(>F)`[1]
    )
  }
  bind_rows(out) %>% mutate(p_adj_bh = p.adjust(p, method = "BH"))
}

pairwise_source_type <- metadata_beta_stats %>%
  group_by(source_carbon, oxygen) %>%
  group_modify(~ {
    pairwise_permanova(bray_dist_stats, .x, "source_type")
  }) %>%
  ungroup()
write_csv(pairwise_source_type, file.path(OUT_DIR, "beta_pairwise_permanova_source_type_by_carbon_oxygen.csv"))

# UniFrac PCoA if the tree is available and compatible.
try({
  wunifrac_dist <- phyloseq::distance(physeq_beta, method = "wunifrac")
  wunifrac_pcoa <- ordinate(physeq_beta, method = "PCoA", distance = wunifrac_dist)
  wunifrac_var <- round(100 * wunifrac_pcoa$values$Relative_eig[1:2], 1)
  wunifrac_df <- as.data.frame(wunifrac_pcoa$vectors[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    rename(PCoA1 = Axis.1, PCoA2 = Axis.2) %>%
    left_join(metadata_beta, by = "sample_id")
  write_csv(wunifrac_df, file.path(OUT_DIR, "beta_weighted_unifrac_pcoa_scores.csv"))
  wunifrac_plot_df <- make_t0_overlay_by_source_type(wunifrac_df)
  write_csv(wunifrac_plot_df, file.path(OUT_DIR, "beta_weighted_unifrac_pcoa_plot_data_with_t0_overlay.csv"))
  wunifrac_plot <- ggplot(wunifrac_plot_df, aes(PCoA1, PCoA2, color = source_carbon, shape = oxygen_plot)) +
    geom_point(size = 2.8, alpha = 0.9) +
    facet_wrap(~ source_type) +
    scale_color_manual(values = c("Inoculum" = "#666666", "Arthrobacter" = "#1b9e77", "Pseudomonas" = "#7570b3")) +
    scale_shape_manual(values = c("Inoculum" = 4, "Aerobic" = 16, "Anaerobic" = 17)) +
    labs(
      x = paste0("PCoA1 (", wunifrac_var[1], "%)"),
      y = paste0("PCoA2 (", wunifrac_var[2], "%)"),
      color = "Source carbon",
      shape = "Oxygen"
    ) +
    theme_recyc()
  save_plot(wunifrac_plot, "03c_beta_weighted_unifrac_pcoa.pdf", width = 11, height = 5.5)
}, silent = FALSE)

# Additional relevant summaries.
top_genera <- bar_df %>%
  group_by(Genus_plot) %>%
  summarise(mean_abundance = mean(Abundance), prevalence = mean(Abundance > 0), .groups = "drop") %>%
  arrange(desc(mean_abundance))
write_csv(top_genera, file.path(OUT_DIR, "top_genera_mean_abundance_prevalence.csv"))

heatmap_df <- bar_df %>%
  filter(Genus_plot != "Other <1%") %>%
  semi_join(top_genera %>% slice_head(n = 30), by = "Genus_plot") %>%
  mutate(group = paste(source_type, source_carbon, oxygen, time_point, sep = " | ")) %>%
  group_by(group, Genus_plot) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop")

genus_heatmap <- ggplot(heatmap_df, aes(x = group, y = fct_reorder(Genus_plot, mean_abundance), fill = mean_abundance)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "magma", name = "Mean %") +
  labs(x = NULL, y = "Genus") +
  theme_recyc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_plot(genus_heatmap, "04_top30_genus_heatmap.pdf", width = 14, height = 8)

heatmap_norm_df <- heatmap_df %>%
  group_by(Genus_plot) %>%
  mutate(
    genus_mean = mean(mean_abundance, na.rm = TRUE),
    genus_sd = sd(mean_abundance, na.rm = TRUE),
    abundance_z = if_else(genus_sd > 0, (mean_abundance - genus_mean) / genus_sd, 0),
    abundance_min = min(mean_abundance, na.rm = TRUE),
    abundance_max = max(mean_abundance, na.rm = TRUE),
    abundance_scaled = if_else(
      abundance_max > abundance_min,
      (mean_abundance - abundance_min) / (abundance_max - abundance_min),
      0
    )
  ) %>%
  ungroup()
write_csv(heatmap_norm_df, file.path(OUT_DIR, "top30_genus_heatmap_normalized_data.csv"))

genus_heatmap_z <- ggplot(heatmap_norm_df, aes(x = group, y = fct_reorder(Genus_plot, mean_abundance), fill = abundance_z)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Row z-score"
  ) +
  labs(x = NULL, y = "Genus") +
  theme_recyc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_plot(genus_heatmap_z, "04b_top30_genus_heatmap_row_zscore.pdf", width = 14, height = 8)

genus_heatmap_scaled <- ggplot(heatmap_norm_df, aes(x = group, y = fct_reorder(Genus_plot, mean_abundance), fill = abundance_scaled)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "viridis", name = "Row-scaled\n0-1") +
  labs(x = NULL, y = "Genus") +
  theme_recyc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_plot(genus_heatmap_scaled, "04c_top30_genus_heatmap_row_scaled.pdf", width = 14, height = 8)

top_genera_no_pyruvate <- bar_df %>%
  filter(
    !is.na(source_carbon),
    source_carbon != "Pyruvate",
    Genus_plot != "Other <1%"
  ) %>%
  group_by(Genus_plot) %>%
  summarise(mean_abundance = mean(Abundance), prevalence = mean(Abundance > 0), .groups = "drop") %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 30)
write_csv(top_genera_no_pyruvate, file.path(OUT_DIR, "top_genera_mean_abundance_prevalence_no_pyruvate.csv"))

heatmap_no_pyruvate_df <- bar_df %>%
  filter(
    !is.na(source_carbon),
    source_carbon != "Pyruvate",
    Genus_plot != "Other <1%"
  ) %>%
  semi_join(top_genera_no_pyruvate, by = "Genus_plot") %>%
  mutate(source_carbon = fct_drop(source_carbon)) %>%
  group_by(source_type, source_carbon, oxygen, time_point, Genus_plot) %>%
  summarise(mean_relative_abundance = mean(Abundance), .groups = "drop") %>%
  mutate(group = paste(source_type, source_carbon, oxygen, time_point, sep = " | ")) %>%
  arrange(source_type, source_carbon, oxygen, time_point) %>%
  mutate(group = factor(group, levels = unique(group)))
write_csv(heatmap_no_pyruvate_df, file.path(OUT_DIR, "top30_genus_heatmap_relative_abundance_no_pyruvate_data.csv"))

genus_heatmap_no_pyruvate <- ggplot(
  heatmap_no_pyruvate_df,
  aes(x = group, y = fct_reorder(Genus_plot, mean_relative_abundance), fill = mean_relative_abundance)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_viridis_c(option = "magma", name = "Mean relative\nabundance (%)") +
  labs(x = NULL, y = "Genus") +
  theme_recyc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_plot(genus_heatmap_no_pyruvate, "04d_top30_genus_heatmap_relative_abundance_no_pyruvate.pdf", width = 14, height = 8)

postprocess_scripts <- c(
  "plot-annotated-no-pyruvate-heatmap.R",
  "plot-annotated-clean-heatmaps.R",
  "plot-timepoint-ordinations.R",
  "plot-alpha-t0-overlay.R",
  "plot-beta-t0-overlay.R",
  "plot-beta-combined-categories.R",
  "dbrda-analysis.R",
  "pseudomonas-t1-anaerobic-diversity-tests.R",
  "arthrobacter-t1-recyc-diversity-tests.R",
  "arthrobacter-t1-recyc-oxygen-specific-diversity-tests.R",
  "arthrobacter-t1-anaerobic-recyc-diversity-tests.R"
)
for (postprocess_script in postprocess_scripts) {
  postprocess_path <- file.path(PROJECT_DIR, "scripts", postprocess_script)
  if (file.exists(postprocess_path)) {
    source(postprocess_path, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_DESEQ2", "false")) %in% c("true", "1", "yes")) {
  deseq2_script <- file.path(PROJECT_DIR, "scripts/deseq2-differential-abundance.R")
  if (file.exists(deseq2_script)) {
    source(deseq2_script, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_CONDITION_DESEQ2", "false")) %in% c("true", "1", "yes")) {
  condition_deseq2_script <- file.path(PROJECT_DIR, "scripts/deseq2-condition-time-contrasts.R")
  if (file.exists(condition_deseq2_script)) {
    source(condition_deseq2_script, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_DESEQ2_CONSISTENCY_PLOTS", "false")) %in% c("true", "1", "yes")) {
  consistency_plot_script <- file.path(PROJECT_DIR, "scripts/plot-deseq2-t5-vs-t1-consistency.R")
  if (file.exists(consistency_plot_script)) {
    source(consistency_plot_script, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_ENDPOINT_DESEQ2", "false")) %in% c("true", "1", "yes")) {
  endpoint_deseq2_script <- file.path(PROJECT_DIR, "scripts/deseq2-endpoint-group-enrichment.R")
  if (file.exists(endpoint_deseq2_script)) {
    source(endpoint_deseq2_script, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_ENDPOINT_DESEQ2_PLOTS", "false")) %in% c("true", "1", "yes")) {
  endpoint_plot_script <- file.path(PROJECT_DIR, "scripts/plot-deseq2-endpoint-group-enrichment.R")
  if (file.exists(endpoint_plot_script)) {
    source(endpoint_plot_script, local = new.env(parent = globalenv()))
  }
}

if (tolower(Sys.getenv("RUN_T1_T0_DESEQ2", "false")) %in% c("true", "1", "yes")) {
  t1_t0_deseq2_script <- file.path(PROJECT_DIR, "scripts/deseq2-t1-vs-t0-inoculum.R")
  if (file.exists(t1_t0_deseq2_script)) {
    source(t1_t0_deseq2_script, local = new.env(parent = globalenv()))
  }
}

write_lines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

message("Done. Outputs written to: ", OUT_DIR)
