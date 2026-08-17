#!/usr/bin/env Rscript

# DESeq2 contrasts for aerobic experiment timepoints.
# Runs pooled aerobic comparisons for:
#   - T1..T5 vs T0 inoculum
#   - T2..T5 vs T1
# Results are filtered to DESeq2 baseMean >= BASEMEAN_MIN.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- Sys.getenv("R_PLOT_DIR", file.path(QIIME_DIR, "R-plots"))
DESEQ_DIR <- file.path(OUT_DIR, "deseq2")
dir.create(DESEQ_DIR, recursive = TRUE, showWarnings = FALSE)

BASEMEAN_MIN <- as.numeric(Sys.getenv("DESEQ_BASEMEAN_MIN", "5"))
if (is.na(BASEMEAN_MIN) || BASEMEAN_MIN < 0) BASEMEAN_MIN <- 5
PADJ_CUTOFF <- as.numeric(Sys.getenv("AEROBIC_TIME_PADJ_CUTOFF", "0.05"))
LFC_CAP <- as.numeric(Sys.getenv("AEROBIC_TIME_LFC_CAP", "10"))

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(DESeq2)
  library(tidyverse)
  library(scales)
})

clean_rank <- function(x, fallback) {
  x <- as.character(x)
  x <- gsub("^[a-z]__", "", x)
  missing <- is.na(x) | x == "" | x == "NA"
  x[missing] <- fallback[missing]
  x
}

theme_aerobic_time <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey94", color = "grey75"),
      legend.position = "right"
    )
}

features_qza <- file.path(QIIME_DIR, "artifacts/table.qza")
taxonomy_qza <- file.path(QIIME_DIR, "artifacts/taxonomy.qza")
metadata_tsv <- file.path(PROJECT_DIR, "metadata/sample-metadata.tsv")
stopifnot(file.exists(features_qza), file.exists(taxonomy_qza), file.exists(metadata_tsv))

physeq <- qza_to_phyloseq(
  features = features_qza,
  taxonomy = taxonomy_qza,
  metadata = metadata_tsv
)

metadata <- data.frame(as(sample_data(physeq), "data.frame"), check.names = FALSE)
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

metadata_for_physeq <- metadata
rownames(metadata_for_physeq) <- metadata_for_physeq$sample_id
metadata_for_physeq$sample_id <- NULL
sample_data(physeq) <- sample_data(metadata_for_physeq)

physeq <- prune_samples(sample_sums(physeq) > 0, physeq)
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

make_tax_annotation <- function(ps, level_name) {
  tax <- data.frame(as(tax_table(ps), "matrix"), check.names = FALSE) %>%
    rownames_to_column("feature_id")
  for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    if (!rank %in% names(tax)) tax[[rank]] <- NA_character_
  }
  tax %>%
    mutate(
      Family_clean = clean_rank(Family, "Unclassified family"),
      Genus_clean = clean_rank(Genus, paste("Unclassified", Family_clean)),
      Genus_clean = if_else(is.na(Genus_clean) | Genus_clean == "Unclassified NA", "Unclassified", Genus_clean),
      feature_level = level_name,
      feature_label = if (level_name == "Genus") Genus_clean else feature_id
    )
}

prepare_genus_phyloseq <- function(ps) {
  ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
  tax <- data.frame(as(tax_table(ps_genus), "matrix"), check.names = FALSE)
  family_clean <- clean_rank(tax$Family, "Unclassified family")
  genus_clean <- clean_rank(tax$Genus, paste("Unclassified", family_clean))
  genus_clean[is.na(genus_clean) | genus_clean == "Unclassified NA"] <- "Unclassified"
  taxa_names(ps_genus) <- make.unique(genus_clean)
  ps_genus
}

time_contrast_grid <- bind_rows(
  tibble(
    comparison_type = "vs_T0_inoculum",
    target_time = paste0("T", 1:5),
    reference_time = "T0",
    reference_label = "T0_inoculum"
  ),
  tibble(
    comparison_type = "vs_T1",
    target_time = paste0("T", 2:5),
    reference_time = "T1",
    reference_label = "T1"
  )
) %>%
  mutate(
    comparison = paste(target_time, "vs", reference_label),
    comparison_label = paste(target_time, "vs", reference_label, "| Aerobic pooled")
  )

source_type_contrast_grid <- tibble(
  comparison_type = "Recyc_vs_Fresh",
  target_time = paste0("T", 1:5),
  reference_time = paste0("T", 1:5),
  reference_label = "Fresh",
  target_label = "Recyc"
) %>%
  mutate(
    comparison = paste(target_label, "vs", reference_label, "at", target_time),
    comparison_label = paste(target_time, "| Recyc vs Fresh | Aerobic")
  )

aerobic_experiment_selector <- function(meta, times) {
  meta$source_type %in% c("Fresh", "Recyc") &
    meta$source_carbon %in% c("Arthrobacter", "Pseudomonas") &
    meta$oxygen == "Aerobic" &
    meta$time_point %in% times
}

reference_selector <- function(meta, reference_label) {
  if (reference_label == "T0_inoculum") {
    meta$source_type == "Inoculum" &
      meta$source_carbon == "Inoculum" &
      meta$time_point == "T0"
  } else if (reference_label == "T1") {
    meta$source_type %in% c("Fresh", "Recyc") &
      meta$source_carbon %in% c("Arthrobacter", "Pseudomonas") &
      meta$oxygen == "Aerobic" &
      meta$time_point == "T1"
  } else {
    stop("Unknown reference label: ", reference_label)
  }
}

mean_normalized_counts <- function(dds, target_level, reference_level) {
  normalized <- counts(dds, normalized = TRUE)
  groups <- as.character(colData(dds)$contrast_group)
  tibble(
    feature_id = rownames(normalized),
    mean_normalized_target = rowMeans(normalized[, groups == target_level, drop = FALSE]),
    mean_normalized_reference = rowMeans(normalized[, groups == reference_level, drop = FALSE])
  )
}

run_single_contrast <- function(ps, feature_level, contrast_row) {
  meta <- data.frame(as(sample_data(ps), "data.frame"), check.names = FALSE)
  keep_samples <- rownames(meta)[
    aerobic_experiment_selector(meta, contrast_row$target_time) |
      reference_selector(meta, contrast_row$reference_label)
  ]

  ps_sub <- prune_samples(keep_samples, ps)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  if (nsamples(ps_sub) < 4 || ntaxa(ps_sub) < 2) return(tibble())

  sub_meta <- data.frame(as(sample_data(ps_sub), "data.frame"), check.names = FALSE)
  sub_meta$contrast_group <- case_when(
    contrast_row$reference_label == "T0_inoculum" &
      as.character(sub_meta$time_point) == "T0" ~ "T0_inoculum",
    contrast_row$reference_label == "T1" &
      as.character(sub_meta$time_point) == "T1" ~ "T1",
    as.character(sub_meta$time_point) == contrast_row$target_time ~ contrast_row$target_time,
    TRUE ~ NA_character_
  )
  sub_meta$contrast_group <- factor(
    sub_meta$contrast_group,
    levels = c(contrast_row$reference_label, contrast_row$target_time)
  )
  rownames(sub_meta) <- sample_names(ps_sub)
  sample_data(ps_sub) <- sample_data(sub_meta)

  group_counts <- table(sample_data(ps_sub)$contrast_group)
  if (any(group_counts[c(contrast_row$reference_label, contrast_row$target_time)] < 2)) return(tibble())

  annotation <- make_tax_annotation(ps_sub, feature_level)
  dds <- phyloseq_to_deseq2(ps_sub, ~ contrast_group)
  dds <- dds[rowSums(counts(dds)) > 0, ]
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)

  result <- results(dds, contrast = c("contrast_group", contrast_row$target_time, contrast_row$reference_label))
  normalized_means <- mean_normalized_counts(dds, contrast_row$target_time, contrast_row$reference_label)

  as.data.frame(result) %>%
    rownames_to_column("feature_id") %>%
    as_tibble() %>%
    left_join(annotation, by = "feature_id") %>%
    left_join(normalized_means, by = "feature_id") %>%
    mutate(
      feature_level = feature_level,
      contrast = paste0(contrast_row$target_time, "_vs_", contrast_row$reference_label),
      comparison_type = contrast_row$comparison_type,
      comparison = contrast_row$comparison,
      comparison_label = contrast_row$comparison_label,
      target_time = contrast_row$target_time,
      reference_time = contrast_row$reference_time,
      reference_label = contrast_row$reference_label,
      experiment_group = "Aerobic pooled",
      n_target = unname(group_counts[contrast_row$target_time]),
      n_reference = unname(group_counts[contrast_row$reference_label]),
      baseMean_cutoff = BASEMEAN_MIN,
      passes_baseMean = !is.na(baseMean) & baseMean >= BASEMEAN_MIN,
      significant_padj_0.05 = !is.na(padj) & padj < PADJ_CUTOFF,
      significant_padj_0.05_lfc1 = !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= 1,
      direction = case_when(
        !passes_baseMean ~ "below_baseMean_cutoff",
        is.na(padj) ~ "not_tested",
        padj < PADJ_CUTOFF & log2FoldChange > 0 ~ paste0("increased_at_", contrast_row$target_time),
        padj < PADJ_CUTOFF & log2FoldChange < 0 ~ paste0("decreased_at_", contrast_row$target_time),
        TRUE ~ "not_significant"
      )
    ) %>%
    filter(passes_baseMean) %>%
    arrange(comparison_type, comparison, padj, desc(abs(log2FoldChange)))
}

run_timepoint_contrasts <- function(ps, feature_level) {
  map_dfr(seq_len(nrow(time_contrast_grid)), function(i) {
    run_single_contrast(ps, feature_level, time_contrast_grid[i, ])
  })
}

run_recyc_vs_fresh_contrasts <- function(ps, feature_level) {
  map_dfr(seq_len(nrow(source_type_contrast_grid)), function(i) {
    contrast_row <- source_type_contrast_grid[i, ]
    meta <- data.frame(as(sample_data(ps), "data.frame"), check.names = FALSE)
    keep_samples <- rownames(meta)[
      meta$source_type %in% c("Fresh", "Recyc") &
        meta$source_carbon %in% c("Arthrobacter", "Pseudomonas") &
        meta$oxygen == "Aerobic" &
        meta$time_point == contrast_row$target_time
    ]

    ps_sub <- prune_samples(keep_samples, ps)
    ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
    if (nsamples(ps_sub) < 4 || ntaxa(ps_sub) < 2) return(tibble())

    sub_meta <- data.frame(as(sample_data(ps_sub), "data.frame"), check.names = FALSE)
    sub_meta$contrast_group <- as.character(sub_meta$source_type)
    sub_meta$contrast_group <- factor(sub_meta$contrast_group, levels = c("Fresh", "Recyc"))
    rownames(sub_meta) <- sample_names(ps_sub)
    sample_data(ps_sub) <- sample_data(sub_meta)

    group_counts <- table(sample_data(ps_sub)$contrast_group)
    if (any(group_counts[c("Fresh", "Recyc")] < 2)) return(tibble())

    annotation <- make_tax_annotation(ps_sub, feature_level)
    dds <- phyloseq_to_deseq2(ps_sub, ~ contrast_group)
    dds <- dds[rowSums(counts(dds)) > 0, ]
    keep <- rowSums(counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    dds <- DESeq(dds, quiet = TRUE)

    result <- results(dds, contrast = c("contrast_group", "Recyc", "Fresh"))
    normalized_means <- mean_normalized_counts(dds, "Recyc", "Fresh")

    as.data.frame(result) %>%
      rownames_to_column("feature_id") %>%
      as_tibble() %>%
      left_join(annotation, by = "feature_id") %>%
      left_join(normalized_means, by = "feature_id") %>%
      mutate(
        feature_level = feature_level,
        contrast = paste0("Recyc_vs_Fresh_at_", contrast_row$target_time),
        comparison_type = contrast_row$comparison_type,
        comparison = contrast_row$comparison,
        comparison_label = contrast_row$comparison_label,
        target_time = contrast_row$target_time,
        reference_time = contrast_row$reference_time,
        reference_label = contrast_row$reference_label,
        target_label = contrast_row$target_label,
        experiment_group = "Aerobic pooled",
        n_target = unname(group_counts["Recyc"]),
        n_reference = unname(group_counts["Fresh"]),
        baseMean_cutoff = BASEMEAN_MIN,
        passes_baseMean = !is.na(baseMean) & baseMean >= BASEMEAN_MIN,
        significant_padj_0.05 = !is.na(padj) & padj < PADJ_CUTOFF,
        significant_padj_0.05_lfc1 = !is.na(padj) & padj < PADJ_CUTOFF & abs(log2FoldChange) >= 1,
        direction = case_when(
          !passes_baseMean ~ "below_baseMean_cutoff",
          is.na(padj) ~ "not_tested",
          padj < PADJ_CUTOFF & log2FoldChange > 0 ~ paste0("increased_at_", contrast_row$target_label),
          padj < PADJ_CUTOFF & log2FoldChange < 0 ~ paste0("decreased_at_", contrast_row$target_label),
          TRUE ~ "not_significant"
        )
      ) %>%
      filter(passes_baseMean) %>%
      arrange(comparison_type, comparison, padj, desc(abs(log2FoldChange)))
  })
}

levels_to_run <- strsplit(Sys.getenv("DESEQ_FEATURE_LEVELS", "Genus,ASV"), ",")[[1]] %>%
  trimws()

all_results <- list()

if ("Genus" %in% levels_to_run) {
  message("Running aerobic timepoint genus DESeq2 contrasts...")
  genus_ps <- prepare_genus_phyloseq(physeq)
  genus_time_results <- run_timepoint_contrasts(genus_ps, "Genus")
  genus_recyc_fresh_results <- run_recyc_vs_fresh_contrasts(genus_ps, "Genus")
  genus_results <- bind_rows(genus_time_results, genus_recyc_fresh_results)
  write_csv(genus_time_results, file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_genus_baseMean5_all.csv"))
  write_csv(
    genus_time_results %>% filter(significant_padj_0.05),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_genus_baseMean5_significant_padj0.05.csv")
  )
  write_csv(
    genus_time_results %>% filter(significant_padj_0.05_lfc1),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_genus_baseMean5_significant_padj0.05_lfc1.csv")
  )
  write_csv(
    genus_time_results %>%
      filter(!is.na(padj)) %>%
      group_by(comparison_type, comparison, comparison_label) %>%
      slice_min(order_by = padj, n = 25, with_ties = FALSE) %>%
      ungroup(),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_genus_baseMean5_top25_by_comparison.csv")
  )
  write_csv(genus_recyc_fresh_results, file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_genus_baseMean5_all.csv"))
  write_csv(
    genus_recyc_fresh_results %>% filter(significant_padj_0.05),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_genus_baseMean5_significant_padj0.05.csv")
  )
  write_csv(
    genus_recyc_fresh_results %>% filter(significant_padj_0.05_lfc1),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_genus_baseMean5_significant_padj0.05_lfc1.csv")
  )
  write_csv(
    genus_recyc_fresh_results %>%
      filter(!is.na(padj)) %>%
      group_by(comparison_type, comparison, comparison_label) %>%
      slice_min(order_by = padj, n = 25, with_ties = FALSE) %>%
      ungroup(),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_genus_baseMean5_top25_by_comparison.csv")
  )
  all_results$Genus <- genus_results
}

if ("ASV" %in% levels_to_run) {
  message("Running aerobic timepoint ASV DESeq2 contrasts...")
  asv_time_results <- run_timepoint_contrasts(physeq, "ASV")
  asv_recyc_fresh_results <- run_recyc_vs_fresh_contrasts(physeq, "ASV")
  asv_results <- bind_rows(asv_time_results, asv_recyc_fresh_results)
  write_csv(asv_time_results, file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_asv_baseMean5_all.csv"))
  write_csv(
    asv_time_results %>% filter(significant_padj_0.05),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_asv_baseMean5_significant_padj0.05.csv")
  )
  write_csv(
    asv_time_results %>% filter(significant_padj_0.05_lfc1),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_asv_baseMean5_significant_padj0.05_lfc1.csv")
  )
  write_csv(
    asv_time_results %>%
      filter(!is.na(padj)) %>%
      group_by(comparison_type, comparison, comparison_label) %>%
      slice_min(order_by = padj, n = 25, with_ties = FALSE) %>%
      ungroup(),
    file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_asv_baseMean5_top25_by_comparison.csv")
  )
  write_csv(asv_recyc_fresh_results, file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_asv_baseMean5_all.csv"))
  write_csv(
    asv_recyc_fresh_results %>% filter(significant_padj_0.05),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_asv_baseMean5_significant_padj0.05.csv")
  )
  write_csv(
    asv_recyc_fresh_results %>% filter(significant_padj_0.05_lfc1),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_asv_baseMean5_significant_padj0.05_lfc1.csv")
  )
  write_csv(
    asv_recyc_fresh_results %>%
      filter(!is.na(padj)) %>%
      group_by(comparison_type, comparison, comparison_label) %>%
      slice_min(order_by = padj, n = 25, with_ties = FALSE) %>%
      ungroup(),
    file.path(DESEQ_DIR, "deseq2_recyc_vs_fresh_aerobic_timepoints_asv_baseMean5_top25_by_comparison.csv")
  )
  all_results$ASV <- asv_results
}

summary_table <- bind_rows(all_results) %>%
  group_by(feature_level, comparison_type, comparison, direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(feature_level, comparison_type, comparison, direction)
write_csv(summary_table, file.path(DESEQ_DIR, "deseq2_aerobic_timepoints_baseMean5_summary.csv"))

message("Done. Aerobic timepoint DESeq2 outputs written to: ", DESEQ_DIR)
