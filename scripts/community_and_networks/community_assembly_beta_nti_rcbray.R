#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(qiime2R)
  library(phyloseq)
  library(tidyverse)
  library(ape)
  library(vegan)
  library(matrixStats)
  library(scales)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
QIIME_DIR <- file.path(PROJECT_DIR, "qiime2-results-local")
OUT_DIR <- file.path(PROJECT_DIR, "outputs/16s_community_assembly")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_PERM_BNTI <- as.integer(Sys.getenv("BNTI_N_PERM", "199"))
N_PERM_RCBRAY <- as.integer(Sys.getenv("RCBRAY_N_PERM", "999"))
MIN_PREVALENCE <- as.integer(Sys.getenv("ASSEMBLY_MIN_PREVALENCE", "3"))
MIN_TOTAL_COUNT <- as.integer(Sys.getenv("ASSEMBLY_MIN_TOTAL_COUNT", "10"))
set.seed(as.integer(Sys.getenv("ASSEMBLY_SEED", "20260505")))

theme_assembly <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey94", color = "grey70"),
      legend.position = "right"
    )
}

clean_metadata <- function(ps) {
  meta <- data.frame(as(sample_data(ps), "data.frame"), check.names = FALSE)
  if ("oxic.level" %in% names(meta) && !"oxic level" %in% names(meta)) {
    names(meta)[names(meta) == "oxic.level"] <- "oxic level"
  }
  meta$sample_id <- rownames(meta)
  meta %>%
    mutate(
      source_type = as.character(necromass_type),
      source_carbon = as.character(necromass_source),
      oxygen = as.character(`oxic level`),
      time_point = as.character(time_point),
      replicate = as.character(replicate)
    )
}

weighted_beta_mntd <- function(rel, phylo_dist) {
  rel <- as.matrix(rel)
  sample_taxa <- lapply(seq_len(nrow(rel)), function(i) which(rel[i, ] > 0))
  min_to_sample <- matrix(0, nrow = ncol(rel), ncol = nrow(rel))
  for (j in seq_along(sample_taxa)) {
    taxa_j <- sample_taxa[[j]]
    min_to_sample[, j] <- rowMins(phylo_dist[, taxa_j, drop = FALSE])
  }
  directional <- rel %*% min_to_sample
  (directional + t(directional)) / 2
}

pair_table_from_matrix <- function(mat, meta, value_name) {
  pairs <- combn(seq_len(nrow(mat)), 2)
  tibble(
    sample_1 = rownames(mat)[pairs[1, ]],
    sample_2 = rownames(mat)[pairs[2, ]],
    !!value_name := mat[pairs[1, ] + (pairs[2, ] - 1) * nrow(mat)]
  ) %>%
    left_join(meta %>% rename_with(~ paste0(.x, "_1"), -sample_id), by = c("sample_1" = "sample_id")) %>%
    left_join(meta %>% rename_with(~ paste0(.x, "_2"), -sample_id), by = c("sample_2" = "sample_id")) %>%
    mutate(
      pair_source_type = case_when(
        source_type_1 == source_type_2 ~ source_type_1,
        TRUE ~ "Fresh-Recyc"
      ),
      pair_source_carbon = case_when(
        source_carbon_1 == source_carbon_2 ~ source_carbon_1,
        TRUE ~ "Mixed carbon"
      ),
      pair_oxygen = case_when(
        oxygen_1 == oxygen_2 ~ oxygen_1,
        TRUE ~ "Mixed oxygen"
      ),
      pair_time = case_when(
        time_point_1 == time_point_2 ~ time_point_1,
        TRUE ~ "Mixed time"
      ),
      same_condition = source_type_1 == source_type_2 &
        source_carbon_1 == source_carbon_2 &
        oxygen_1 == oxygen_2 &
        time_point_1 == time_point_2,
      same_regime = source_type_1 == source_type_2
    )
}

classify_process <- function(beta_nti, rc_bray) {
  case_when(
    beta_nti > 2 ~ "Variable selection",
    beta_nti < -2 ~ "Homogeneous selection",
    abs(beta_nti) <= 2 & rc_bray > 0.95 ~ "Dispersal limitation",
    abs(beta_nti) <= 2 & rc_bray < -0.95 ~ "Homogenizing dispersal",
    abs(beta_nti) <= 2 ~ "Drift/undominated",
    TRUE ~ NA_character_
  )
}

pairwise_path <- file.path(OUT_DIR, "beta_nti_rcbray_pairwise.csv")
metadata_path <- file.path(OUT_DIR, "assembly_sample_metadata_used.csv")
using_existing_pairwise <- file.exists(pairwise_path) && !tolower(Sys.getenv("ASSEMBLY_FORCE_RERUN", "false")) %in% c("true", "1", "yes")
if (using_existing_pairwise) {
  message("Using existing pairwise table: ", pairwise_path)
  pair_df <- read_csv(pairwise_path, show_col_types = FALSE) %>%
    mutate(pair_source_type = factor(pair_source_type, levels = c("Fresh", "Recyc", "Fresh-Recyc")))
} else {
  ps <- qza_to_phyloseq(
    features = file.path(QIIME_DIR, "artifacts/table.qza"),
    tree = file.path(QIIME_DIR, "artifacts/rooted-tree.qza"),
    metadata = file.path(PROJECT_DIR, "metadata/sample-metadata.tsv")
  )

  meta_all <- clean_metadata(ps)
  rownames(meta_all) <- meta_all$sample_id
  sample_data(ps) <- sample_data(meta_all %>% select(-sample_id))

  ps <- subset_samples(
    ps,
    source_type %in% c("Fresh", "Recyc") &
      source_carbon %in% c("Pseudomonas", "Arthrobacter") &
      oxygen %in% c("Aerobic", "Anaerobic") &
      time_point %in% paste0("T", 1:5)
  )
  ps <- prune_samples(sample_sums(ps) > 0, ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)

  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)

  taxa_keep <- colSums(otu) >= MIN_TOTAL_COUNT & colSums(otu > 0) >= MIN_PREVALENCE
  ps <- prune_taxa(colnames(otu)[taxa_keep], ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- otu[rowSums(otu) > 0, colSums(otu) > 0, drop = FALSE]

  meta <- data.frame(as(sample_data(ps), "data.frame"), check.names = FALSE) %>%
    rownames_to_column("sample_id") %>%
    mutate(
      source_type = factor(source_type, levels = c("Fresh", "Recyc")),
      source_carbon = factor(source_carbon, levels = c("Arthrobacter", "Pseudomonas")),
      oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic")),
      time_point = factor(time_point, levels = paste0("T", 1:5))
    )
  otu <- otu[meta$sample_id, , drop = FALSE]
  rel <- sweep(otu, 1, rowSums(otu), "/")
  rel[is.na(rel)] <- 0

  tree <- phy_tree(ps)
  tree <- keep.tip(tree, colnames(otu))
  phylo_dist <- cophenetic.phylo(tree)
  phylo_dist <- phylo_dist[colnames(otu), colnames(otu)]

  message("Samples: ", nrow(otu), "; ASVs retained: ", ncol(otu))
  message("Computing observed betaMNTD...")
  obs_bmntd <- weighted_beta_mntd(rel, phylo_dist)
  rownames(obs_bmntd) <- colnames(obs_bmntd) <- rownames(rel)

  message("Running beta-NTI null model with ", N_PERM_BNTI, " taxa-shuffle permutations...")
  null_sum <- matrix(0, nrow = nrow(obs_bmntd), ncol = ncol(obs_bmntd))
  null_sumsq <- matrix(0, nrow = nrow(obs_bmntd), ncol = ncol(obs_bmntd))
  for (i in seq_len(N_PERM_BNTI)) {
    perm <- sample(seq_len(ncol(phylo_dist)))
    dist_perm <- phylo_dist[perm, perm]
    dimnames(dist_perm) <- dimnames(phylo_dist)
    null_i <- weighted_beta_mntd(rel, dist_perm)
    null_sum <- null_sum + null_i
    null_sumsq <- null_sumsq + null_i^2
    if (i %% 25 == 0 || i == N_PERM_BNTI) message("  beta-NTI permutation ", i, "/", N_PERM_BNTI)
  }
  null_mean <- null_sum / N_PERM_BNTI
  null_var <- pmax((null_sumsq - (null_sum^2 / N_PERM_BNTI)) / (N_PERM_BNTI - 1), 0)
  null_sd <- sqrt(null_var)
  beta_nti <- (obs_bmntd - null_mean) / null_sd
  beta_nti[!is.finite(beta_nti)] <- NA_real_
  diag(beta_nti) <- NA_real_

  message("Running RCbray null model with ", N_PERM_RCBRAY, " richness/occurrence-weighted permutations...")
  obs_bray <- as.matrix(vegdist(rel, method = "bray"))
  taxa_occurrence <- colSums(otu > 0)
  taxa_weights <- taxa_occurrence / sum(taxa_occurrence)
  sample_richness <- rowSums(otu > 0)
  less_count <- matrix(0, nrow = nrow(otu), ncol = nrow(otu), dimnames = list(rownames(otu), rownames(otu)))
  tie_count <- less_count

  for (i in seq_len(N_PERM_RCBRAY)) {
    null_counts <- matrix(0, nrow = nrow(otu), ncol = ncol(otu), dimnames = dimnames(otu))
    for (s in seq_len(nrow(otu))) {
      present <- which(otu[s, ] > 0)
      selected <- sample(seq_len(ncol(otu)), size = length(present), replace = FALSE, prob = taxa_weights)
      null_counts[s, selected] <- sample(otu[s, present], size = length(present), replace = FALSE)
    }
    null_rel <- sweep(null_counts, 1, rowSums(null_counts), "/")
    null_rel[is.na(null_rel)] <- 0
    null_bray <- as.matrix(vegdist(null_rel, method = "bray"))
    less_count <- less_count + (null_bray < obs_bray)
    tie_count <- tie_count + (abs(null_bray - obs_bray) < .Machine$double.eps^0.5)
    if (i %% 100 == 0 || i == N_PERM_RCBRAY) message("  RCbray permutation ", i, "/", N_PERM_RCBRAY)
  }
  rc_bray <- 2 * ((less_count + 0.5 * tie_count) / N_PERM_RCBRAY) - 1
  diag(rc_bray) <- NA_real_

  pair_df <- pair_table_from_matrix(beta_nti, meta, "beta_nti") %>%
    select(-starts_with("source_type_"), -starts_with("source_carbon_"), -starts_with("oxygen_"), -starts_with("time_point_"), -starts_with("replicate_")) %>%
    left_join(
      pair_table_from_matrix(rc_bray, meta, "rc_bray") %>%
        select(sample_1, sample_2, rc_bray),
      by = c("sample_1", "sample_2")
    ) %>%
    left_join(
      pair_table_from_matrix(obs_bray, meta, "bray_curtis") %>%
        select(sample_1, sample_2, bray_curtis),
      by = c("sample_1", "sample_2")
    ) %>%
    mutate(
      assembly_process = classify_process(beta_nti, rc_bray),
      pair_source_type = factor(pair_source_type, levels = c("Fresh", "Recyc", "Fresh-Recyc"))
    )

  write_csv(pair_df, pairwise_path)

  metadata_used <- meta %>%
    mutate(total_counts_retained = rowSums(otu), richness_retained = sample_richness)
  write_csv(metadata_used, metadata_path)
}

process_summary <- pair_df %>%
  filter(!is.na(assembly_process)) %>%
  dplyr::count(pair_source_type, assembly_process, name = "n_pairs") %>%
  group_by(pair_source_type) %>%
  mutate(frac_pairs = n_pairs / sum(n_pairs)) %>%
  ungroup()
write_csv(process_summary, file.path(OUT_DIR, "assembly_process_summary_by_source_type_pair.csv"))

process_summary_same_condition <- pair_df %>%
  filter(same_condition, !is.na(assembly_process)) %>%
  dplyr::count(pair_source_type, pair_source_carbon, pair_oxygen, pair_time, assembly_process, name = "n_pairs") %>%
  group_by(pair_source_type, pair_source_carbon, pair_oxygen, pair_time) %>%
  mutate(frac_pairs = n_pairs / sum(n_pairs)) %>%
  ungroup()
write_csv(process_summary_same_condition, file.path(OUT_DIR, "assembly_process_summary_same_condition_replicates.csv"))

source_tests <- bind_rows(
  pair_df %>%
    filter(pair_source_type %in% c("Fresh", "Recyc")) %>%
    summarise(
      metric = "beta_nti",
      fresh_median = median(beta_nti[pair_source_type == "Fresh"], na.rm = TRUE),
      recyc_median = median(beta_nti[pair_source_type == "Recyc"], na.rm = TRUE),
      wilcox_p = wilcox.test(beta_nti ~ pair_source_type)$p.value
    ),
  pair_df %>%
    filter(pair_source_type %in% c("Fresh", "Recyc"), abs(beta_nti) <= 2) %>%
    summarise(
      metric = "rc_bray_stochastic_pairs",
      fresh_median = median(rc_bray[pair_source_type == "Fresh"], na.rm = TRUE),
      recyc_median = median(rc_bray[pair_source_type == "Recyc"], na.rm = TRUE),
      wilcox_p = wilcox.test(rc_bray ~ pair_source_type)$p.value
    )
) %>%
  mutate(p_adj_bh = p.adjust(wilcox_p, method = "BH"))
write_csv(source_tests, file.path(OUT_DIR, "fresh_vs_recyc_beta_nti_rcbray_tests.csv"))

summary_stats <- pair_df %>%
  group_by(pair_source_type) %>%
  summarise(
    n_pairs = n(),
    median_beta_nti = median(beta_nti, na.rm = TRUE),
    mean_beta_nti = mean(beta_nti, na.rm = TRUE),
    frac_abs_beta_nti_lt_2 = mean(abs(beta_nti) <= 2, na.rm = TRUE),
    frac_homogeneous_selection = mean(assembly_process == "Homogeneous selection", na.rm = TRUE),
    frac_variable_selection = mean(assembly_process == "Variable selection", na.rm = TRUE),
    frac_stochastic_undominated = mean(assembly_process == "Drift/undominated", na.rm = TRUE),
    frac_dispersal_limitation = mean(assembly_process == "Dispersal limitation", na.rm = TRUE),
    frac_homogenizing_dispersal = mean(assembly_process == "Homogenizing dispersal", na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_stats, file.path(OUT_DIR, "assembly_summary_stats_by_source_type_pair.csv"))

run_info <- tibble(
  parameter = c("n_samples", "n_asvs_retained", "min_prevalence", "min_total_count", "beta_nti_permutations", "rcbray_permutations", "seed", "resumed_from_pairwise_table"),
  value = c(
    length(unique(c(pair_df$sample_1, pair_df$sample_2))),
    if (exists("otu")) ncol(otu) else NA_integer_,
    MIN_PREVALENCE,
    MIN_TOTAL_COUNT,
    N_PERM_BNTI,
    N_PERM_RCBRAY,
    Sys.getenv("ASSEMBLY_SEED", "20260505"),
    using_existing_pairwise
  )
)
write_csv(run_info, file.path(OUT_DIR, "assembly_run_info.csv"))

process_plot <- ggplot(process_summary, aes(pair_source_type, frac_pairs, fill = assembly_process)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = "Pairwise comparisons", fill = "Assembly process") +
  theme_assembly()
ggsave(file.path(OUT_DIR, "assembly_process_stacked_bar_by_source_type.pdf"), process_plot, width = 7, height = 5)
ggsave(file.path(OUT_DIR, "assembly_process_stacked_bar_by_source_type.png"), process_plot, width = 7, height = 5, dpi = 300)

bnti_plot <- ggplot(pair_df %>% filter(pair_source_type %in% c("Fresh", "Recyc")), aes(pair_source_type, beta_nti, fill = pair_source_type)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey40") +
  geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.15) +
  scale_fill_manual(values = c("Fresh" = "#2f7fc1", "Recyc" = "#d7433f")) +
  labs(x = NULL, y = expression(beta*"-NTI")) +
  theme_assembly() +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "beta_nti_fresh_vs_recyc_violin.pdf"), bnti_plot, width = 6, height = 5)
ggsave(file.path(OUT_DIR, "beta_nti_fresh_vs_recyc_violin.png"), bnti_plot, width = 6, height = 5, dpi = 300)

rc_plot <- ggplot(pair_df %>% filter(pair_source_type %in% c("Fresh", "Recyc"), abs(beta_nti) <= 2), aes(pair_source_type, rc_bray, fill = pair_source_type)) +
  geom_hline(yintercept = c(-0.95, 0.95), linetype = "dashed", color = "grey40") +
  geom_violin(trim = TRUE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.18, outlier.alpha = 0.15) +
  scale_fill_manual(values = c("Fresh" = "#2f7fc1", "Recyc" = "#d7433f")) +
  labs(x = NULL, y = "RCbray among |beta-NTI| <= 2 pairs") +
  theme_assembly() +
  theme(legend.position = "none")
ggsave(file.path(OUT_DIR, "rcbray_stochastic_pairs_fresh_vs_recyc_violin.pdf"), rc_plot, width = 6, height = 5)
ggsave(file.path(OUT_DIR, "rcbray_stochastic_pairs_fresh_vs_recyc_violin.png"), rc_plot, width = 6, height = 5, dpi = 300)

message("Done. Outputs written to: ", OUT_DIR)
