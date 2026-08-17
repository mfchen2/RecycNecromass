#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(vegan)
  library(ape)
  library(scales)
})

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
PICRUST_DIR <- Sys.getenv("PICRUST_DIR", file.path(PROJECT_DIR, "outputs/picrust2_ko_predictions"))
MEASURED_KO <- Sys.getenv(
  "MEASURED_KO",
  file.path(PROJECT_DIR, "outputs/manuscript_story_analysis/metaG_timepoint_KO_enrichment/timepoint_ko_sample_normalized_abundance.csv")
)
OUT_DIR <- file.path(PROJECT_DIR, "outputs/picrust2_vs_measured_ko")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(as.integer(Sys.getenv("PICRUST_COMPARE_SEED", "20260505")))
time_levels <- c("T2", "T3", "T5")

find_picrust_ko <- function(picrust_dir) {
  candidates <- c(
    file.path(picrust_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv.gz"),
    file.path(picrust_dir, "KO_metagenome_out", "pred_metagenome_unstrat.tsv"),
    file.path(picrust_dir, "pred_metagenome_unstrat.tsv.gz"),
    file.path(picrust_dir, "pred_metagenome_unstrat.tsv")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "Could not find PICRUSt2 KO prediction table. Expected one of:\n",
      paste(candidates, collapse = "\n")
    )
  }
  existing[[1]]
}

make_key <- function(source_type, source_carbon, oxygen, timepoint) {
  paste(source_type, source_carbon, oxygen, timepoint, sep = "__")
}

pcoa_scores <- function(dist_obj, max_axes = 5) {
  n <- attr(dist_obj, "Size")
  k <- max(1, min(max_axes, n - 1))
  fit <- cmdscale(dist_obj, k = k, eig = TRUE, add = TRUE)
  scores <- as.matrix(fit$points)
  keep <- seq_len(min(ncol(scores), k))
  scores[, keep, drop = FALSE]
}

read_picrust_predicted <- function() {
  ko_path <- find_picrust_ko(PICRUST_DIR)
  pred <- read_tsv(ko_path, show_col_types = FALSE, progress = FALSE)
  names(pred)[1] <- "ko"
  pred_long <- pred %>%
    pivot_longer(-ko, names_to = "sample_id", values_to = "predicted_abundance") %>%
    mutate(
      ko = str_replace(ko, "^ko:", ""),
      predicted_abundance = as.numeric(predicted_abundance)
    )

  meta <- read_tsv(file.path(PROJECT_DIR, "metadata/sample-metadata.tsv"), comment = "", show_col_types = FALSE) %>%
    rename(
      sample_id = `#SampleID`,
      source_type = necromass_type,
      source_carbon = necromass_source,
      oxygen = `oxic level`,
      timepoint = time_point
    ) %>%
    filter(sample_id != "#q2:types") %>%
    filter(
      source_type %in% c("Fresh", "Recyc"),
      source_carbon %in% c("Pseudomonas", "Arthrobacter"),
      oxygen %in% c("Aerobic", "Anaerobic"),
      timepoint %in% time_levels
    ) %>%
    mutate(condition_key = make_key(source_type, source_carbon, oxygen, timepoint))

  pred_long %>%
    inner_join(meta, by = "sample_id") %>%
    group_by(condition_key, source_type, source_carbon, oxygen, timepoint, ko) %>%
    summarise(predicted_abundance = mean(predicted_abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(condition_key) %>%
    mutate(predicted_fraction = predicted_abundance / sum(predicted_abundance, na.rm = TRUE)) %>%
    ungroup()
}

read_measured_ko <- function() {
  read_csv(MEASURED_KO, show_col_types = FALSE) %>%
    filter(
      source_type %in% c("Fresh", "Recyc"),
      source_carbon %in% c("Pseudomonas", "Arthrobacter"),
      oxygen %in% c("Aerobic", "Anaerobic"),
      timepoint %in% time_levels
    ) %>%
    mutate(condition_key = make_key(source_type, source_carbon, oxygen, timepoint)) %>%
    group_by(condition_key, source_type, source_carbon, oxygen, timepoint, ko, ko_definition) %>%
    summarise(measured_fraction = mean(fraction, na.rm = TRUE), .groups = "drop") %>%
    group_by(condition_key) %>%
    mutate(measured_fraction = measured_fraction / sum(measured_fraction, na.rm = TRUE)) %>%
    ungroup()
}

safe_mantel_procrustes <- function(pred_dist, meas_dist, meta, partition, level) {
  if (nrow(meta) < 4) {
    return(tibble(
      partition = partition, level = level, n_units = nrow(meta),
      mantel_r = NA_real_, mantel_p = NA_real_,
      protest_m2 = NA_real_, protest_r = NA_real_, protest_p = NA_real_,
      note = "Skipped: fewer than 4 units"
    ))
  }
  pred_sub <- as.dist(as.matrix(pred_dist)[meta$condition_key, meta$condition_key])
  meas_sub <- as.dist(as.matrix(meas_dist)[meta$condition_key, meta$condition_key])
  mantel_res <- mantel(pred_sub, meas_sub, method = "spearman", permutations = 999)
  ord_pred <- pcoa_scores(pred_sub)
  ord_meas <- pcoa_scores(meas_sub)
  k <- min(ncol(ord_pred), ncol(ord_meas))
  ord_pred <- ord_pred[, seq_len(k), drop = FALSE]
  ord_meas <- ord_meas[, seq_len(k), drop = FALSE]
  protest_res <- protest(ord_pred, ord_meas, permutations = 999, symmetric = TRUE)
  tibble(
    partition = partition,
    level = level,
    n_units = nrow(meta),
    mantel_r = unname(mantel_res$statistic),
    mantel_p = mantel_res$signif,
    protest_m2 = protest_res$ss,
    protest_r = protest_res$t0,
    protest_p = protest_res$signif,
    note = NA_character_
  )
}

pred_long <- read_picrust_predicted()
meas_long <- read_measured_ko()

common_kos <- intersect(unique(pred_long$ko), unique(meas_long$ko))
matched_keys <- intersect(unique(pred_long$condition_key), unique(meas_long$condition_key))

pred_mat <- pred_long %>%
  filter(ko %in% common_kos, condition_key %in% matched_keys) %>%
  select(condition_key, ko, predicted_fraction) %>%
  pivot_wider(names_from = ko, values_from = predicted_fraction, values_fill = 0) %>%
  arrange(condition_key) %>%
  column_to_rownames("condition_key") %>%
  as.matrix()

meas_mat <- meas_long %>%
  filter(ko %in% common_kos, condition_key %in% matched_keys) %>%
  select(condition_key, ko, measured_fraction) %>%
  pivot_wider(names_from = ko, values_from = measured_fraction, values_fill = 0) %>%
  arrange(condition_key) %>%
  column_to_rownames("condition_key") %>%
  as.matrix()

pred_mat <- pred_mat[sort(matched_keys), sort(common_kos), drop = FALSE]
meas_mat <- meas_mat[sort(matched_keys), sort(common_kos), drop = FALSE]

matched_meta <- pred_long %>%
  distinct(condition_key, source_type, source_carbon, oxygen, timepoint) %>%
  filter(condition_key %in% rownames(pred_mat)) %>%
  arrange(condition_key)

write_csv(matched_meta, file.path(OUT_DIR, "matched_condition_units.csv"))
write_csv(tibble(ko = common_kos), file.path(OUT_DIR, "common_kos_used.csv"))

pred_dist <- vegdist(pred_mat, method = "bray")
meas_dist <- vegdist(meas_mat, method = "bray")

partition_results <- bind_rows(
  safe_mantel_procrustes(pred_dist, meas_dist, matched_meta, "all_matched_conditions", "all"),
  map_dfr(sort(unique(matched_meta$source_type)), function(x) {
    safe_mantel_procrustes(pred_dist, meas_dist, matched_meta %>% filter(source_type == x), "within_source_type", x)
  }),
  map_dfr(sort(unique(matched_meta$source_carbon)), function(x) {
    safe_mantel_procrustes(pred_dist, meas_dist, matched_meta %>% filter(source_carbon == x), "within_source_carbon", x)
  }),
  map_dfr(sort(unique(matched_meta$oxygen)), function(x) {
    safe_mantel_procrustes(pred_dist, meas_dist, matched_meta %>% filter(oxygen == x), "within_oxygen", x)
  })
)
write_csv(partition_results, file.path(OUT_DIR, "picrust2_measured_ko_mantel_procrustes.csv"))

pair_ids <- combn(rownames(pred_mat), 2, simplify = FALSE)
pair_df <- map_dfr(pair_ids, function(pair) {
  tibble(
    condition_1 = pair[1],
    condition_2 = pair[2],
    predicted_bray = as.matrix(pred_dist)[pair[1], pair[2]],
    measured_bray = as.matrix(meas_dist)[pair[1], pair[2]]
  )
}) %>%
  left_join(matched_meta %>% rename_with(~paste0(.x, "_1"), -condition_key), by = c("condition_1" = "condition_key")) %>%
  left_join(matched_meta %>% rename_with(~paste0(.x, "_2"), -condition_key), by = c("condition_2" = "condition_key")) %>%
  mutate(
    pair_source_type = if_else(source_type_1 == source_type_2, source_type_1, "Fresh-Recyc"),
    pair_source_carbon = if_else(source_carbon_1 == source_carbon_2, source_carbon_1, "Mixed carbon"),
    pair_oxygen = if_else(oxygen_1 == oxygen_2, oxygen_1, "Mixed oxygen")
  )
write_csv(pair_df, file.path(OUT_DIR, "picrust2_measured_pairwise_bray.csv"))

ko_correlations <- map_dfr(colnames(pred_mat), function(ko_id) {
  pred <- pred_mat[, ko_id]
  meas <- meas_mat[, ko_id]
  tibble(
    ko = ko_id,
    n_units = length(pred),
    mean_predicted_fraction = mean(pred),
    mean_measured_fraction = mean(meas),
    spearman_rho = suppressWarnings(cor(pred, meas, method = "spearman")),
    pearson_r = suppressWarnings(cor(pred, meas, method = "pearson")),
    spearman_p = suppressWarnings(cor.test(pred, meas, method = "spearman", exact = FALSE)$p.value),
    pearson_p = suppressWarnings(cor.test(pred, meas, method = "pearson")$p.value)
  )
}) %>%
  left_join(meas_long %>% distinct(ko, ko_definition), by = "ko") %>%
  mutate(
    spearman_q = p.adjust(spearman_p, method = "BH"),
    pearson_q = p.adjust(pearson_p, method = "BH"),
    abs_spearman_rho = abs(spearman_rho)
  ) %>%
  arrange(desc(abs_spearman_rho))
write_csv(ko_correlations, file.path(OUT_DIR, "ko_predictability_correlations.csv"))
write_csv(
  ko_correlations %>% filter(!is.na(spearman_q), spearman_q < 0.05) %>% arrange(desc(spearman_rho)),
  file.path(OUT_DIR, "ko_well_predicted_spearman_q0.05.csv")
)
write_csv(
  ko_correlations %>% filter(is.na(spearman_q) | spearman_q >= 0.05) %>% arrange(desc(mean_measured_fraction)),
  file.path(OUT_DIR, "ko_poorly_predicted_or_uncertain.csv")
)

pair_plot <- ggplot(pair_df, aes(predicted_bray, measured_bray, color = pair_source_type, shape = pair_oxygen)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Fresh" = "#2f7fc1", "Recyc" = "#d7433f", "Fresh-Recyc" = "#666666")) +
  labs(x = "PICRUSt2-predicted KO Bray-Curtis", y = "Measured KO Bray-Curtis", color = "Pair source", shape = "Pair oxygen") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_DIR, "picrust2_vs_measured_ko_pairwise_bray.png"), pair_plot, width = 7.5, height = 5.8, dpi = 300)
ggsave(file.path(OUT_DIR, "picrust2_vs_measured_ko_pairwise_bray.pdf"), pair_plot, width = 7.5, height = 5.8)

top_plot_df <- ko_correlations %>%
  filter(is.finite(spearman_rho)) %>%
  slice_max(abs_spearman_rho, n = 30) %>%
  mutate(
    ko_label = if_else(is.na(ko_definition), ko, paste0(ko, " | ", str_trunc(ko_definition, 42))),
    ko_label = reorder(ko_label, spearman_rho)
  )
ko_plot <- ggplot(top_plot_df, aes(spearman_rho, ko_label, fill = spearman_q < 0.05)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, color = "grey45") +
  scale_fill_manual(values = c("TRUE" = "#3f7f3f", "FALSE" = "#b8b8b8"), labels = c("FALSE" = "q >= 0.05", "TRUE" = "q < 0.05")) +
  labs(x = "Spearman correlation: predicted versus measured KO", y = NULL, fill = NULL) +
  theme_bw(base_size = 9)
ggsave(file.path(OUT_DIR, "top_ko_predictability_correlations.png"), ko_plot, width = 8.5, height = 8, dpi = 300)
ggsave(file.path(OUT_DIR, "top_ko_predictability_correlations.pdf"), ko_plot, width = 8.5, height = 8)

message("Matched condition units: ", nrow(matched_meta))
message("Common KOs: ", length(common_kos))
message("Outputs written to: ", OUT_DIR)
