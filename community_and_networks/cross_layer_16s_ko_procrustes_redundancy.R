suppressPackageStartupMessages({
  library(biomformat)
  library(vegan)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(tibble)
})

out_dir <- "outputs/manuscript_story_analysis/cross_layer_16s_ko_redundancy"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260505)

time_levels <- c("T2", "T3", "T5")

make_key <- function(source_type, source_carbon, oxygen, timepoint) {
  paste(source_type, source_carbon, oxygen, timepoint, sep = "__")
}

bray_one <- function(x, y) {
  denom <- sum(x + y, na.rm = TRUE)
  if (denom == 0) return(NA_real_)
  sum(abs(x - y), na.rm = TRUE) / denom
}

read_16s_condition_centroids <- function() {
  biom <- read_biom("qiime2-results-local/exports/feature-table/feature-table.biom")
  mat <- as.matrix(biom_data(biom))
  mat <- t(mat)

  meta <- read_tsv("metadata/sample-metadata.tsv", comment = "", show_col_types = FALSE) %>%
    rename(
      sample_id = `#SampleID`,
      oxygen = `oxic level`,
      source_type = necromass_type,
      source_carbon = necromass_source,
      timepoint = time_point
    ) %>%
    filter(sample_id != "#q2:types") %>%
    filter(
      source_type %in% c("Fresh", "Recyc"),
      source_carbon %in% c("Pseudomonas", "Arthrobacter"),
      oxygen %in% c("Aerobic", "Anaerobic"),
      timepoint %in% time_levels
    ) %>%
    mutate(
      timepoint = factor(timepoint, levels = time_levels),
      condition_key = make_key(source_type, source_carbon, oxygen, timepoint)
    )

  mat <- mat[meta$sample_id, , drop = FALSE]
  mat <- sweep(mat, 1, rowSums(mat), "/")
  mat[is.na(mat)] <- 0

  centroid_meta <- meta %>%
    distinct(condition_key, source_type, source_carbon, oxygen, timepoint) %>%
    arrange(source_type, source_carbon, oxygen, timepoint)

  centroid_mat <- map_dfr(centroid_meta$condition_key, function(key) {
    ids <- meta$sample_id[meta$condition_key == key]
    as_tibble(t(colMeans(mat[ids, , drop = FALSE])), .name_repair = "minimal")
  }) %>%
    as.matrix()

  colnames(centroid_mat) <- colnames(mat)
  rownames(centroid_mat) <- centroid_meta$condition_key

  list(mat = centroid_mat, meta = centroid_meta)
}

read_ko_profiles <- function() {
  ko_long <- read_csv(
    "outputs/manuscript_story_analysis/metaG_timepoint_KO_enrichment/timepoint_ko_sample_normalized_abundance.csv",
    show_col_types = FALSE
  ) %>%
    filter(
      source_type %in% c("Fresh", "Recyc"),
      source_carbon %in% c("Pseudomonas", "Arthrobacter"),
      oxygen %in% c("Aerobic", "Anaerobic"),
      timepoint %in% time_levels
    ) %>%
    mutate(
      timepoint = factor(timepoint, levels = time_levels),
      condition_key = make_key(source_type, source_carbon, oxygen, timepoint)
    )

  meta <- ko_long %>%
    distinct(sample_id, condition_key, source_type, source_carbon, oxygen, timepoint, total_ko_gene_hits) %>%
    arrange(source_type, source_carbon, oxygen, timepoint)

  mat <- ko_long %>%
    select(condition_key, ko, fraction) %>%
    distinct() %>%
    pivot_wider(names_from = ko, values_from = fraction, values_fill = 0) %>%
    arrange(match(condition_key, meta$condition_key)) %>%
    column_to_rownames("condition_key") %>%
    as.matrix()

  list(mat = mat, meta = meta)
}

pcoa_scores <- function(dist_obj, max_axes = 5) {
  n <- attr(dist_obj, "Size")
  k <- max(1, min(max_axes, n - 1))
  fit <- cmdscale(dist_obj, k = k, eig = TRUE, add = TRUE)
  scores <- as.matrix(fit$points)
  pos_axes <- which(fit$eig > 1e-9)
  keep <- seq_len(min(ncol(scores), max(2, length(pos_axes)), k))
  scores[, keep, drop = FALSE]
}

safe_cross_layer <- function(dist16, distko, meta, partition, level) {
  if (nrow(meta) < 4) {
    return(tibble(
      partition = partition, level = level, n_units = nrow(meta),
      mantel_r = NA_real_, mantel_p = NA_real_,
      protest_m2 = NA_real_, protest_r = NA_real_, protest_p = NA_real_,
      note = "Skipped: fewer than 4 matched units"
    ))
  }

  d16 <- as.dist(as.matrix(dist16)[meta$condition_key, meta$condition_key])
  dko <- as.dist(as.matrix(distko)[meta$condition_key, meta$condition_key])

  mantel_res <- mantel(d16, dko, method = "spearman", permutations = 999)

  ord16 <- pcoa_scores(d16)
  ordko <- pcoa_scores(dko)
  k <- min(ncol(ord16), ncol(ordko))
  ord16 <- ord16[, seq_len(k), drop = FALSE]
  ordko <- ordko[, seq_len(k), drop = FALSE]
  rownames(ord16) <- rownames(ordko) <- meta$condition_key

  protest_res <- protest(ord16, ordko, permutations = 999, symmetric = TRUE)

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

run_partitions <- function(dist16, distko, meta) {
  bind_rows(
    safe_cross_layer(dist16, distko, meta, "all_matched_conditions", "all"),
    map_dfr(sort(unique(meta$source_type)), function(x) {
      safe_cross_layer(dist16, distko, meta %>% filter(source_type == x), "within_necromass_origin", x)
    }),
    map_dfr(sort(unique(meta$source_carbon)), function(x) {
      safe_cross_layer(dist16, distko, meta %>% filter(source_carbon == x), "within_necromass_source_organism", x)
    }),
    map_dfr(sort(unique(meta$oxygen)), function(x) {
      safe_cross_layer(dist16, distko, meta %>% filter(oxygen == x), "within_oxygen", x)
    }),
    meta %>%
      distinct(source_type, oxygen) %>%
      pmap_dfr(function(source_type, oxygen) {
        src_type <- source_type
        oxy <- oxygen
        safe_cross_layer(
          dist16, distko,
          meta %>% filter(.data$source_type == src_type, .data$oxygen == oxy),
          "within_necromass_origin_x_oxygen",
          paste(src_type, oxy, sep = " | ")
        )
      }),
    meta %>%
      distinct(source_carbon, oxygen) %>%
      pmap_dfr(function(source_carbon, oxygen) {
        src_carbon <- source_carbon
        oxy <- oxygen
        safe_cross_layer(
          dist16, distko,
          meta %>% filter(.data$source_carbon == src_carbon, .data$oxygen == oxy),
          "within_source_organism_x_oxygen",
          paste(src_carbon, oxy, sep = " | ")
        )
      })
  )
}

paired_distance_table <- function(dist16, distko, meta, cell_cols, cell_label) {
  d16m <- as.matrix(dist16)
  dkom <- as.matrix(distko)

  meta %>%
    group_by(across(all_of(cell_cols))) %>%
    group_modify(~{
      ids <- .x$condition_key
      if (length(ids) < 3) return(tibble())
      pairs <- combn(ids, 2, simplify = FALSE)
      map_dfr(pairs, function(pair) {
        tibble(
          unit_1 = pair[1],
          unit_2 = pair[2],
          beta_16s = d16m[pair[1], pair[2]],
          beta_ko = dkom[pair[1], pair[2]],
          beta_diff_ko_minus_16s = beta_ko - beta_16s
        )
      })
    }) %>%
    ungroup() %>%
    mutate(cell_definition = cell_label, .before = 1)
}

centroid_dispersion_table <- function(dist16, distko, meta, cell_cols, cell_label) {
  calc_layer <- function(dist_obj, layer_name) {
    pcoa <- pcoa_scores(dist_obj, max_axes = 10)
    pcoa <- as_tibble(pcoa, rownames = "condition_key")

    pcoa %>%
      left_join(meta, by = "condition_key") %>%
      group_by(across(all_of(cell_cols))) %>%
      group_modify(~{
        coord_cols <- setdiff(names(.x), names(meta))
        coord_cols <- setdiff(coord_cols, cell_cols)
        coords <- as.matrix(.x[, coord_cols, drop = FALSE])
        centroid <- colMeans(coords)
        tibble(
          condition_key = .x$condition_key,
          distance_to_cell_centroid = sqrt(rowSums(sweep(coords, 2, centroid)^2))
        )
      }) %>%
      ungroup() %>%
      mutate(layer = layer_name)
  }

  bind_rows(calc_layer(dist16, "16S"), calc_layer(distko, "KO")) %>%
    mutate(cell_definition = cell_label, .before = 1)
}

summarise_pairwise_tests <- function(pair_df, cell_cols) {
  pair_df %>%
    group_by(cell_definition, across(all_of(cell_cols))) %>%
    summarise(
      n_units = n_distinct(c(unit_1, unit_2)),
      n_pairs = n(),
      mean_beta_16s = mean(beta_16s),
      mean_beta_ko = mean(beta_ko),
      median_beta_16s = median(beta_16s),
      median_beta_ko = median(beta_ko),
      mean_diff_ko_minus_16s = mean(beta_diff_ko_minus_16s),
      median_diff_ko_minus_16s = median(beta_diff_ko_minus_16s),
      wilcox_p = if_else(
        n_pairs >= 3 && sd(beta_diff_ko_minus_16s) > 0,
        wilcox.test(beta_ko, beta_16s, paired = TRUE, exact = FALSE, alternative = "less")$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    group_by(cell_definition) %>%
    mutate(wilcox_q = p.adjust(wilcox_p, method = "BH")) %>%
    ungroup() %>%
    mutate(
      ko_lower_than_16s = mean_diff_ko_minus_16s < 0,
      significant_lower_at_0.05 = ko_lower_than_16s & wilcox_p < 0.05
    )
}

summarise_dispersion_tests <- function(disp_df, cell_cols) {
  wide <- disp_df %>%
    select(cell_definition, all_of(cell_cols), condition_key, layer, distance_to_cell_centroid) %>%
    pivot_wider(names_from = layer, values_from = distance_to_cell_centroid) %>%
    rename(disp_16s = `16S`, disp_ko = KO) %>%
    mutate(disp_diff_ko_minus_16s = disp_ko - disp_16s)

  wide %>%
    group_by(cell_definition, across(all_of(cell_cols))) %>%
    summarise(
      n_units = n(),
      mean_disp_16s = mean(disp_16s),
      mean_disp_ko = mean(disp_ko),
      median_disp_16s = median(disp_16s),
      median_disp_ko = median(disp_ko),
      mean_diff_ko_minus_16s = mean(disp_diff_ko_minus_16s),
      wilcox_p = if_else(
        n_units >= 3 && sd(disp_diff_ko_minus_16s) > 0,
        wilcox.test(disp_ko, disp_16s, paired = TRUE, exact = FALSE, alternative = "less")$p.value,
        NA_real_
      ),
      .groups = "drop"
    ) %>%
    group_by(cell_definition) %>%
    mutate(wilcox_q = p.adjust(wilcox_p, method = "BH")) %>%
    ungroup() %>%
    mutate(
      ko_lower_than_16s = mean_diff_ko_minus_16s < 0,
      significant_lower_at_0.05 = ko_lower_than_16s & wilcox_p < 0.05
    )
}

dat16 <- read_16s_condition_centroids()
datko <- read_ko_profiles()

matched_keys <- intersect(rownames(dat16$mat), rownames(datko$mat))

matched_meta <- datko$meta %>%
  filter(condition_key %in% matched_keys) %>%
  arrange(source_type, source_carbon, oxygen, timepoint) %>%
  mutate(
    condition_key = as.character(condition_key),
    timepoint = factor(timepoint, levels = time_levels)
  )

mat16 <- dat16$mat[matched_meta$condition_key, , drop = FALSE]
matko <- datko$mat[matched_meta$condition_key, , drop = FALSE]

dist16 <- vegdist(mat16, method = "bray")
distko <- vegdist(matko, method = "bray")

matched_meta_out <- matched_meta %>%
  select(condition_key, sample_id, source_type, source_carbon, oxygen, timepoint, total_ko_gene_hits)

write_csv(matched_meta_out, file.path(out_dir, "matched_condition_time_units.csv"))

cross_layer_results <- run_partitions(dist16, distko, matched_meta)
write_csv(cross_layer_results, file.path(out_dir, "procrustes_mantel_partitioned_results.csv"))

pair_origin_oxygen <- paired_distance_table(
  dist16, distko, matched_meta,
  c("source_type", "oxygen"),
  "necromass_origin_x_oxygen"
)
pair_source_oxygen <- paired_distance_table(
  dist16, distko, matched_meta,
  c("source_carbon", "oxygen"),
  "source_organism_x_oxygen"
)
pair_all <- bind_rows(pair_origin_oxygen, pair_source_oxygen)
write_csv(pair_all, file.path(out_dir, "paired_beta_distances_within_cells.csv"))

pair_tests <- bind_rows(
  summarise_pairwise_tests(pair_origin_oxygen, c("source_type", "oxygen")),
  summarise_pairwise_tests(pair_source_oxygen, c("source_carbon", "oxygen"))
)
write_csv(pair_tests, file.path(out_dir, "functional_redundancy_pairwise_beta_tests.csv"))

disp_origin_oxygen <- centroid_dispersion_table(
  dist16, distko, matched_meta,
  c("source_type", "oxygen"),
  "necromass_origin_x_oxygen"
)
disp_source_oxygen <- centroid_dispersion_table(
  dist16, distko, matched_meta,
  c("source_carbon", "oxygen"),
  "source_organism_x_oxygen"
)
disp_all <- bind_rows(disp_origin_oxygen, disp_source_oxygen)
write_csv(disp_all, file.path(out_dir, "centroid_dispersion_within_cells_long.csv"))

disp_tests <- bind_rows(
  summarise_dispersion_tests(disp_origin_oxygen, c("source_type", "oxygen")),
  summarise_dispersion_tests(disp_source_oxygen, c("source_carbon", "oxygen"))
)
write_csv(disp_tests, file.path(out_dir, "functional_redundancy_centroid_dispersion_tests.csv"))

ord16 <- pcoa_scores(dist16, max_axes = 4) %>%
  as_tibble(rownames = "condition_key") %>%
  rename_with(~paste0("PCoA16S_", seq_along(.x)), -condition_key)
ordko <- pcoa_scores(distko, max_axes = 4) %>%
  as_tibble(rownames = "condition_key") %>%
  rename_with(~paste0("PCoAKO_", seq_along(.x)), -condition_key)

ord_out <- matched_meta %>%
  select(condition_key, source_type, source_carbon, oxygen, timepoint) %>%
  left_join(ord16, by = "condition_key") %>%
  left_join(ordko, by = "condition_key")
write_csv(ord_out, file.path(out_dir, "matched_16s_ko_pcoa_scores.csv"))

p_scatter <- pair_all %>%
  mutate(
    cell = case_when(
      cell_definition == "necromass_origin_x_oxygen" ~ paste(source_type, oxygen, sep = " | "),
      TRUE ~ paste(source_carbon, oxygen, sep = " | ")
    )
  ) %>%
  ggplot(aes(x = beta_16s, y = beta_ko)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey45") +
  geom_point(alpha = 0.75, size = 1.9) +
  facet_wrap(~cell_definition + cell, scales = "free") +
  theme_bw(base_size = 9) +
  labs(
    x = "16S Bray-Curtis distance",
    y = "KO Bray-Curtis distance",
    title = "Matched condition-pair beta diversity: 16S versus KO"
  )

ggsave(file.path(out_dir, "paired_beta_16s_vs_ko_scatter.png"), p_scatter, width = 10.5, height = 7, dpi = 300)
ggsave(file.path(out_dir, "paired_beta_16s_vs_ko_scatter.pdf"), p_scatter, width = 10.5, height = 7)

p_summary <- pair_tests %>%
  mutate(
    cell = if_else(
      cell_definition == "necromass_origin_x_oxygen",
      paste(source_type, oxygen, sep = " | "),
      paste(source_carbon, oxygen, sep = " | ")
    )
  ) %>%
  select(cell_definition, cell, mean_beta_16s, mean_beta_ko) %>%
  pivot_longer(c(mean_beta_16s, mean_beta_ko), names_to = "layer", values_to = "mean_beta") %>%
  mutate(layer = recode(layer, mean_beta_16s = "16S", mean_beta_ko = "KO")) %>%
  ggplot(aes(x = cell, y = mean_beta, fill = layer)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  facet_wrap(~cell_definition, scales = "free_x") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "Mean within-cell Bray-Curtis distance", fill = NULL)

ggsave(file.path(out_dir, "functional_redundancy_mean_beta_by_cell.png"), p_summary, width = 9.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "functional_redundancy_mean_beta_by_cell.pdf"), p_summary, width = 9.5, height = 4.8)

message("Matched condition-time units: ", nrow(matched_meta))
message("Missing 16S condition keys: ", length(setdiff(rownames(datko$mat), rownames(dat16$mat))))
message("Missing KO condition keys: ", length(setdiff(rownames(dat16$mat), rownames(datko$mat))))
message("Wrote outputs to: ", out_dir)
