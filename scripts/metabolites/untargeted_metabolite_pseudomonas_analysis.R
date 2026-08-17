#!/usr/bin/env Rscript

# Untargeted metabolite analysis for Recyc necromass Pseudomonas workbook.
# Raw workbook is left untouched; all derived tables and plots are written to OUT_DIR.

options(stringsAsFactors = FALSE)

PROJECT_DIR <- Sys.getenv("PROJECT_DIR", "/Users/mingfeichen/Recyc_necromass")
INPUT_XLSX <- Sys.getenv(
  "INPUT_XLSX",
  "/Users/mingfeichen/Recyc_necromass_untargeted_positive_pseudomonas.xlsx"
)
OUT_DIR <- Sys.getenv(
  "OUT_DIR",
  file.path(PROJECT_DIR, "outputs", "untargeted_metabolites_pseudomonas")
)

RT_MIN <- as.numeric(Sys.getenv("RT_MIN", "0.8"))
RT_MAX <- as.numeric(Sys.getenv("RT_MAX", "17.5"))
CONTROL_FOLD <- as.numeric(Sys.getenv("CONTROL_FOLD", "10"))
UPSET_THRESHOLD <- as.numeric(Sys.getenv("UPSET_THRESHOLD", "1000000"))
AEROBIC_ONLY <- tolower(Sys.getenv("AEROBIC_ONLY", "TRUE")) %in% c("1", "true", "t", "yes", "y")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(vegan)
  library(UpSetR)
  library(patchwork)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_metab <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey94", color = "grey70"),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = rel(1.15))
    )
}

oxygen_colors <- c("Aerobic" = "#2C7FB8", "Anaerobic" = "#D95F02", "Control" = "#6A51A3")
source_fills <- c("fresh" = "#F8766D", "recyc" = "#00BFC4", "na" = "#BDBDBD")
source_shapes <- c("fresh" = 21, "recyc" = 24, "na" = 22)
day_shapes <- c("d0" = 21, "d2" = 22, "d6" = 23, "d11" = 24, "d27" = 25)

if (AEROBIC_ONLY) {
  oxygen_colors <- c("Aerobic" = "#2C7FB8", "Control" = "#6A51A3")
}

save_plot <- function(plot, filename, width = 8, height = 6, dpi = 300) {
  pdf_path <- file.path(OUT_DIR, paste0(filename, ".pdf"))
  png_path <- file.path(OUT_DIR, paste0(filename, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, units = "in")
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf = pdf_path, png = png_path))
}

parse_sample_names <- function(cols) {
  tibble(column = cols) %>%
    mutate(
      sample_id = str_remove(column, "\\.mzML Peak height$"),
      sample_id = str_remove(sample_id, "^\\d+_"),
      sample_id = str_remove(sample_id, "__.*$"),
      replicate = as.integer(str_match(sample_id, "_([0-9]+)$")[, 2]),
      condition = str_remove(sample_id, "_[0-9]+$")
    ) %>%
    separate(
      condition,
      into = c("treatment", "carbon", "organism", "oxygen_raw", "source", "day", "inoculation"),
      sep = "-",
      remove = FALSE,
      fill = "right",
      extra = "merge"
    ) %>%
    mutate(
      is_control = treatment == "ExCtrl" | carbon == "noC" | str_detect(condition, "ExCtrl|noC"),
      is_experimental = !is_control,
      oxygen = case_when(
        oxygen_raw == "A" ~ "Aerobic",
        oxygen_raw == "Ana" ~ "Anaerobic",
        TRUE ~ "Control"
      ),
      day = factor(day, levels = c("d0", "d2", "d6", "d11", "d27"), ordered = TRUE),
      source = factor(source, levels = c("fresh", "recyc", "na")),
      group_id = paste(treatment, carbon, organism, oxygen_raw, source, day, inoculation, sep = "-"),
      comparison_key = paste(treatment, carbon, organism, oxygen_raw, source, day, sep = "-"),
      plot_label = paste(oxygen, as.character(source), as.character(day), inoculation, replicate, sep = " | ")
    )
}

feature_ids <- function(df) {
  df %>%
    transmute(
      feature_id = paste0("ID", `row ID`, "_mz", signif(`row m/z`, 6), "_rt", signif(`row retention time`, 4)),
      `row ID`,
      `row m/z`,
      `row retention time`
    )
}

safe_t_test <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  if (sd(x) == 0 && sd(y) == 0) return(ifelse(mean(x) == mean(y), 1, NA_real_))
  out <- tryCatch(t.test(x, y)$p.value, error = function(e) NA_real_)
  out
}

ordination_input <- function(values, sample_meta, sample_filter) {
  selected_meta <- sample_meta %>%
    filter({{ sample_filter }}) %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . } %>%
    mutate(
      day = factor(as.character(day), levels = c("d0", "d2", "d6", "d11", "d27"), ordered = TRUE),
      source = factor(as.character(source), levels = c("fresh", "recyc", "na")),
      oxygen = factor(oxygen, levels = c("Aerobic", "Anaerobic", "Control")),
      inoculation = factor(inoculation, levels = c("non", "inoc"))
    )

  if (nrow(selected_meta) < 3 || nrow(values) < 2) return(NULL)

  mat <- t(as.matrix(values[, selected_meta$column, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- log1p(mat)
  mat <- mat[, colSums(mat) > 0, drop = FALSE]
  if (ncol(mat) < 2) return(NULL)

  list(meta = selected_meta, mat = mat, dist = vegdist(mat, method = "bray"))
}

tidy_adonis <- function(adonis_result, sheet_name, model_name, formula_label, by_label) {
  as.data.frame(adonis_result) %>%
    tibble::rownames_to_column("term") %>%
    as_tibble() %>%
    filter(!term %in% c("Total")) %>%
    transmute(
      sheet = sheet_name,
      model = model_name,
      formula = formula_label,
      by = by_label,
      term,
      df = Df,
      sum_squares = SumOfSqs,
      r2 = R2,
      f_statistic = F,
      p_value = `Pr(>F)`
    )
}

safe_adonis <- function(dist_obj, meta, formula_obj, sheet_name, model_name, formula_label, by_label) {
  out <- tryCatch(
    adonis2(formula_obj, data = meta, permutations = 999, by = by_label),
    error = function(e) {
      message("PERMANOVA skipped for ", sheet_name, " / ", model_name, " / ", formula_label, ": ", e$message)
      NULL
    }
  )
  if (is.null(out)) return(tibble())
  tidy_adonis(out, sheet_name, model_name, formula_label, by_label)
}

safe_dispersion <- function(dist_obj, meta, variables, sheet_name, model_name) {
  bind_rows(lapply(variables, function(var_name) {
    group <- droplevels(factor(meta[[var_name]]))
    if (nlevels(group) < 2 || any(table(group) < 2)) return(tibble())

    bd <- tryCatch(betadisper(dist_obj, group), error = function(e) NULL)
    if (is.null(bd)) return(tibble())
    test <- tryCatch(permutest(bd, permutations = 999), error = function(e) NULL)
    if (is.null(test)) return(tibble())

    tibble(
      sheet = sheet_name,
      model = model_name,
      variable = var_name,
      df = test$tab[1, "Df"],
      f_statistic = test$tab[1, "F"],
      p_value = test$tab[1, "Pr(>F)"]
    )
  }))
}

make_permanova <- function(values, sample_meta, sheet_name) {
  all_input <- ordination_input(
    values,
    sample_meta,
    is_experimental &
      ((inoculation == "inoc" & as.character(day) != "d0") |
        (inoculation == "non" & as.character(day) == "d0"))
  )

  inoc_input <- ordination_input(
    values,
    sample_meta,
    is_experimental & inoculation == "inoc" & as.character(day) != "d0"
  )

  permanova_tables <- list()
  dispersion_tables <- list()

  if (!is.null(all_input)) {
    meta <- all_input$meta %>%
      mutate(
        source = droplevels(source),
        day = droplevels(day),
        inoculation = droplevels(inoculation)
      )
    dist_obj <- all_input$dist
    formula_all <- if (AEROBIC_ONLY) {
      dist_obj ~ source + day + inoculation
    } else {
      dist_obj ~ source + oxygen + day + inoculation
    }
    permanova_tables$all_additive <- safe_adonis(
      dist_obj,
      meta,
      formula_all,
      sheet_name,
      "ordination_selected",
      if (AEROBIC_ONLY) "Bray ~ source + day + inoculation" else "Bray ~ source + oxygen + day + inoculation",
      "margin"
    )
    dispersion_tables$all <- safe_dispersion(
      dist_obj,
      meta,
      if (AEROBIC_ONLY) c("source", "day", "inoculation") else c("source", "oxygen", "day", "inoculation"),
      sheet_name,
      "ordination_selected"
    )
  }

  if (!is.null(inoc_input)) {
    meta <- inoc_input$meta %>%
      mutate(
        source = droplevels(source),
        day = droplevels(day)
      )
    dist_obj <- inoc_input$dist
    formula_inoc <- if (AEROBIC_ONLY) {
      dist_obj ~ source + day
    } else {
      dist_obj ~ source + oxygen + day
    }
    permanova_tables$inoc_additive <- safe_adonis(
      dist_obj,
      meta,
      formula_inoc,
      sheet_name,
      "inoculated_only",
      if (AEROBIC_ONLY) "Bray ~ source + day" else "Bray ~ source + oxygen + day",
      "margin"
    )
    formula_inoc_terms <- if (AEROBIC_ONLY) {
      dist_obj ~ source * day
    } else {
      dist_obj ~ source * oxygen * day
    }
    permanova_tables$inoc_interaction <- safe_adonis(
      dist_obj,
      meta,
      formula_inoc_terms,
      sheet_name,
      "inoculated_only_interactions",
      if (AEROBIC_ONLY) "Bray ~ source * day" else "Bray ~ source * oxygen * day",
      "terms"
    )
    dispersion_tables$inoc <- safe_dispersion(
      dist_obj,
      meta,
      if (AEROBIC_ONLY) c("source", "day") else c("source", "oxygen", "day"),
      sheet_name,
      "inoculated_only"
    )
  }

  permanova <- bind_rows(permanova_tables)
  dispersion <- bind_rows(dispersion_tables)

  if (nrow(permanova) > 0) {
    write_csv(permanova, file.path(OUT_DIR, paste0(sheet_name, "_permanova_results.csv")))

    plot_df <- permanova %>%
      filter(
        model %in% c("ordination_selected", "inoculated_only"),
        !term %in% c("Residual")
      ) %>%
      mutate(
        term = factor(term, levels = if (AEROBIC_ONLY) c("source", "day", "inoculation") else c("source", "oxygen", "day", "inoculation")),
        p_label = if_else(is.na(p_value), "p = NA", paste0("p = ", signif(p_value, 2)))
      )

    if (nrow(plot_df) > 0) {
      p <- ggplot(plot_df, aes(term, r2, fill = term)) +
        geom_col(width = 0.65, show.legend = FALSE) +
        geom_text(aes(label = p_label), vjust = -0.25, size = 2.7) +
        facet_wrap(~ model, scales = "free_x") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
        labs(
          title = paste(sheet_name, "PERMANOVA variable contributions"),
          subtitle = "Bray-Curtis distances on log1p filtered peak heights; bars show marginal R2",
          x = NULL,
          y = "PERMANOVA R2"
        ) +
        theme_metab(base_size = 9) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
      save_plot(p, paste0(sheet_name, "_permanova_variable_r2"), width = 8.5, height = 4.8)
    }
  }

  if (nrow(dispersion) > 0) {
    write_csv(dispersion, file.path(OUT_DIR, paste0(sheet_name, "_betadisper_results.csv")))
  }

  list(permanova = permanova, dispersion = dispersion)
}

make_bray_pcoa <- function(values, sample_meta, sheet_name) {
  selected_meta <- sample_meta %>%
    filter(
      is_experimental,
      (inoculation == "inoc" & as.character(day) != "d0") |
        (inoculation == "non" & as.character(day) == "d0")
    ) %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . }

  if (nrow(selected_meta) < 3 || nrow(values) < 2) return(tibble())

  mat <- t(as.matrix(values[, selected_meta$column, drop = FALSE]))
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  mat <- log1p(mat)
  keep_features <- colSums(mat) > 0
  mat <- mat[, keep_features, drop = FALSE]
  if (ncol(mat) < 2) return(tibble())

  dist_obj <- vegdist(mat, method = "bray")
  pcoa <- cmdscale(dist_obj, k = 2, eig = TRUE)
  positive_eig <- pcoa$eig[pcoa$eig > 0]
  pct <- if (length(positive_eig) >= 2) {
    round(100 * positive_eig[1:2] / sum(positive_eig), 1)
  } else {
    c(NA_real_, NA_real_)
  }

  scores <- selected_meta %>%
    mutate(
      PCoA1 = pcoa$points[, 1],
      PCoA2 = pcoa$points[, 2],
      sheet = sheet_name
    )

  if (AEROBIC_ONLY) {
    p <- ggplot(scores, aes(PCoA1, PCoA2)) +
      geom_point(
        aes(color = source, shape = day),
        size = 3.5,
        stroke = 0.9,
        alpha = 0.95
      ) +
      scale_color_manual(
        values = c("fresh" = source_fills["fresh"], "recyc" = source_fills["recyc"]),
        labels = c("fresh" = "Fresh", "recyc" = "Recycled"),
        drop = FALSE
      ) +
      scale_shape_manual(values = day_shapes, drop = FALSE) +
      labs(
        title = paste(sheet_name, "aerobic filtered features: Bray-Curtis PCoA"),
        x = paste0("PCoA1 (", pct[1], "%)"),
        y = paste0("PCoA2 (", pct[2], "%)"),
        color = "Source",
        shape = "Day"
      ) +
      theme_metab()
  } else {
    p <- ggplot(scores, aes(PCoA1, PCoA2)) +
      geom_point(
        aes(color = oxygen, shape = day, fill = source),
        size = 3.5,
        stroke = 0.9,
        alpha = 0.95
      ) +
      scale_color_manual(values = oxygen_colors, drop = FALSE) +
      scale_shape_manual(values = day_shapes, drop = FALSE) +
      scale_fill_manual(
        values = source_fills,
        drop = TRUE,
        guide = guide_legend(override.aes = list(shape = 21, color = "grey25", size = 4))
      ) +
      labs(
        title = paste(sheet_name, "filtered features: Bray-Curtis PCoA"),
        x = paste0("PCoA1 (", pct[1], "%)"),
        y = paste0("PCoA2 (", pct[2], "%)"),
        color = "Aeration",
        shape = "Day",
        fill = "Source"
      ) +
      theme_metab()
  }

  save_plot(p, paste0(sheet_name, "_bray_pcoa_selected_samples"), width = 8.5, height = 6.2)
  scores
}

volcano_classes <- function(log2fc, p_value, up_label, down_label) {
  case_when(
    !is.na(p_value) & p_value < 0.05 & log2fc >= 1 ~ up_label,
    !is.na(p_value) & p_value < 0.05 & log2fc <= -1 ~ down_label,
    TRUE ~ "Not significant"
  )
}

volcano_counts <- function(stats, up_label, down_label) {
  counts <- stats %>%
    count(volcano_class, name = "n") %>%
    tidyr::complete(
      volcano_class = c(up_label, down_label, "Not significant"),
      fill = list(n = 0)
    )

  tibble(
    increase = counts$n[match(up_label, counts$volcano_class)],
    decrease = counts$n[match(down_label, counts$volcano_class)],
    ns = counts$n[match("Not significant", counts$volcano_class)]
  )
}

volcano_count_label <- function(stats, up_label, down_label, up_short = "increase", down_short = "decrease") {
  counts <- volcano_counts(stats, up_label, down_label)
  paste0(
    up_short, ": ", counts$increase,
    " | ", down_short, ": ", counts$decrease,
    " | ns: ", counts$ns
  )
}

make_single_volcano_plot <- function(stats, sheet_name, comparison_label, x_label, title_suffix,
                                     filename, up_label, down_label, up_short, down_short) {
  label_df <- stats %>%
    filter(!is.na(p_value), is.finite(log2FC)) %>%
    arrange(p_value, desc(abs(log2FC))) %>%
    slice_head(n = 10) %>%
    mutate(label = paste0("ID", `row ID`, "\nrt ", signif(`row retention time`, 3)))

  p <- ggplot(stats, aes(log2FC, neg_log10_p)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60", linewidth = 0.35) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.35) +
    geom_point(aes(color = volcano_class), alpha = 0.65, size = 1.5) +
    geom_text_repel(
      data = label_df,
      aes(label = label),
      size = 2.6,
      max.overlaps = 20,
      min.segment.length = 0
    ) +
    scale_color_manual(
      values = c(
        up_label = "#B2182B",
        down_label = "#2166AC",
        "Not significant" = "grey55"
      ) %>% setNames(c(up_label, down_label, "Not significant"))
    ) +
    labs(
      title = paste(sheet_name, comparison_label, title_suffix, sep = ": "),
      subtitle = volcano_count_label(stats, up_label, down_label, up_short, down_short),
      x = x_label,
      y = "-log10 p-value",
      color = NULL
    ) +
    theme_metab()

  save_plot(p, filename, width = 7.2, height = 6)
}

make_volcanoes <- function(values, feature_table, sample_meta, sheet_name) {
  old_volcano_files <- list.files(
    OUT_DIR,
    pattern = paste0("^", sheet_name, "_volcano_.*\\.(pdf|png)$"),
    full.names = TRUE
  )
  if (length(old_volcano_files) > 0) unlink(old_volcano_files)

  # Media controls only exist for aerobic -non samples. Compare each inoculated
  # group, including anaerobic groups, to the same source/day aerobic -non group.
  non_groups <- sample_meta %>%
    filter(is_experimental, inoculation == "non", oxygen_raw == "A") %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . } %>%
    distinct(non_key = comparison_key)

  comparison_pairs <- sample_meta %>%
    filter(is_experimental, inoculation == "inoc") %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . } %>%
    distinct(treatment, carbon, organism, oxygen_raw, oxygen, source, day, inoc_key = comparison_key) %>%
    mutate(
      non_key = paste(treatment, carbon, organism, "A", source, day, sep = "-"),
      comparison = paste0(inoc_key, "_inoc_vs_", non_key, "_non"),
      comparison_label = paste0(inoc_key, " inoc vs ", non_key, " non")
    ) %>%
    inner_join(non_groups, by = "non_key") %>%
    arrange(day, source, oxygen_raw)

  if (nrow(comparison_pairs) == 0) return(tibble())

  all_stats <- vector("list", nrow(comparison_pairs))
  names(all_stats) <- comparison_pairs$comparison

  for (pair_idx in seq_len(nrow(comparison_pairs))) {
    pair <- comparison_pairs[pair_idx, ]
    inoc_cols <- sample_meta %>%
      filter(comparison_key == pair$inoc_key, inoculation == "inoc") %>%
      arrange(replicate) %>%
      pull(column)
    non_cols <- sample_meta %>%
      filter(comparison_key == pair$non_key, inoculation == "non") %>%
      arrange(replicate) %>%
      pull(column)

    stats <- feature_table %>%
      mutate(
        inoc_mean = rowMeans(values[, inoc_cols, drop = FALSE], na.rm = TRUE),
        non_mean = rowMeans(values[, non_cols, drop = FALSE], na.rm = TRUE),
        log2FC = log2((inoc_mean + 1) / (non_mean + 1)),
        p_value = vapply(seq_len(nrow(values)), function(i) {
          safe_t_test(
            log2(as.numeric(values[i, inoc_cols, drop = TRUE]) + 1),
            log2(as.numeric(values[i, non_cols, drop = TRUE]) + 1)
          )
        }, numeric(1)),
        comparison = pair$comparison,
        inoc_group = pair$inoc_key,
        non_group = pair$non_key,
        non_reference = "aerobic -non",
        inoc_oxygen = pair$oxygen,
        inoc_source = as.character(pair$source),
        inoc_day = as.character(pair$day),
        neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
        significant = !is.na(p_value) & p_value < 0.05 & abs(log2FC) >= 1,
      volcano_class = volcano_classes(log2FC, p_value, "Higher in inoc", "Higher in aerobic non")
    ) %>%
      group_by(comparison) %>%
      mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
      ungroup()

    all_stats[[pair$comparison]] <- stats

    file_key <- str_replace_all(pair$comparison, "[^A-Za-z0-9]+", "_")
    make_single_volcano_plot(
      stats,
      sheet_name,
      pair$comparison_label,
      "log2FC inoc / aerobic non",
      "inoc vs same-day aerobic non",
      paste0(sheet_name, "_volcano_", file_key),
      "Higher in inoc",
      "Higher in aerobic non",
      "increase",
      "decrease"
    )
  }

  bind_rows(all_stats)
}

make_combined_volcano <- function(volcano_stats, sheet_name) {
  if (nrow(volcano_stats) == 0) return(invisible(NULL))

  plot_df <- volcano_stats %>%
    mutate(
      inoc_day = factor(inoc_day, levels = c("d2", "d6", "d11", "d27"), ordered = TRUE),
      inoc_oxygen = factor(inoc_oxygen, levels = if (AEROBIC_ONLY) c("Aerobic") else c("Aerobic", "Anaerobic")),
      inoc_source = factor(inoc_source, levels = c("fresh", "recyc")),
      volcano_class = case_when(
        !is.na(p_value) & p_value < 0.05 & log2FC >= 1 ~ "Higher in inoc",
        !is.na(p_value) & p_value < 0.05 & log2FC <= -1 ~ "Higher in aerobic non",
        TRUE ~ "Not significant"
      )
    )

  count_labels <- plot_df %>%
    group_by(inoc_source, inoc_oxygen, inoc_day) %>%
    summarise(
      label = volcano_count_label(pick(everything()), "Higher in inoc", "Higher in aerobic non", "inc", "dec"),
      .groups = "drop"
    )

  p <- ggplot(plot_df, aes(log2FC, neg_log10_p)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.25) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.25) +
    geom_point(aes(color = volcano_class), alpha = 0.45, size = 0.55) +
    geom_text(
      data = count_labels,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.04,
      vjust = 1.25,
      size = 2.1,
      color = "grey20"
    ) +
    { if (AEROBIC_ONLY) facet_grid(inoc_source ~ inoc_day) else facet_grid(inoc_source ~ inoc_oxygen + inoc_day) } +
    scale_color_manual(
      values = c(
        "Higher in inoc" = "#B2182B",
        "Higher in aerobic non" = "#2166AC",
        "Not significant" = "grey68"
      )
    ) +
    labs(
      title = paste(sheet_name, "combined volcano plots"),
      subtitle = "Each inoculated group is compared with the matching aerobic -non media control",
      x = "log2FC inoc / aerobic non",
      y = "-log10 p-value",
      color = NULL
    ) +
    theme_metab(base_size = 8.5) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 7.5),
      axis.text = element_text(size = 6.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2)
    )

  save_plot(p, paste0(sheet_name, "_combined_volcano_all_timepoints_aerobic_anaerobic"), width = 18, height = 7.5)
  invisible(NULL)
}

make_d0_volcanoes <- function(values, feature_table, sample_meta, sheet_name) {
  old_d0_volcano_files <- list.files(
    OUT_DIR,
    pattern = paste0("^", sheet_name, "_d0_volcano_.*\\.(pdf|png)$"),
    full.names = TRUE
  )
  if (length(old_d0_volcano_files) > 0) unlink(old_d0_volcano_files)

  d0_groups <- sample_meta %>%
    filter(
      is_experimental,
      treatment == "TxCtrl",
      inoculation == "non",
      oxygen_raw == "A",
      as.character(day) == "d0"
    ) %>%
    distinct(d0_key = comparison_key)

  comparison_pairs <- sample_meta %>%
    filter(is_experimental, inoculation == "inoc", as.character(day) != "d0") %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . } %>%
    distinct(treatment, carbon, organism, oxygen_raw, oxygen, source, day, inoc_key = comparison_key) %>%
    mutate(
      d0_key = paste("TxCtrl", carbon, organism, "A", source, "d0", sep = "-"),
      comparison = paste0(inoc_key, "_inoc_vs_", d0_key, "_non"),
      comparison_label = paste0(inoc_key, " inoc vs ", d0_key, " non")
    ) %>%
    inner_join(d0_groups, by = "d0_key") %>%
    arrange(day, source, oxygen_raw)

  if (nrow(comparison_pairs) == 0) return(tibble())

  all_stats <- vector("list", nrow(comparison_pairs))
  names(all_stats) <- comparison_pairs$comparison

  for (pair_idx in seq_len(nrow(comparison_pairs))) {
    pair <- comparison_pairs[pair_idx, ]
    inoc_cols <- sample_meta %>%
      filter(comparison_key == pair$inoc_key, inoculation == "inoc") %>%
      arrange(replicate) %>%
      pull(column)
    d0_cols <- sample_meta %>%
      filter(comparison_key == pair$d0_key, inoculation == "non") %>%
      arrange(replicate) %>%
      pull(column)

    stats <- feature_table %>%
      mutate(
        inoc_mean = rowMeans(values[, inoc_cols, drop = FALSE], na.rm = TRUE),
        d0_mean = rowMeans(values[, d0_cols, drop = FALSE], na.rm = TRUE),
        log2FC = log2((inoc_mean + 1) / (d0_mean + 1)),
        p_value = vapply(seq_len(nrow(values)), function(i) {
          safe_t_test(
            log2(as.numeric(values[i, inoc_cols, drop = TRUE]) + 1),
            log2(as.numeric(values[i, d0_cols, drop = TRUE]) + 1)
          )
        }, numeric(1)),
        comparison = pair$comparison,
        inoc_group = pair$inoc_key,
        d0_group = pair$d0_key,
        reference = "aerobic d0 -non",
        inoc_oxygen = pair$oxygen,
        inoc_source = as.character(pair$source),
        inoc_day = as.character(pair$day),
        neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
        significant = !is.na(p_value) & p_value < 0.05 & abs(log2FC) >= 1,
        volcano_class = volcano_classes(log2FC, p_value, "Higher than d0", "Lower than d0")
      ) %>%
      group_by(comparison) %>%
      mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%
      ungroup()

    all_stats[[pair$comparison]] <- stats

    file_key <- str_replace_all(pair$comparison, "[^A-Za-z0-9]+", "_")
    make_single_volcano_plot(
      stats,
      sheet_name,
      pair$comparison_label,
      "log2FC inoc / d0",
      "inoc vs d0",
      paste0(sheet_name, "_d0_volcano_", file_key),
      "Higher than d0",
      "Lower than d0",
      "increase",
      "decrease"
    )
  }

  bind_rows(all_stats)
}

make_combined_d0_volcano <- function(d0_volcano_stats, sheet_name) {
  if (nrow(d0_volcano_stats) == 0) return(invisible(NULL))

  plot_df <- d0_volcano_stats %>%
    mutate(
      inoc_day = factor(inoc_day, levels = c("d2", "d6", "d11", "d27"), ordered = TRUE),
      inoc_oxygen = factor(inoc_oxygen, levels = if (AEROBIC_ONLY) c("Aerobic") else c("Aerobic", "Anaerobic")),
      inoc_source = factor(inoc_source, levels = c("fresh", "recyc")),
      volcano_class = volcano_classes(log2FC, p_value, "Higher than d0", "Lower than d0")
    )

  count_labels <- plot_df %>%
    group_by(inoc_source, inoc_oxygen, inoc_day) %>%
    summarise(
      label = volcano_count_label(pick(everything()), "Higher than d0", "Lower than d0", "inc", "dec"),
      .groups = "drop"
    )

  p <- ggplot(plot_df, aes(log2FC, neg_log10_p)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.25) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.25) +
    geom_point(aes(color = volcano_class), alpha = 0.45, size = 0.55) +
    geom_text(
      data = count_labels,
      aes(x = Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = 1.04,
      vjust = 1.25,
      size = 2.1,
      color = "grey20"
    ) +
    { if (AEROBIC_ONLY) facet_grid(inoc_source ~ inoc_day) else facet_grid(inoc_source ~ inoc_oxygen + inoc_day) } +
    scale_color_manual(
      values = c(
        "Higher than d0" = "#B2182B",
        "Lower than d0" = "#2166AC",
        "Not significant" = "grey68"
      )
    ) +
    labs(
      title = paste(sheet_name, "combined d0-reference volcano plots"),
      subtitle = "Each inoculated group is compared with the matching same-source aerobic d0 -non control",
      x = "log2FC inoc / d0",
      y = "-log10 p-value",
      color = NULL
    ) +
    theme_metab(base_size = 8.5) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(size = 7.5),
      axis.text = element_text(size = 6.5),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2)
    )

  save_plot(p, paste0(sheet_name, "_combined_d0_volcano_all_timepoints_aerobic_anaerobic"), width = 18, height = 7.5)
  invisible(NULL)
}

make_upset <- function(values, feature_table, sample_meta, sheet_name) {
  selected_groups <- sample_meta %>%
    filter(
      is_experimental,
      (inoculation == "inoc" & as.character(day) != "d0") |
        (inoculation == "non" & as.character(day) == "d0")
    ) %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . } %>%
    arrange(day, oxygen, source, inoculation, replicate)

  if (nrow(selected_groups) == 0) return(tibble())

  group_order <- selected_groups %>%
    distinct(group_id, day, oxygen, source, inoculation) %>%
    arrange(day, oxygen, source, inoculation) %>%
    pull(group_id)

  presence <- feature_table
  for (grp in group_order) {
    cols <- selected_groups %>% filter(group_id == grp) %>% pull(column)
    presence[[grp]] <- as.integer(apply(values[, cols, drop = FALSE], 1, max, na.rm = TRUE) > UPSET_THRESHOLD)
  }

  upset_matrix <- presence %>% select(all_of(group_order))
  if (ncol(upset_matrix) > 1 && sum(as.matrix(upset_matrix)) > 0) {
    pdf(file.path(OUT_DIR, paste0(sheet_name, "_upset_presence_gt_1e6.pdf")), width = 13, height = 7)
    print(
      upset(
        as.data.frame(upset_matrix),
        sets = rev(group_order),
        nsets = length(group_order),
        nintersects = 35,
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Feature intersections",
        sets.x.label = "Features > 1e6 in >=1 replicate"
      )
    )
    dev.off()

    png(file.path(OUT_DIR, paste0(sheet_name, "_upset_presence_gt_1e6.png")), width = 3900, height = 2100, res = 300)
    print(
      upset(
        as.data.frame(upset_matrix),
        sets = rev(group_order),
        nsets = length(group_order),
        nintersects = 35,
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Feature intersections",
        sets.x.label = "Features > 1e6 in >=1 replicate"
      )
    )
    dev.off()
  }

  presence
}

make_extra_plots <- function(values, feature_table, sample_meta, sheet_name) {
  max_experimental <- apply(values[, sample_meta$is_experimental, drop = FALSE], 1, max, na.rm = TRUE)
  rt_mz <- feature_table %>%
    mutate(max_experimental = max_experimental)

  p_rt <- ggplot(rt_mz, aes(`row retention time`, `row m/z`, color = log10(max_experimental + 1))) +
    geom_point(alpha = 0.7, size = 1.4) +
    scale_color_viridis_c(option = "magma") +
    labs(
      title = paste(sheet_name, "selected feature RT/m/z distribution"),
      x = "Retention time",
      y = "m/z",
      color = "log10 max\nexperimental"
    ) +
    theme_metab()
  save_plot(p_rt, paste0(sheet_name, "_rt_mz_selected_features"), width = 7.5, height = 5.8)

  selected_meta <- sample_meta %>%
    filter(
      is_experimental,
      (inoculation == "inoc" & as.character(day) != "d0") |
        (inoculation == "non" & as.character(day) == "d0")
    ) %>%
    { if (AEROBIC_ONLY) filter(., oxygen == "Aerobic") else . }

  long_selected <- values %>%
    mutate(feature_id = feature_table$feature_id) %>%
    pivot_longer(all_of(selected_meta$column), names_to = "column", values_to = "peak_height") %>%
    left_join(selected_meta, by = "column")

  feature_counts <- long_selected %>%
    group_by(group_id, day, oxygen, source, inoculation) %>%
    summarise(features_gt_1e6 = sum(peak_height > UPSET_THRESHOLD, na.rm = TRUE), .groups = "drop") %>%
    mutate(group_id = factor(group_id, levels = unique(group_id[order(day, oxygen, source)])))

  p_counts <- ggplot(feature_counts, aes(group_id, features_gt_1e6, fill = oxygen)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = oxygen_colors, drop = FALSE) +
    labs(
      title = paste(sheet_name, "features above 1e6 by selected group"),
      x = NULL,
      y = "Feature-replicate count > 1e6",
      fill = "Aeration"
    ) +
    theme_metab()
  save_plot(p_counts, paste0(sheet_name, "_feature_counts_gt_1e6_by_group"), width = 9, height = 7)

  sample_totals <- long_selected %>%
    group_by(column, group_id, day, oxygen, source, inoculation, replicate) %>%
    summarise(total_intensity = sum(peak_height, na.rm = TRUE), .groups = "drop")

  p_totals <- ggplot(sample_totals, aes(day, log10(total_intensity + 1))) +
    geom_boxplot(aes(fill = source), outlier.shape = NA, alpha = 0.55) +
    geom_point(aes(color = oxygen, shape = source), position = position_jitter(width = 0.12, height = 0), size = 2.4) +
    facet_wrap(~ oxygen, nrow = 1) +
    scale_color_manual(values = oxygen_colors, drop = FALSE) +
    scale_fill_manual(values = source_fills, drop = FALSE) +
    labs(
      title = paste(sheet_name, "sample total selected-feature intensity"),
      x = "Day",
      y = "log10 summed peak height",
      fill = "Source",
      shape = "Source",
      color = "Aeration"
    ) +
    theme_metab()
  save_plot(p_totals, paste0(sheet_name, "_selected_sample_total_intensity"), width = 8.5, height = 5)

  condition_means <- long_selected %>%
    group_by(feature_id, group_id, day, oxygen, source, inoculation) %>%
    summarise(mean_peak = mean(peak_height, na.rm = TRUE), .groups = "drop")

  top_features <- condition_means %>%
    group_by(feature_id) %>%
    summarise(var_log = var(log10(mean_peak + 1), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(var_log)) %>%
    slice_head(n = 50) %>%
    pull(feature_id)

  heatmap_df <- condition_means %>%
    filter(feature_id %in% top_features) %>%
    group_by(feature_id) %>%
    mutate(z = as.numeric(scale(log10(mean_peak + 1)))) %>%
    ungroup() %>%
    mutate(
      feature_id = factor(feature_id, levels = rev(top_features)),
      group_id = factor(group_id, levels = unique(group_id[order(day, oxygen, source)]))
    )

  p_heat <- ggplot(heatmap_df, aes(group_id, feature_id, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", na.value = "grey90") +
    labs(
      title = paste(sheet_name, "top 50 variable selected features"),
      x = NULL,
      y = NULL,
      fill = "row z-score"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )
  save_plot(p_heat, paste0(sheet_name, "_top50_variable_feature_heatmap"), width = 12, height = 9)

  list(feature_counts = feature_counts, sample_totals = sample_totals)
}

analyze_sheet <- function(sheet_name) {
  message("Analyzing sheet: ", sheet_name)
  raw <- read_excel(INPUT_XLSX, sheet = sheet_name)
  sample_cols <- setdiff(names(raw), c("row ID", "row m/z", "row retention time"))
  sample_meta <- parse_sample_names(sample_cols)

  values_all <- raw %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), as.numeric))

  control_cols <- sample_meta %>% filter(is_control) %>% pull(column)
  experimental_cols <- sample_meta %>% filter(is_experimental) %>% pull(column)
  analysis_meta <- sample_meta
  if (AEROBIC_ONLY) {
    analysis_meta <- sample_meta %>% filter(oxygen == "Aerobic")
  }

  max_control <- apply(values_all[, control_cols, drop = FALSE], 1, max, na.rm = TRUE)
  max_experimental <- apply(values_all[, experimental_cols, drop = FALSE], 1, max, na.rm = TRUE)
  max_control_for_ratio <- pmax(max_control, 1)

  metrics <- raw %>%
    select(`row ID`, `row m/z`, `row retention time`) %>%
    mutate(
      max_control = max_control,
      max_experimental = max_experimental,
      experimental_control_ratio = max_experimental / max_control_for_ratio,
      pass_rt = `row retention time` >= RT_MIN & `row retention time` <= RT_MAX,
      pass_control_fold = max_experimental > CONTROL_FOLD * max_control_for_ratio,
      pass_filter = pass_rt & pass_control_fold
    )

  filtered_raw <- raw %>%
    mutate(
      max_control = max_control,
      max_experimental = max_experimental,
      experimental_control_ratio = max_experimental / max_control_for_ratio
    ) %>%
    filter(metrics$pass_filter)

  filtered_values <- filtered_raw %>%
    select(all_of(sample_cols)) %>%
    mutate(across(everything(), as.numeric))

  fids <- feature_ids(filtered_raw)
  filtered_feature_table <- filtered_raw %>%
    select(`row ID`, `row m/z`, `row retention time`, max_control, max_experimental, experimental_control_ratio) %>%
    mutate(feature_id = fids$feature_id, .before = 1)

  write_csv(analysis_meta, file.path(OUT_DIR, paste0(sheet_name, "_sample_metadata.csv")))
  write_csv(metrics, file.path(OUT_DIR, paste0(sheet_name, "_filtering_metrics_all_features.csv")))
  write_csv(filtered_feature_table, file.path(OUT_DIR, paste0(sheet_name, "_filtered_feature_metrics.csv")))
  write_csv(filtered_raw, file.path(OUT_DIR, paste0(sheet_name, "_filtered_peak_heights.csv")))

  pcoa_scores <- make_bray_pcoa(filtered_values, analysis_meta, sheet_name)
  permanova_results <- make_permanova(filtered_values, analysis_meta, sheet_name)
  volcano_stats <- make_volcanoes(filtered_values, filtered_feature_table, analysis_meta, sheet_name)
  d0_volcano_stats <- make_d0_volcanoes(filtered_values, filtered_feature_table, analysis_meta, sheet_name)
  upset_presence <- make_upset(filtered_values, filtered_feature_table, analysis_meta, sheet_name)
  extras <- make_extra_plots(filtered_values, filtered_feature_table, analysis_meta, sheet_name)
  make_combined_volcano(volcano_stats, sheet_name)
  make_combined_d0_volcano(d0_volcano_stats, sheet_name)

  if (nrow(pcoa_scores) > 0) {
    write_csv(pcoa_scores, file.path(OUT_DIR, paste0(sheet_name, "_bray_pcoa_scores.csv")))
  }
  if (nrow(volcano_stats) > 0) {
    write_csv(volcano_stats, file.path(OUT_DIR, paste0(sheet_name, "_volcano_stats.csv")))
  }
  if (nrow(d0_volcano_stats) > 0) {
    write_csv(d0_volcano_stats, file.path(OUT_DIR, paste0(sheet_name, "_d0_volcano_stats.csv")))
  }
  if (nrow(upset_presence) > 0) {
    write_csv(upset_presence, file.path(OUT_DIR, paste0(sheet_name, "_upset_presence_matrix_gt_1e6.csv")))
  }
  write_csv(extras$feature_counts, file.path(OUT_DIR, paste0(sheet_name, "_feature_counts_gt_1e6_by_group.csv")))
  write_csv(extras$sample_totals, file.path(OUT_DIR, paste0(sheet_name, "_selected_sample_total_intensity.csv")))

  tibble(
    sheet = sheet_name,
    raw_features = nrow(raw),
    rt_filtered_features = sum(metrics$pass_rt, na.rm = TRUE),
    selected_features = sum(metrics$pass_filter, na.rm = TRUE),
    control_columns = length(control_cols),
    experimental_columns = length(experimental_cols),
    pcoa_samples = nrow(pcoa_scores),
    permanova_terms = nrow(permanova_results$permanova),
    volcano_comparisons = n_distinct(volcano_stats$comparison),
    d0_volcano_comparisons = n_distinct(d0_volcano_stats$comparison),
    upset_groups = max(0, ncol(upset_presence) - ncol(filtered_feature_table))
  )
}

sheet_names <- excel_sheets(INPUT_XLSX)
summary_tbl <- bind_rows(lapply(sheet_names, analyze_sheet))
write_csv(summary_tbl, file.path(OUT_DIR, "filtering_summary.csv"))

combined_permanova <- bind_rows(lapply(sheet_names, function(sheet_name) {
  path <- file.path(OUT_DIR, paste0(sheet_name, "_permanova_results.csv"))
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else tibble()
}))
if (nrow(combined_permanova) > 0) {
  write_csv(combined_permanova, file.path(OUT_DIR, "combined_permanova_results.csv"))
}

combined_betadisper <- bind_rows(lapply(sheet_names, function(sheet_name) {
  path <- file.path(OUT_DIR, paste0(sheet_name, "_betadisper_results.csv"))
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else tibble()
}))
if (nrow(combined_betadisper) > 0) {
  write_csv(combined_betadisper, file.path(OUT_DIR, "combined_betadisper_results.csv"))
}

p_summary <- ggplot(summary_tbl, aes(sheet, selected_features, fill = sheet)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = selected_features), vjust = -0.35, size = 3.8) +
  labs(
    title = "Filtered feature counts by ionization tab",
    subtitle = paste0(
      "RT ", RT_MIN, "-", RT_MAX,
      "; experimental max > ", CONTROL_FOLD, "x control max"
    ),
    x = NULL,
    y = "Selected features"
  ) +
  theme_metab()
save_plot(p_summary, "filtering_summary_selected_features", width = 6, height = 4.5)

message("Wrote outputs to: ", OUT_DIR)
