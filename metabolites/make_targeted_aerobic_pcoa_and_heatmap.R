#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(ggplot2)
  library(vegan)
  library(grid)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(svglite)
  library(ragg)
})

PROJECT_DIR <- "/Users/mingfeichen/Recyc_necromass"
PSEU_PATH <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_filtered", "Pseudo_targeted_metabolites_051926_filtered_3x.csv")
ARTH_PATH <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_filtered", "Arthrobacter_targeted_metabolites_filtered_3x.csv")

PSEU_OUT <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_pseudomonas")
ARTH_OUT <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_arthrobacter")
HEAT_OUT <- file.path(PROJECT_DIR, "outputs", "targeted_metabolites_composite")

dir.create(PSEU_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(ARTH_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(HEAT_OUT, recursive = TRUE, showWarnings = FALSE)

BLUE <- "#2C7FB8"
GREEN <- "#2CA25F"
SOURCE_FILL <- c(fresh = "#FFFFFF", recyc = "#2C7FB8")
SOURCE_FILL_D0 <- c(fresh = "#FFFFFF", recyc = "#2CA25F")
DAY_SHAPES <- c(d0 = 22, d2 = 21, d6 = 24, d11 = 23, d27 = 25)
DAY_ORDER <- c("d0", "d2", "d6", "d11", "d27")

CLASS_PALETTE <- c(
  "amino acid/deriv." = "#8CB9E8",
  "nucleic acid constituent" = "#E88F8C",
  "organic acid" = "#C9B4DE",
  "vitamin/cofactor" = "#F3C86A",
  "carbohydrate" = "#93D08F",
  "other" = "#D9D9D9"
)

theme_set(
  theme_classic(base_size = 11.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.48, colour = "black"),
      axis.ticks = element_line(linewidth = 0.48, colour = "black"),
      axis.title = element_text(size = 12.5),
      axis.text = element_text(size = 10.5, colour = "black"),
      strip.text = element_text(size = 11.2, face = "bold"),
      legend.title = element_text(size = 11.3, face = "bold"),
      legend.text = element_text(size = 10.2),
      plot.title = element_text(size = 13.5, face = "bold"),
      plot.subtitle = element_text(size = 10.0),
      plot.caption = element_text(size = 9.2, colour = "#555555"),
      panel.grid = element_blank()
    )
)

normalize_label <- function(x) {
  x <- str_squish(as.character(x))
  x[str_to_lower(x) == "malic acid peak 1"] <- "Malic acid"
  x
}

parse_sample_meta <- function(df) {
  cols <- names(df)[grepl("^\\d+_.+_\\d+$", names(df))]
  meta <- lapply(cols, function(col) {
    m <- str_match(col, "^\\d+_(.+)_(\\d+)$")
    sample_name <- m[, 2]
    replicate <- as.integer(m[, 3])
    parts <- str_split(sample_name, "-", simplify = TRUE)
    if (ncol(parts) < 7) {
      return(NULL)
    }
    if (parts[1, 1] == "TxCtrl") {
      type <- "d0"
    } else if (parts[1, 1] == "sup") {
      type <- "sup"
    } else {
      return(NULL)
    }
    tibble(
      column = col,
      sample_name = sample_name,
      replicate = replicate,
      type = type,
      organism = parts[1, 3],
      oxygen = parts[1, 4],
      source = parts[1, 5],
      day = parts[1, 6],
      inoculation = parts[1, 7]
    )
  }) %>% bind_rows()

  if (nrow(meta) == 0) {
    stop("No usable targeted sample columns found.")
  }
  meta
}

pcoa_from_bray <- function(mat) {
  dist_obj <- vegan::vegdist(mat, method = "bray")
  fit <- stats::cmdscale(dist_obj, k = 2, eig = TRUE)
  eig <- fit$eig
  pos <- eig[eig > 0]
  if (length(pos) < 2) {
    stop("Not enough positive eigenvalues for 2D PCoA.")
  }
  pct <- round(100 * pos[1:2] / sum(pos), 1)
  list(coords = as.data.frame(fit$points), pct = pct)
}

load_pcoa_panel <- function(path, organism_code, include_d0 = TRUE) {
  df <- read_csv(path, show_col_types = FALSE)
  meta <- parse_sample_meta(df) %>%
    filter(organism == organism_code, oxygen == "A")

  if (include_d0) {
    meta <- meta %>% filter(type == "d0" | (type == "sup" & inoculation == "inoc"))
  } else {
    meta <- meta %>% filter(type == "sup", inoculation == "inoc", day != "d0")
  }

  if (nrow(meta) == 0) {
    stop("No aerobic samples available for ", organism_code)
  }

  meta <- meta %>%
    mutate(
      source = factor(source, levels = c("fresh", "recyc")),
      day = factor(day, levels = DAY_ORDER),
      state = if_else(type == "d0", "d0 baseline", "Aerobic"),
      fill_col = if_else(type == "d0",
                         if_else(source == "fresh", SOURCE_FILL_D0["fresh"], SOURCE_FILL_D0["recyc"]),
                         if_else(source == "fresh", SOURCE_FILL["fresh"], SOURCE_FILL["recyc"]))
    ) %>%
    arrange(type, source, day, replicate)

  cols <- meta$column
  mat <- df[, cols, drop = FALSE] %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.matrix()
  mat[is.na(mat)] <- 0
  mat <- t(log1p(mat))
  keep <- colSums(mat) > 0
  mat <- mat[, keep, drop = FALSE]
  if (ncol(mat) < 2) {
    stop("Not enough nonzero features in ", path)
  }

  ord <- pcoa_from_bray(mat)
  meta$PCoA1 <- ord$coords[, 1]
  meta$PCoA2 <- ord$coords[, 2]
  meta$pct <- ord$pct[1]
  attr(meta, "pct") <- ord$pct
  meta
}

plot_pcoa <- function(scores, title, out_base, include_d0 = TRUE) {
  pct <- attr(scores, "pct")

  plot_df <- scores
  if (!include_d0) {
    plot_df <- plot_df %>% filter(type == "sup", inoculation == "inoc", day != "d0")
  }

  plot_df <- plot_df %>%
    mutate(
      day_plot = if_else(type == "d0", "d0", as.character(day)),
      day_plot = factor(day_plot, levels = DAY_ORDER)
    )

  p <- ggplot(plot_df, aes(PCoA1, PCoA2)) +
    geom_hline(yintercept = 0, colour = "#D0D0D0", linewidth = 0.45, linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "#D0D0D0", linewidth = 0.45, linetype = "dashed") +
    geom_point(
      aes(shape = day_plot, colour = state, fill = fill_col),
      size = 4.3,
      stroke = 1.05,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = c("Aerobic" = BLUE, "d0 baseline" = GREEN),
      breaks = if (include_d0) c("Aerobic", "d0 baseline") else c("Aerobic"),
      drop = FALSE
    ) +
    scale_fill_identity(guide = "none") +
    scale_shape_manual(values = DAY_SHAPES, drop = FALSE) +
    labs(
      title = title,
      x = sprintf("PCoA1 (%.1f%% of variance)", pct[1]),
      y = sprintf("PCoA2 (%.1f%% of variance)", pct[2]),
      caption = "Open symbols = fresh; filled symbols = recycled."
    ) +
    guides(
      colour = guide_legend(title = "Color", override.aes = list(shape = 22, fill = "white", size = 4.5)),
      shape = guide_legend(title = "Day", override.aes = list(colour = "black", fill = "white", size = 4.5))
    ) +
    theme(
      legend.position = "right",
      legend.justification = c(0, 1),
      legend.box = "vertical",
      legend.box.just = "left",
      legend.box.spacing = unit(0.28, "cm"),
      plot.title = element_text(hjust = 0.5, margin = margin(b = 6)),
      plot.caption = element_text(hjust = 0.5, margin = margin(t = 7))
    )

  w <- 8.8 / 25.4
  h <- 6.4 / 25.4
  ggsave(paste0(out_base, ".png"), p, width = 9.2, height = 6.8, units = "in", dpi = 320, bg = "white")
  ggsave(paste0(out_base, ".pdf"), p, width = 9.2, height = 6.8, units = "in", device = cairo_pdf, bg = "white")
  svglite::svglite(paste0(out_base, ".svg"), width = 9.2, height = 6.8)
  print(p)
  dev.off()

  write_csv(
    plot_df %>% select(column, sample_name, replicate, type, organism, oxygen, source, day, inoculation, PCoA1, PCoA2, state),
    paste0(out_base, "_scores.csv")
  )
}

first_non_other <- function(x) {
  x <- x[!is.na(x) & x != "other"]
  if (length(x) == 0) "other" else x[1]
}

metabolite_class <- function(label) {
  s <- str_to_lower(normalize_label(label))
  class <- rep("other", length(s))

  class[str_detect(s, paste(c(
    "alanine", "arginine", "asparagine", "aspartic", "cysteine", "glutamic",
    "glutamine", "glycine", "histidine", "isoleucine", "leucine", "lysine",
    "methionine", "phenylalanine", "proline", "serine", "threonine",
    "tryptophan", "tyrosine", "valine", "ornithine", "citrulline",
    "beta-alanine", "dl-2-aminobutyric acid", "aminocaproic acid",
    "4-guanidinobutanoic acid", "n-epsilon-acetyl-l-lysine",
    "deoxycarnitine", "creatine", "creatinine", "choline", "betaine",
    "tyramine", "kynurenine", "pipecolic", "aminobutyric"
  ), collapse = "|"))] <- "amino acid/deriv."

  class[str_detect(s, paste(c(
    "adenosine", "guanosine", "cytidine", "uridine", "thymidine", "inosine",
    "deoxycytidine", "deoxyadenosine", "2-deoxyguanosine", "2'-deoxycytidine",
    "5-methylcytidine", "5-methyluridine", "5-methylthioadenosine",
    "2'-cyclic amp", "3'-cmp", "3-cmp", "2',3'-cyclic amp",
    "guanosine 3',5'-cyclic monophosphate", "adenosine 5'-monophosphoric acid",
    "uridine-5-monophosphoric acid", "thymidine 5'-monophosphoric acid",
    "cytidine-2',3'-cyclic mono-phosphoric acid", "2-deoxyadenosine",
    "2-deoxycytidine", "2-deoxyguanosine", "adenine", "guanine", "uracil",
    "thymine", "cytosine", "xanthine", "hypoxanthine"
  ), collapse = "|"))] <- "nucleic acid constituent"

  class[str_detect(s, paste(c(
    "malic acid", "fumaric acid", "pyruvic acid", "succinic acid", "citric acid",
    "gluconic acid", "2-oxovaleric acid", "alpha-ketoglutaric acid",
    "alpha-hydroxyisobutyric acid", "2-hydroxybutyric acid", "lactic acid",
    "ribose-5-phosphate", "glycerol 2-phosphoric acid",
    "alpha-glucose 1-phosphoric acid", "alpha-glucose 1-phosphate",
    "glyceric acid", "2-quinolinecarboxylic acid",
    "4-hydroxy-2-quinolinecarboxylic acid", "pterin",
    "3-dehydroshikimic acid", "4-hydroxyphenylacetic acid",
    "2-hydroxyphenylacetic acid", "4-acetamidobutanoic acid", "azelaic acid",
    "suberic acid", "maleic acid", "3-amino-5-hydroxybenzoic acid",
    "mono-methyl glutaric acid"
  ), collapse = "|"))] <- "organic acid"

  class[str_detect(s, paste(c(
    "glucose", "mannitol", "sorbitol", "xylose", "ribose", "trehalose",
    "palatinose", "myo-inositol", "inositol", "fructose", "beta-d-glucose",
    "mannose", "lactulose", "n-acetyl-galactosamine", "glycerol"
  ), collapse = "|"))] <- "carbohydrate"

  class[str_detect(s, paste(c(
    "riboflavin", "nicotinic acid", "nicotinamide", "pyridoxamine", "pyridoxal",
    "pyridoxic acid", "folic acid", "thiamine", "biotin", "pantothenic acid",
    "5-methyltetrahydrofolate", "2-hydroxypyridine", "4-pyridoxic acid"
  ), collapse = "|"))] <- "vitamin/cofactor"

  class
}

read_aerobic_inoc_results <- function(path, organism) {
  read_csv(path, show_col_types = FALSE) %>%
    mutate(label = normalize_label(label)) %>%
    mutate(
      organism = organism,
      source = str_match(comparison, "^(fresh|recyc)-A-(d\\d+)-inoc_vs_d0$")[, 2],
      day = str_match(comparison, "^(fresh|recyc)-A-(d\\d+)-inoc_vs_d0$")[, 3],
      p_value = pval
    ) %>%
    filter(!is.na(source), !is.na(day)) %>%
    mutate(
      source = factor(source, levels = c("fresh", "recyc")),
      day = factor(day, levels = c("d2", "d6", "d11", "d27")),
      comparison_label = paste(organism, if_else(source == "fresh", "Fresh", "Recycled"), as.character(day))
    )
}

make_aerobic_heatmap <- function(top_n = 40, out_tag = "top40") {
  pseudo <- read_aerobic_inoc_results(file.path(PSEU_OUT, "inoc_vs_d0_results.csv"), "Pseudomonas")
  arthro <- read_aerobic_inoc_results(file.path(ARTH_OUT, "inoc_vs_d0_results_arthro.csv"), "Arthrobacter")

  long_df <- bind_rows(pseudo, arthro) %>%
    transmute(
      label = label,
      formula = formula,
      polarity = polarity,
      class = metabolite_class(label),
      comparison_label = comparison_label,
      log2FC = log2FC,
      p_value = p_value
    )

  comparison_levels <- c(
    paste("Pseudomonas", "Fresh", c("d2", "d6", "d11", "d27")),
    paste("Pseudomonas", "Recycled", c("d2", "d6", "d11", "d27")),
    paste("Arthrobacter", "Fresh", c("d2", "d6", "d11", "d27")),
    paste("Arthrobacter", "Recycled", c("d2", "d6", "d11", "d27"))
  )

  meta <- long_df %>%
    group_by(label) %>%
    summarise(
      formula = first(na.omit(formula)),
      polarity = first(na.omit(polarity)),
      class = first_non_other(class),
      .groups = "drop"
    )

  stats_df <- long_df %>%
    group_by(comparison_label) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH"),
      star = if_else(!is.na(q_value) & q_value < 0.05, "*", "")
    ) %>%
    ungroup()

  wide <- stats_df %>%
    mutate(comparison_label = factor(comparison_label, levels = comparison_levels)) %>%
    group_by(label, comparison_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = comparison_label, values_from = log2FC)

  wide_star <- stats_df %>%
    mutate(comparison_label = factor(comparison_label, levels = comparison_levels)) %>%
    group_by(label, comparison_label) %>%
    summarise(star = first(star), .groups = "drop") %>%
    pivot_wider(names_from = comparison_label, values_from = star, values_fill = "")

  heat_tbl <- meta %>%
    left_join(wide, by = "label")
  star_tbl <- meta %>%
    left_join(wide_star, by = "label")

  heat_mat <- heat_tbl %>%
    select(all_of(comparison_levels)) %>%
    as.matrix()
  rownames(heat_mat) <- heat_tbl$label

  row_metric <- apply(heat_mat, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(-Inf)
    max(abs(x))
  })
  keep <- order(row_metric, decreasing = TRUE)
  if (is.null(top_n)) {
    keep <- keep
  } else {
    top_n <- min(top_n, length(keep))
    keep <- keep[seq_len(top_n)]
  }

  plot_mat <- heat_mat[keep, , drop = FALSE]
  plot_meta <- heat_tbl[keep, , drop = FALSE]
  star_mat <- star_tbl[keep, , drop = FALSE] %>%
    select(all_of(comparison_levels)) %>%
    as.matrix()
  plot_mat[is.na(plot_mat)] <- 0

  column_split <- factor(
    c(
      rep("Pseudomonas | Fresh", 4),
      rep("Pseudomonas | Recycled", 4),
      rep("Arthrobacter | Fresh", 4),
      rep("Arthrobacter | Recycled", 4)
    ),
    levels = c(
      "Pseudomonas | Fresh",
      "Pseudomonas | Recycled",
      "Arthrobacter | Fresh",
      "Arthrobacter | Recycled"
    )
  )
  top_info <- tibble(
    Organism = factor(
      rep(c("Pseudomonas", "Pseudomonas", "Arthrobacter", "Arthrobacter"), each = 4),
      levels = c("Pseudomonas", "Arthrobacter")
    ),
    Source = factor(
      rep(c("Fresh", "Recycled", "Fresh", "Recycled"), each = 4),
      levels = c("Fresh", "Recycled")
    ),
    Day = factor(
      rep(c("d2", "d6", "d11", "d27"), times = 4),
      levels = c("d2", "d6", "d11", "d27")
    )
  )
  top_colors <- list(
    Organism = c(Pseudomonas = "#4C78A8", Arthrobacter = "#F58518"),
    Source = c(Fresh = "#8DBCE3", Recycled = "#F0B37E"),
    Day = c(d2 = "#D0E1F2", d6 = "#9EC5E3", d11 = "#6FA8D6", d27 = "#3182BD")
  )
  top_anno <- HeatmapAnnotation(
    df = top_info,
    col = top_colors,
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
    simple_anno_size = unit(4, "mm")
  )
  row_anno <- rowAnnotation(
    Class = factor(plot_meta$class, levels = names(CLASS_PALETTE)),
    col = list(Class = CLASS_PALETTE),
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
    simple_anno_size = unit(4.2, "mm")
  )
  col_fun <- colorRamp2(c(-4, -2, 0, 2, 4), c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"))
  rownames(plot_mat) <- plot_meta$label
  colnames(plot_mat) <- rep(c("d2", "d6", "d11", "d27"), times = 4)
  rownames(star_mat) <- plot_meta$label
  colnames(star_mat) <- rep(c("d2", "d6", "d11", "d27"), times = 4)

  ht <- Heatmap(
    plot_mat,
    name = "log2FC\nrelative to d0",
    col = col_fun,
    na_col = "#F2F2F2",
    top_annotation = top_anno,
    left_annotation = row_anno,
    column_split = column_split,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    cluster_rows = TRUE,
    row_names_gp = gpar(fontsize = 7.2),
    row_names_max_width = unit(8.2, "cm"),
    column_names_gp = gpar(fontsize = 8.0),
    column_names_rot = 0,
    column_gap = unit(2.5, "mm"),
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    show_heatmap_legend = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.25),
    border = TRUE,
    column_title = "Aerobic incubation-driven metabolite change relative to d0",
    column_title_gp = gpar(fontsize = 11.5, fontface = "bold"),
    heatmap_legend_param = list(
      title = "log2FC",
      title_gp = gpar(fontsize = 8.8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7.4),
      at = c(-4, -2, 0, 2, 4),
      labels = c("-4", "-2", "0", "2", "4")
    ),
    cell_fun = function(j, i, x, y, w, h, fill) {
      star <- star_mat[i, j]
      if (!is.na(star) && nzchar(star)) {
        grid.text(star, x, y, gp = gpar(fontsize = 8.8, fontface = "bold", col = "black"))
      }
    }
  )

  class_lgd <- Legend(
    title = "Class",
    labels = names(CLASS_PALETTE),
    legend_gp = gpar(fill = CLASS_PALETTE),
    nrow = 1,
    by_row = TRUE,
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.4)
  )
  organism_lgd <- Legend(
    title = "Organism",
    labels = c("Pseudomonas", "Arthrobacter"),
    legend_gp = gpar(fill = c("#4C78A8", "#F58518")),
    nrow = 1,
    by_row = TRUE,
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.4)
  )
  source_lgd <- Legend(
    title = "Source",
    labels = c("Fresh", "Recycled"),
    legend_gp = gpar(fill = c("#8DBCE3", "#F0B37E")),
    nrow = 1,
    by_row = TRUE,
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.4)
  )
  day_lgd <- Legend(
    title = "Day",
    labels = c("d2", "d6", "d11", "d27"),
    legend_gp = gpar(fill = c("#D0E1F2", "#9EC5E3", "#6FA8D6", "#3182BD")),
    nrow = 1,
    by_row = TRUE,
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.4)
  )
  heatmap_lgd <- Legend(
    title = "log2FC",
    col_fun = col_fun,
    at = c(-4, -2, 0, 2, 4),
    labels = c("-4", "-2", "0", "2", "4"),
    direction = "horizontal",
    title_gp = gpar(fontsize = 8.8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7.4)
  )
  bottom_legends <- packLegend(
    class_lgd,
    organism_lgd,
    source_lgd,
    day_lgd,
    heatmap_lgd,
    direction = "horizontal",
    gap = unit(7, "mm"),
    max_width = unit(15.4, "in")
  )

  source_data <- heat_tbl[keep, ] %>%
    mutate(row_metric = row_metric[keep])
  write_csv(source_data, file.path(HEAT_OUT, "aerobic_t0_logfc_heatmap_source_data.csv"))

  base <- file.path(HEAT_OUT, paste0("aerobic_t0_logfc_pseudo_arthro_heatmap_", out_tag))
  save_heatmap <- function(filename, device = c("png", "pdf", "svg"), width = 16.2, height = 9.6) {
    device <- match.arg(device)
    if (device == "png") {
      ragg::agg_png(filename, width = width, height = height, units = "in", res = 320, background = "white")
    } else if (device == "pdf") {
      grDevices::cairo_pdf(filename, width = width, height = height, family = "Arial")
    } else {
      svglite::svglite(filename, width = width, height = height)
    }
    grid.newpage()
    legend_h <- grid::convertHeight(grid::grobHeight(bottom_legends@grob), "in", valueOnly = TRUE)
    legend_h <- min(max(legend_h + 0.10, 0.95), height * 0.33)
    pushViewport(viewport(
      layout = grid.layout(
        nrow = 2,
        heights = unit.c(unit(height - legend_h, "in"), unit(legend_h, "in"))
      )
    ))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    draw(
      ht,
      newpage = FALSE,
      show_heatmap_legend = FALSE,
      show_annotation_legend = FALSE
    )
    upViewport()
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid.draw(bottom_legends)
    upViewport(2)
    dev.off()
  }
  save_heatmap(paste0(base, ".png"), device = "png")
  save_heatmap(paste0(base, ".pdf"), device = "pdf")
  save_heatmap(paste0(base, ".svg"), device = "svg")
}

main <- function() {
  pseudo_scores <- load_pcoa_panel(PSEU_PATH, "Pseu", include_d0 = TRUE)
  arthro_scores <- load_pcoa_panel(ARTH_PATH, "Arth", include_d0 = TRUE)

  plot_pcoa(pseudo_scores, "Pseudomonas targeted metabolite PCoA", file.path(PSEU_OUT, "pcoa_targeted_metabolites_aerobic"))
  plot_pcoa(arthro_scores, "Arthrobacter targeted metabolite PCoA", file.path(ARTH_OUT, "pcoa_targeted_metabolites_aerobic"))
  make_aerobic_heatmap(top_n = 40, out_tag = "top40")
  make_aerobic_heatmap(top_n = NULL, out_tag = "all")
}

main()
