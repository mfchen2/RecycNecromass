suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(broom)
})

set.seed(20260505)

out_dir <- "outputs/manuscript_story_analysis/mag_functional_guild_redundancy"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_div <- function(x, y) ifelse(y == 0, NA_real_, x / y)
safe_log2_ratio <- function(x, y, pseudocount = 0.5) log2((x + pseudocount) / (y + pseudocount))

parse_group <- function(df, group_col = "group") {
  df %>%
    mutate(
      source_type = case_when(
        str_detect(.data[[group_col]], "^fresh") ~ "Fresh",
        str_detect(.data[[group_col]], "^recyc") ~ "Recyc",
        TRUE ~ NA_character_
      ),
      source_carbon = case_when(
        str_detect(.data[[group_col]], "pseudo") ~ "Pseudomonas",
        str_detect(.data[[group_col]], "arthro") ~ "Arthrobacter",
        TRUE ~ NA_character_
      ),
      oxygen = case_when(
        str_detect(.data[[group_col]], "_aer$") ~ "Aerobic",
        str_detect(.data[[group_col]], "_ana$") ~ "Anaerobic",
        TRUE ~ NA_character_
      )
    )
}

mg_to_time <- function(mg_num) {
  case_when(
    mg_num >= 1 & mg_num <= 8 ~ "T2",
    mg_num >= 9 & mg_num <= 16 ~ "T3",
    mg_num >= 17 & mg_num <= 24 ~ "T5",
    TRUE ~ NA_character_
  )
}

time_to_day <- function(tp) as.numeric(recode(tp, T2 = 6, T3 = 11, T5 = 27))

clean_rank <- function(x, prefix) {
  x %>%
    str_remove(paste0("^", prefix, "__")) %>%
    na_if("") %>%
    replace_na("Unclassified")
}

shannon <- function(x) {
  x <- table(x)
  p <- as.numeric(x) / sum(x)
  -sum(p * log(p))
}

raw_matrix_to_long <- function(path, group_name) {
  raw <- read_tsv(path, show_col_types = FALSE)
  parts <- str_split_fixed(raw$KO_gene_definition, fixed("|"), 3)

  raw %>%
    mutate(
      gene_name = parts[, 1],
      ko = parts[, 2],
      ko_definition = parts[, 3]
    ) %>%
    pivot_longer(
      cols = -c(KO_gene_definition, gene_name, ko, ko_definition),
      names_to = "sample_bin",
      values_to = "present"
    ) %>%
    filter(present > 0) %>%
    mutate(
      sample_name = str_replace(sample_bin, "__bin.*$", ""),
      group = group_name
    ) %>%
    select(KO_gene_definition, gene_name, ko, ko_definition, sample_bin, present, sample_name, group)
}

denitrification_kos <- tribble(
  ~module, ~step, ~gene_symbol, ~ko, ~description,
  "denitrification", "nitrate_to_nitrite", "narG", "K00370", "membrane nitrate reductase alpha / nitrate-nitrite oxidoreductase alpha",
  "denitrification", "nitrate_to_nitrite", "napA_like", "K02567", "periplasmic/cytochrome nitrate reductase catalytic subunit proxy",
  "denitrification", "nitrite_to_NO", "nirK_nirS_like", "K00368", "NO-forming nitrite reductase",
  "denitrification", "nitrite_to_NO", "nirK_nirS_like", "K15864", "NO-forming nitrite reductase / hydroxylamine reductase",
  "denitrification", "NO_to_N2O", "norB", "K04561", "nitric oxide reductase subunit B",
  "denitrification", "N2O_to_N2", "nosZ", "K00376", "nitrous-oxide reductase"
)

dnra_kos <- tribble(
  ~module, ~step, ~gene_symbol, ~ko, ~description,
  "DNRA", "nitrite_to_ammonium", "nrfA", "K03385", "cytochrome c nitrite reductase",
  "DNRA", "nitrite_to_ammonium", "nrfH_like", "K15876", "cytochrome c nitrite reductase small subunit",
  "DNRA", "accessory", "nrfG", "K04018", "formate-dependent nitrite reductase complex subunit NrfG"
)

fermentation_kos <- tribble(
  ~module, ~step, ~gene_symbol, ~ko, ~description,
  "fermentation", "alcohol", "adh", "K00001", "alcohol dehydrogenase",
  "fermentation", "lactate", "ldh", "K00016", "L-lactate dehydrogenase",
  "fermentation", "lactate", "lldD", "K00101", "cytochrome L-lactate dehydrogenase",
  "fermentation", "lactate", "dld", "K00102", "cytochrome D-lactate dehydrogenase",
  "fermentation", "formate", "fdh", "K00123", "formate dehydrogenase major subunit",
  "fermentation", "formate", "fdh", "K00124", "formate dehydrogenase Fe-S subunit",
  "fermentation", "formate", "fdh", "K00126", "formate dehydrogenase delta subunit",
  "fermentation", "formate", "fdh", "K00127", "formate dehydrogenase gamma subunit",
  "fermentation", "pyruvate_formate", "pflB", "K00656", "pyruvate formate-lyase",
  "fermentation", "acetate", "ackA", "K00925", "acetate kinase",
  "fermentation", "pyruvate_formate", "pflA", "K04069", "pyruvate formate-lyase activating enzyme",
  "fermentation", "ethanol", "adhE", "K04072", "acetaldehyde/alcohol dehydrogenase",
  "fermentation", "formate_transport", "focA", "K06212", "formate transporter",
  "fermentation", "formate_transport", "oxlT", "K08177", "oxalate/formate antiporter",
  "fermentation", "formate", "fdnH", "K08349", "formate dehydrogenase-N beta subunit",
  "fermentation", "formate", "fdnI", "K08350", "formate dehydrogenase-N gamma subunit"
)

pathway_ref <- read_tsv("Summarized KoFam/pathway_ko_reference.tsv", show_col_types = FALSE) %>%
  distinct(pathway_category, pathway_name, ko)

cazyme_peptidase_kos <- pathway_ref %>%
  mutate(
    module = case_when(
      pathway_category %in% c("Microbial necromass", "Amino sugar metabolism", "Plant polymer degradation", "EPS_biofilm") ~ "CAZyme_like_carbon_hydrolase",
      pathway_category == "Protein necromass" ~ "peptidase_peptide_uptake",
      TRUE ~ NA_character_
    ),
    step = pathway_name,
    gene_symbol = pathway_name,
    description = paste(pathway_category, pathway_name, sep = " | ")
  ) %>%
  filter(!is.na(module)) %>%
  select(module, step, gene_symbol, ko, description)

functional_reference <- bind_rows(
  denitrification_kos,
  dnra_kos,
  fermentation_kos,
  cazyme_peptidase_kos
) %>%
  distinct(module, step, gene_symbol, ko, description)

write_csv(functional_reference, file.path(out_dir, "functional_ko_reference_used.csv"))

ko_long_pre <- read_tsv(
  "Summarized KoFam/pathway_summary_outputs/ko_summary/all_groups_ko_presence_long.tsv",
  show_col_types = FALSE
)

raw_augments <- c(
  fresh_arthro_ana = "Summarized KoFam/fresh_arthro_ana_kofam_matrix.tsv"
)

raw_augmented <- bind_rows(lapply(names(raw_augments), function(group_name) {
  path <- raw_augments[[group_name]]
  if (!file.exists(path)) return(tibble())
  raw_matrix_to_long(path, group_name)
}))

ko_long <- bind_rows(
  ko_long_pre %>% filter(group != "fresh_arthro_ana"),
  raw_augmented
) %>%
  filter(str_detect(sample_bin, "MG[0-9]+__bin\\.")) %>%
  mutate(
    assembler = str_extract(sample_bin, "^[^-]+"),
    mg_num = as.integer(str_match(sample_bin, "MG([0-9]+)")[, 2]),
    timepoint = mg_to_time(mg_num),
    day_num = time_to_day(timepoint),
    mag_id = sample_bin
  ) %>%
  parse_group("group")

tax <- read_excel("Summarized_taxa/Summary SK.xlsx") %>%
  mutate(
    group = paste(Necromass_type, necromass_organisms, oxygen, sep = "_"),
    mag_id = paste0(sample_ID, "__", bins),
    assembler = str_extract(sample_ID, "^[^-]+"),
    mg_num = as.integer(str_match(sample_ID, "MG([0-9]+)")[, 2]),
    timepoint = mg_to_time(mg_num),
    day_num = time_to_day(timepoint),
    phylum = clean_rank(Phylum, "p"),
    class = clean_rank(Class, "c"),
    order = clean_rank(Order, "o"),
    family = clean_rank(Family, "f"),
    genus = clean_rank(Genus, "g"),
    species = clean_rank(Species, "s")
  ) %>%
  parse_group("group") %>%
  select(mag_id, group, source_type, source_carbon, oxygen, assembler, mg_num, timepoint, day_num,
         phylum, class, order, family, genus, species)

mag_meta <- tax %>%
  semi_join(ko_long %>% distinct(mag_id), by = "mag_id")

availability <- tibble(
  input = c(
    "MAG-level KOfam presence",
    "MAG taxonomy",
    "fresh_arthro_ana KO matrix",
    "fresh_pseudo_aer taxonomy-linked KO matrix",
    "MAG abundance/coverage across samples",
    "nosZ clade I/II annotation"
  ),
  status = c(
    paste0("available for ", n_distinct(ko_long$mag_id), " taxonomy-linkable MAGs"),
    paste0("available for ", n_distinct(tax$mag_id), " MAGs"),
    if_else(nrow(raw_augmented) > 0, "augmented from raw matrix", "missing"),
    "not taxonomy-linkable in current KO long table; sample_bin IDs are Refined_*",
    "not found",
    "not found; KOfam resolves nosZ as K00376 only"
  ),
  implication = c(
    "MAG guild/completeness tests can be run on taxonomy-linkable MAGs.",
    "Taxonomic redundancy tests can be run.",
    "Fresh Arthrobacter anaerobic can be included.",
    "Fresh Pseudomonas aerobic excluded from MAG-level taxonomic redundancy unless bin ID map is provided.",
    "Bloom tracking is limited to recovered MAG/bin counts by group/timepoint, not true relative abundance.",
    "nir/nosZ ratios can be calculated, but nosZ clade I versus II cannot be separated."
  )
)
write_csv(availability, file.path(out_dir, "mag_analysis_input_availability.csv"))

mag_ko <- ko_long %>%
  semi_join(tax, by = "mag_id") %>%
  group_by(mag_id, group, source_type, source_carbon, oxygen, assembler, mg_num, timepoint, day_num, ko) %>%
  summarise(
    ko_definition = first(ko_definition),
    n_gene_hits = n(),
    .groups = "drop"
  )

mag_totals <- mag_ko %>%
  group_by(mag_id) %>%
  summarise(total_ko_gene_hits = sum(n_gene_hits), n_distinct_kos = n_distinct(ko), .groups = "drop")

mag_function_long <- mag_ko %>%
  inner_join(functional_reference, by = "ko", relationship = "many-to-many") %>%
  left_join(mag_meta, by = c("mag_id", "group", "source_type", "source_carbon", "oxygen", "assembler", "mg_num", "timepoint", "day_num")) %>%
  left_join(mag_totals, by = "mag_id") %>%
  mutate(percent_of_ko_hits = 100 * n_gene_hits / total_ko_gene_hits)

write_csv(mag_function_long, file.path(out_dir, "mag_functional_hits_long.csv"))

has_ko <- function(kos, set) as.integer(any(kos %in% set))
count_ko <- function(kos, set) sum(kos %in% set)

mag_flags <- mag_ko %>%
  group_by(mag_id, group, source_type, source_carbon, oxygen, assembler, mg_num, timepoint, day_num) %>%
  summarise(kos = list(unique(ko)), .groups = "drop") %>%
  left_join(mag_meta, by = c("mag_id", "group", "source_type", "source_carbon", "oxygen", "assembler", "mg_num", "timepoint", "day_num")) %>%
  left_join(mag_totals, by = "mag_id") %>%
  rowwise() %>%
  mutate(
    has_narG = has_ko(kos, "K00370"),
    has_napA_like = has_ko(kos, "K02567"),
    has_nitrate_reduction = as.integer(has_narG == 1 || has_napA_like == 1),
    has_nirK_nirS_like = has_ko(kos, c("K00368", "K15864")),
    has_norB = has_ko(kos, "K04561"),
    has_nosZ = has_ko(kos, "K00376"),
    has_nrfA = has_ko(kos, "K03385"),
    has_dnra_module = has_ko(kos, c("K03385", "K15876")),
    has_fermentation = has_ko(kos, fermentation_kos$ko),
    n_fermentation_kos = count_ko(kos, fermentation_kos$ko),
    has_cazyme_like = has_ko(kos, cazyme_peptidase_kos$ko[cazyme_peptidase_kos$module == "CAZyme_like_carbon_hydrolase"]),
    n_cazyme_like_kos = count_ko(kos, cazyme_peptidase_kos$ko[cazyme_peptidase_kos$module == "CAZyme_like_carbon_hydrolase"]),
    has_peptidase_or_peptide = has_ko(kos, cazyme_peptidase_kos$ko[cazyme_peptidase_kos$module == "peptidase_peptide_uptake"]),
    n_peptidase_or_peptide_kos = count_ko(kos, cazyme_peptidase_kos$ko[cazyme_peptidase_kos$module == "peptidase_peptide_uptake"]),
    denitrification_steps_present = has_nitrate_reduction + has_nirK_nirS_like + has_norB + has_nosZ,
    denitrification_completeness = denitrification_steps_present / 4,
    full_denitrifier = as.integer(denitrification_steps_present == 4),
    partial_denitrifier_no_nosZ = as.integer(has_nitrate_reduction == 1 && has_nirK_nirS_like == 1 && has_norB == 1 && has_nosZ == 0),
    nitrite_reducer_no_nosZ = as.integer(has_nirK_nirS_like == 1 && has_nosZ == 0),
    n_nir_like_kos = count_ko(kos, c("K00368", "K15864")),
    n_nosZ_kos = count_ko(kos, "K00376")
  ) %>%
  ungroup() %>%
  select(-kos)

write_csv(mag_flags, file.path(out_dir, "mag_functional_guild_flags_and_completeness.csv"))

sample_summary <- mag_flags %>%
  group_by(group, source_type, source_carbon, oxygen, timepoint, day_num) %>%
  summarise(
    n_mags = n(),
    n_full_denitrifiers = sum(full_denitrifier),
    frac_full_denitrifiers = n_full_denitrifiers / n_mags,
    n_partial_denit_no_nosZ = sum(partial_denitrifier_no_nosZ),
    frac_partial_denit_no_nosZ = n_partial_denit_no_nosZ / n_mags,
    n_nitrite_reducer_no_nosZ = sum(nitrite_reducer_no_nosZ),
    frac_nitrite_reducer_no_nosZ = n_nitrite_reducer_no_nosZ / n_mags,
    n_dnra_nrfA = sum(has_nrfA),
    frac_dnra_nrfA = n_dnra_nrfA / n_mags,
    n_fermentation = sum(has_fermentation),
    frac_fermentation = n_fermentation / n_mags,
    n_cazyme_like = sum(has_cazyme_like),
    frac_cazyme_like = n_cazyme_like / n_mags,
    n_peptidase_or_peptide = sum(has_peptidase_or_peptide),
    frac_peptidase_or_peptide = n_peptidase_or_peptide / n_mags,
    mean_denitrification_completeness = mean(denitrification_completeness),
    n_nir_like_mags = sum(has_nirK_nirS_like),
    n_nosZ_mags = sum(has_nosZ),
    nir_to_nosZ_mag_ratio = safe_div(n_nir_like_mags, n_nosZ_mags),
    .groups = "drop"
  )
write_csv(sample_summary, file.path(out_dir, "sample_mag_guild_summary.csv"))

treatment_summary <- mag_flags %>%
  group_by(source_type, source_carbon, oxygen) %>%
  summarise(
    n_mags = n(),
    n_samples = n_distinct(paste(group, timepoint)),
    n_full_denitrifiers = sum(full_denitrifier),
    frac_full_denitrifiers = n_full_denitrifiers / n_mags,
    n_partial_denit_no_nosZ = sum(partial_denitrifier_no_nosZ),
    frac_partial_denit_no_nosZ = n_partial_denit_no_nosZ / n_mags,
    n_dnra_nrfA = sum(has_nrfA),
    frac_dnra_nrfA = n_dnra_nrfA / n_mags,
    n_fermentation = sum(has_fermentation),
    frac_fermentation = n_fermentation / n_mags,
    n_cazyme_like = sum(has_cazyme_like),
    frac_cazyme_like = n_cazyme_like / n_mags,
    n_peptidase_or_peptide = sum(has_peptidase_or_peptide),
    frac_peptidase_or_peptide = n_peptidase_or_peptide / n_mags,
    mean_denitrification_completeness = mean(denitrification_completeness),
    n_nir_like_mags = sum(has_nirK_nirS_like),
    n_nosZ_mags = sum(has_nosZ),
    nir_to_nosZ_mag_ratio = safe_div(n_nir_like_mags, n_nosZ_mags),
    .groups = "drop"
  )
write_csv(treatment_summary, file.path(out_dir, "treatment_mag_guild_summary.csv"))

guild_membership_long <- mag_flags %>%
  transmute(
    mag_id, source_type, source_carbon, oxygen, timepoint, day_num,
    phylum, class, order, family, genus, species,
    `full denitrification` = full_denitrifier == 1,
    `partial denitrification no nosZ` = partial_denitrifier_no_nosZ == 1,
    `nitrite reducer no nosZ` = nitrite_reducer_no_nosZ == 1,
    `DNRA nrfA` = has_nrfA == 1,
    fermentation = has_fermentation == 1,
    `CAZyme-like` = has_cazyme_like == 1,
    `peptidase/peptide uptake` = has_peptidase_or_peptide == 1
  ) %>%
  pivot_longer(
    cols = c(`full denitrification`, `partial denitrification no nosZ`, `nitrite reducer no nosZ`,
             `DNRA nrfA`, fermentation, `CAZyme-like`, `peptidase/peptide uptake`),
    names_to = "guild",
    values_to = "in_guild"
  ) %>%
  filter(in_guild)

tax_diversity <- guild_membership_long %>%
  group_by(guild) %>%
  summarise(
    n_mags = n(),
    n_phyla = n_distinct(phylum),
    n_classes = n_distinct(class),
    n_orders = n_distinct(order),
    n_families = n_distinct(family),
    n_genera = n_distinct(genus),
    genus_shannon = shannon(genus),
    effective_genera = exp(genus_shannon),
    top_genera = paste(names(sort(table(genus), decreasing = TRUE))[1:min(8, n_distinct(genus))], collapse = "; "),
    .groups = "drop"
  )

null_tax_diversity <- map_dfr(unique(guild_membership_long$guild), function(g) {
  ids <- guild_membership_long %>% filter(guild == g) %>% pull(mag_id) %>% unique()
  n <- length(ids)
  pool <- mag_flags
  if (n < 2 || n > nrow(pool)) return(tibble())
  map_dfr(seq_len(999), function(i) {
    sampled <- sample(seq_len(nrow(pool)), n, replace = FALSE)
    tibble(
      guild = g,
      perm = i,
      null_n_genera = n_distinct(pool$genus[sampled]),
      null_genus_shannon = shannon(pool$genus[sampled]),
      null_effective_genera = exp(null_genus_shannon)
    )
  })
})

tax_redundancy_tests <- tax_diversity %>%
  left_join(
    null_tax_diversity %>%
      group_by(guild) %>%
      summarise(
        null_mean_n_genera = mean(null_n_genera),
        null_mean_effective_genera = mean(null_effective_genera),
        p_more_genus_rich_than_random = (sum(null_n_genera >= tax_diversity$n_genera[match(first(guild), tax_diversity$guild)]) + 1) / (n() + 1),
        p_less_genus_rich_than_random = (sum(null_n_genera <= tax_diversity$n_genera[match(first(guild), tax_diversity$guild)]) + 1) / (n() + 1),
        .groups = "drop"
      ),
    by = "guild"
  ) %>%
  mutate(
    broad_taxonomic_redundancy = n_genera >= 3 & n_phyla >= 2,
    interpretation = case_when(
      broad_taxonomic_redundancy & p_less_genus_rich_than_random > 0.05 ~ "taxonomically broad guild",
      broad_taxonomic_redundancy ~ "broad but less rich than random MAG draw",
      TRUE ~ "taxonomically narrow or rare guild"
    )
  )

write_csv(guild_membership_long, file.path(out_dir, "mag_guild_membership_long.csv"))
write_csv(tax_diversity, file.path(out_dir, "guild_taxonomic_diversity_observed.csv"))
write_csv(null_tax_diversity, file.path(out_dir, "guild_taxonomic_diversity_null_permutations.csv"))
write_csv(tax_redundancy_tests, file.path(out_dir, "guild_taxonomic_redundancy_tests.csv"))

fresh_recyc_ana_fisher <- mag_flags %>%
  filter(oxygen == "Anaerobic", source_type %in% c("Fresh", "Recyc")) %>%
  select(mag_id, source_type, source_carbon, full_denitrifier, partial_denitrifier_no_nosZ,
         nitrite_reducer_no_nosZ, has_nrfA, has_fermentation, has_cazyme_like, has_peptidase_or_peptide) %>%
  pivot_longer(
    cols = c(full_denitrifier, partial_denitrifier_no_nosZ, nitrite_reducer_no_nosZ,
             has_nrfA, has_fermentation, has_cazyme_like, has_peptidase_or_peptide),
    names_to = "feature",
    values_to = "present"
  ) %>%
  group_by(feature) %>%
  group_modify(~{
    tab <- table(.x$source_type, .x$present > 0)
    if (nrow(tab) < 2 || ncol(tab) < 2) {
      return(tibble(p.value = NA_real_, odds_ratio = NA_real_))
    }
    ft <- fisher.test(tab)
    tibble(p.value = ft$p.value, odds_ratio = unname(ft$estimate))
  }) %>%
  ungroup() %>%
  mutate(q.value = p.adjust(p.value, method = "BH"))
write_csv(fresh_recyc_ana_fisher, file.path(out_dir, "fresh_vs_recycled_anaerobic_mag_guild_fisher_tests.csv"))

fresh_recyc_ana_sample_tests <- sample_summary %>%
  filter(oxygen == "Anaerobic") %>%
  select(source_type, source_carbon, timepoint, frac_full_denitrifiers, frac_partial_denit_no_nosZ,
         frac_nitrite_reducer_no_nosZ, frac_dnra_nrfA, mean_denitrification_completeness,
         nir_to_nosZ_mag_ratio) %>%
  pivot_longer(
    cols = c(frac_full_denitrifiers, frac_partial_denit_no_nosZ, frac_nitrite_reducer_no_nosZ,
             frac_dnra_nrfA, mean_denitrification_completeness, nir_to_nosZ_mag_ratio),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric) %>%
  group_modify(~{
    df <- .x %>% filter(!is.na(value))
    if (n_distinct(df$source_type) < 2 || sd(df$value) == 0) {
      return(tibble(p.value = NA_real_, fresh_mean = mean(df$value[df$source_type == "Fresh"], na.rm = TRUE),
                    recyc_mean = mean(df$value[df$source_type == "Recyc"], na.rm = TRUE)))
    }
    tibble(
      p.value = wilcox.test(value ~ source_type, data = df, exact = FALSE)$p.value,
      fresh_mean = mean(df$value[df$source_type == "Fresh"], na.rm = TRUE),
      recyc_mean = mean(df$value[df$source_type == "Recyc"], na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  mutate(q.value = p.adjust(p.value, method = "BH"),
         fresh_minus_recyc = fresh_mean - recyc_mean)
write_csv(fresh_recyc_ana_sample_tests, file.path(out_dir, "fresh_vs_recycled_anaerobic_sample_level_tests.csv"))

nir_nos_treatment <- treatment_summary %>%
  mutate(
    nosZ_clade_I_mags = NA_integer_,
    nosZ_clade_II_mags = NA_integer_,
    nosZ_clade_note = "KOfam K00376 does not distinguish nosZ clade I versus clade II"
  ) %>%
  select(source_type, source_carbon, oxygen, n_mags, n_nir_like_mags, n_nosZ_mags,
         nir_to_nosZ_mag_ratio, nosZ_clade_I_mags, nosZ_clade_II_mags, nosZ_clade_note)
write_csv(nir_nos_treatment, file.path(out_dir, "nir_like_to_nosZ_ratios_by_treatment.csv"))

phys_fcm <- read_csv("outputs/manuscript_story_analysis/physiology_fcm_gctcd_ph7/fcm_log10_cell_count_summary_ph7.csv", show_col_types = FALSE) %>%
  filter(timepoint %in% c("T2", "T3", "T5")) %>%
  transmute(source_type, source_carbon, oxygen, timepoint,
            log10_cells_mean = mean, log10_cells_se = se)

phys_co2 <- read_csv("outputs/manuscript_story_analysis/physiology_fcm_gctcd_ph7/gctcd_co2_summary_ph7.csv", show_col_types = FALSE) %>%
  filter(timepoint %in% c("T2", "T3", "T5")) %>%
  transmute(source_type, source_carbon, oxygen, timepoint,
            co2_rate_since_start_mean, co2_rate_interval_mean)

ic <- read_csv("outputs/manuscript_story_analysis/ic_nitrate_nitrite/ic_summary_by_condition_timepoint.csv", show_col_types = FALSE) %>%
  filter(timepoint %in% c("T2", "T3", "T5"), pH == 7, source_type %in% c("Fresh", "Recyc")) %>%
  transmute(source_type, source_carbon, oxygen, timepoint,
            nitrate_mM_mean, nitrite_mM_mean, nox_mM_mean)

guild_phys <- sample_summary %>%
  left_join(phys_fcm, by = c("source_type", "source_carbon", "oxygen", "timepoint")) %>%
  left_join(phys_co2, by = c("source_type", "source_carbon", "oxygen", "timepoint")) %>%
  left_join(ic, by = c("source_type", "source_carbon", "oxygen", "timepoint"))
write_csv(guild_phys, file.path(out_dir, "sample_guild_summary_linked_to_physiology.csv"))

phys_correlations <- guild_phys %>%
  select(frac_full_denitrifiers, frac_partial_denit_no_nosZ, frac_nitrite_reducer_no_nosZ,
         frac_dnra_nrfA, frac_fermentation, frac_cazyme_like, frac_peptidase_or_peptide,
         mean_denitrification_completeness, nir_to_nosZ_mag_ratio,
         log10_cells_mean, co2_rate_since_start_mean, co2_rate_interval_mean,
         nitrate_mM_mean, nitrite_mM_mean, nox_mM_mean) %>%
  pivot_longer(
    cols = c(frac_full_denitrifiers, frac_partial_denit_no_nosZ, frac_nitrite_reducer_no_nosZ,
             frac_dnra_nrfA, frac_fermentation, frac_cazyme_like, frac_peptidase_or_peptide,
             mean_denitrification_completeness, nir_to_nosZ_mag_ratio),
    names_to = "guild_metric",
    values_to = "guild_value"
  ) %>%
  pivot_longer(
    cols = c(log10_cells_mean, co2_rate_since_start_mean, co2_rate_interval_mean,
             nitrate_mM_mean, nitrite_mM_mean, nox_mM_mean),
    names_to = "physiology_metric",
    values_to = "physiology_value"
  ) %>%
  filter(!is.na(guild_value), !is.na(physiology_value)) %>%
  group_by(guild_metric, physiology_metric) %>%
  summarise(
    n = n(),
    spearman_rho = suppressWarnings(cor(guild_value, physiology_value, method = "spearman")),
    p.value = if_else(n >= 4 && sd(guild_value) > 0 && sd(physiology_value) > 0,
                      suppressWarnings(cor.test(guild_value, physiology_value, method = "spearman", exact = FALSE)$p.value),
                      NA_real_),
    .groups = "drop"
  ) %>%
  group_by(physiology_metric) %>%
  mutate(q.value = p.adjust(p.value, method = "BH")) %>%
  ungroup() %>%
  arrange(p.value)
write_csv(phys_correlations, file.path(out_dir, "guild_physiology_spearman_correlations.csv"))

bloom_candidates <- mag_flags %>%
  filter(oxygen == "Anaerobic") %>%
  group_by(source_type, source_carbon, timepoint, phylum, class, order, family, genus) %>%
  summarise(
    n_recovered_mags = n(),
    n_full_denitrifiers = sum(full_denitrifier),
    n_partial_denit_no_nosZ = sum(partial_denitrifier_no_nosZ),
    n_dnra_nrfA = sum(has_nrfA),
    n_fermentation = sum(has_fermentation),
    n_cazyme_like = sum(has_cazyme_like),
    n_peptidase_or_peptide = sum(has_peptidase_or_peptide),
    example_mags = paste(head(mag_id, 5), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(source_carbon, timepoint, desc(n_recovered_mags))
write_csv(bloom_candidates, file.path(out_dir, "anaerobic_recovered_mag_bloom_candidates_by_taxon.csv"))

p_guild <- treatment_summary %>%
  select(source_type, source_carbon, oxygen, frac_full_denitrifiers, frac_partial_denit_no_nosZ,
         frac_dnra_nrfA, frac_fermentation, frac_cazyme_like, frac_peptidase_or_peptide) %>%
  pivot_longer(starts_with("frac_"), names_to = "guild_metric", values_to = "fraction_mags") %>%
  mutate(
    guild_metric = recode(
      guild_metric,
      frac_full_denitrifiers = "full denitrification",
      frac_partial_denit_no_nosZ = "partial denit no nosZ",
      frac_dnra_nrfA = "DNRA nrfA",
      frac_fermentation = "fermentation",
      frac_cazyme_like = "CAZyme-like",
      frac_peptidase_or_peptide = "peptidase/peptide"
    ),
    treatment = paste(source_type, source_carbon, oxygen, sep = " | ")
  ) %>%
  ggplot(aes(x = treatment, y = fraction_mags, fill = oxygen)) +
  geom_col(width = 0.72) +
  facet_wrap(~guild_metric, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
  labs(x = NULL, y = "Fraction of recovered MAGs", fill = "Oxygen",
       title = "MAG functional guild membership by treatment")
ggsave(file.path(out_dir, "01_mag_guild_fraction_by_treatment.png"), p_guild, width = 11, height = 8, dpi = 300)
ggsave(file.path(out_dir, "01_mag_guild_fraction_by_treatment.pdf"), p_guild, width = 11, height = 8)

p_completeness <- sample_summary %>%
  filter(oxygen == "Anaerobic") %>%
  ggplot(aes(x = timepoint, y = mean_denitrification_completeness, color = source_type, group = source_type)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~source_carbon) +
  theme_bw(base_size = 10) +
  labs(x = "Timepoint", y = "Mean MAG denitrification completeness", color = "Necromass",
       title = "Anaerobic denitrification completeness in recovered MAGs")
ggsave(file.path(out_dir, "02_anaerobic_denitrification_completeness_timecourse.png"), p_completeness, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(out_dir, "02_anaerobic_denitrification_completeness_timecourse.pdf"), p_completeness, width = 7.5, height = 4.5)

p_nir_nos <- treatment_summary %>%
  ggplot(aes(x = paste(source_type, source_carbon, sep = " | "), y = nir_to_nosZ_mag_ratio, fill = oxygen)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = NULL, y = "MAG count ratio: nir-like / nosZ", fill = "Oxygen",
       title = "nirK/nirS-like to nosZ ratio by treatment")
ggsave(file.path(out_dir, "03_nir_like_to_nosZ_ratio_by_treatment.png"), p_nir_nos, width = 8, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "03_nir_like_to_nosZ_ratio_by_treatment.pdf"), p_nir_nos, width = 8, height = 4.8)

message("Wrote MAG guild outputs to: ", out_dir)
