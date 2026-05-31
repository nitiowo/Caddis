library(tidyverse)
library(Biostrings)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

aln_file <- file.path(root, "nav_selection/outputs/alignments/nav_protein_alignment_trimmed.fasta")
rep_meta_file <- file.path(root, "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_metadata.csv")
rnaseq_file <- file.path(root, "data/metadata/MASTER_caddis_rnaseq_samples.csv")
exclude_file <- file.path(root, "nav_selection/config/nav_primary_exclude_samples.csv")
site_context_file <- file.path(root, "nav_selection/outputs/selection/nav_site_context.csv")

sel_dir <- file.path(root, "nav_selection/outputs/selection")
dir.create(sel_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Functions ----

norm_sample_id <- function(value) {
  sample_num <- str_match(str_trim(as.character(value)), "(?i)CF[-_ ]*0*([0-9]+)")[, 2]
  ifelse(is.na(sample_num), NA_character_, sprintf("CF%03d", as.integer(sample_num)))
}

as_bool <- function(value) {
  toupper(as.character(value)) == "TRUE"
}

calc_dist <- function(seq_one, seq_two) {
  aa <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
  called <- seq_one %in% aa & seq_two %in% aa
  n_compared <- sum(called)

  if (n_compared == 0) {
    return(tibble(n_compared = 0L, p_distance = NA_real_))
  }

  tibble(
    n_compared = n_compared,
    p_distance = sum(seq_one[called] != seq_two[called]) / n_compared
  )
}

read_burden <- function(path, source_name) {
  if (!file.exists(path)) {
    return(tibble())
  }
  read_csv(path, show_col_types = FALSE) %>%
    mutate(source = source_name)
}

# ---- Sample Metadata ----

aln <- readAAStringSet(aln_file)
names(aln) <- sub("^(\\S+).*", "\\1", names(aln))
aln_mat <- as.matrix(aln)

rep_meta <- read_csv(rep_meta_file, show_col_types = FALSE)

rnaseq <- read_csv(rnaseq_file, show_col_types = FALSE) %>%
  mutate(
    sample_id = norm_sample_id(SAMPLE),
    species_label = paste(str_trim(`Genbank genus`), str_trim(`Genbank species`)),
    focal = str_trim(`Genbank species`) %in% c("flavastellus", "externus")
  ) %>%
  select(sample_id, species_label, focal, Location)

if (file.exists(exclude_file)) {
  exclude_tbl <- read_csv(exclude_file, show_col_types = FALSE) %>%
    mutate(sample_id = norm_sample_id(sample_id))
} else {
  exclude_tbl <- tibble(sample_id = character(), reason = character())
}

sample_ann <- rep_meta %>%
  select(sample_id, seq_id, orf_type) %>%
  left_join(rnaseq, by = "sample_id") %>%
  filter(seq_id %in% rownames(aln_mat)) %>%
  mutate(
    focal = as_bool(focal),
    genus = word(species_label, 1),
    excluded_primary = sample_id %in% exclude_tbl$sample_id
  ) %>%
  arrange(sample_id)

audit_tbl <- sample_ann %>%
  mutate(
    called_aa = rowSums(aln_mat[seq_id, , drop = FALSE] %in% strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]),
    alignment_width = ncol(aln_mat),
    occupancy = called_aa / alignment_width
  ) %>%
  left_join(exclude_tbl, by = "sample_id")

write_csv(audit_tbl, file.path(sel_dir, "nav_primary_sequence_audit.csv"))

# ---- Pairwise protein distances ----

primary_ann <- sample_ann %>%
  filter(!excluded_primary)

if (nrow(primary_ann) >= 2) {
  pair_index <- combn(seq_len(nrow(primary_ann)), 2)

  pair_tbl <- map_dfr(seq_len(ncol(pair_index)), function(pair_num) {
    row_a <- pair_index[1, pair_num]
    row_b <- pair_index[2, pair_num]
    meta_a <- primary_ann[row_a, ]
    meta_b <- primary_ann[row_b, ]
    dist_tbl <- calc_dist(aln_mat[meta_a$seq_id, ], aln_mat[meta_b$seq_id, ])

    tibble(
      sample_a = meta_a$sample_id,
      sample_b = meta_b$sample_id,
      species_a = meta_a$species_label,
      species_b = meta_b$species_label,
      focal_a = meta_a$focal,
      focal_b = meta_b$focal,
      genus_a = meta_a$genus,
      genus_b = meta_b$genus,
      n_compared = dist_tbl$n_compared,
      p_distance = dist_tbl$p_distance
    )
  })

  pair_classes <- bind_rows(
    pair_tbl %>% filter(focal_a & focal_b) %>% mutate(pair_class = "Within LiFl"),
    pair_tbl %>% filter(focal_a != focal_b, genus_a == "Limnephilus", genus_b == "Limnephilus") %>% mutate(pair_class = "LiFl to other Limnephilus"),
    pair_tbl %>% filter(focal_a != focal_b) %>% mutate(pair_class = "LiFl to all non-LiFl"),
    pair_tbl %>% filter(!focal_a, !focal_b, genus_a == "Limnephilus", genus_b == "Limnephilus") %>% mutate(pair_class = "Non-LiFl Limnephilus"),
    pair_tbl %>% filter(!focal_a, !focal_b) %>% mutate(pair_class = "All non-LiFl")
  )
} else {
  pair_classes <- tibble()
}

pair_summary <- pair_classes %>%
  filter(!is.na(p_distance)) %>%
  group_by(pair_class) %>%
  summarise(
    n_pairs = n(),
    mean_p_distance = mean(p_distance),
    median_p_distance = median(p_distance),
    mean_compared_sites = mean(n_compared),
    .groups = "drop"
  )

write_csv(pair_classes, file.path(sel_dir, "nav_pairwise_protein_distances.csv"))
write_csv(pair_summary, file.path(sel_dir, "nav_pairwise_protein_distance_summary.csv"))

# ---- Derived Burden Summaries ----

burden_tbl <- bind_rows(
  read_burden(file.path(sel_dir, "ortholog_polarized_burden_by_seq.csv"), "deep_ortholog"),
  read_burden(file.path(sel_dir, "outgroup_polarized_burden_by_seq.csv"), "caddis_outgroup")
)

if (nrow(burden_tbl) > 0 && "sample_id" %in% names(burden_tbl)) {
  burden_tbl <- burden_tbl %>%
    left_join(sample_ann %>% select(sample_id, excluded_primary), by = "sample_id") %>%
    filter(!coalesce(excluded_primary, FALSE))

  burden_summary <- burden_tbl %>%
    group_by(source, focal) %>%
    summarise(
      n_seq = n(),
      mean_derived = mean(n_derived, na.rm = TRUE),
      median_derived = median(n_derived, na.rm = TRUE),
      mean_total_changes = mean(n_total_changes, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  burden_summary <- tibble()
}

write_csv(burden_summary, file.path(sel_dir, "nav_burden_support_summary.csv"))

# ---- Support for mutations from deep orthologs and caddis samples ----

subs_file <- file.path(sel_dir, "ortholog_polarized_substitutions.csv")
if (!file.exists(subs_file)) {
  subs_file <- file.path(sel_dir, "outgroup_polarized_substitutions.csv")
}

if (file.exists(subs_file)) {
  subs_tbl <- read_csv(subs_file, show_col_types = FALSE)

  if (!"domain" %in% names(subs_tbl)) {
    subs_tbl$domain <- NA_character_
  }
  if (!"in_ttx_window" %in% names(subs_tbl)) {
    subs_tbl$in_ttx_window <- NA
  }

  taxon_support <- subs_tbl %>%
    left_join(sample_ann %>% select(seq_id, excluded_primary), by = "seq_id") %>%
    filter(!coalesce(excluded_primary, FALSE), polarization == "Derived") %>%
    group_by(species_label, focal, lim_pos, domain, in_ttx_window) %>%
    summarise(n_sequences = n_distinct(seq_id), n_derived_calls = n(), .groups = "drop")
} else {
  taxon_support <- tibble()
}

write_csv(taxon_support, file.path(sel_dir, "nav_taxon_derived_support.csv"))

# ---- Selection direction  ----

if (file.exists(site_context_file)) {
  site_context <- read_csv(site_context_file, show_col_types = FALSE)

  if (all(c("fel_p_value", "fel_direction") %in% names(site_context))) {
    fel_summary <- site_context %>%
      filter(!is.na(fel_p_value), fel_p_value <= 0.05) %>%
      count(method = "FEL", direction = fel_direction, name = "n_sites")
  } else {
    fel_summary <- tibble()
  }

  if ("meme_p_value" %in% names(site_context)) {
    meme_summary <- site_context %>%
      filter(!is.na(meme_p_value), meme_p_value <= 0.05) %>%
      summarise(method = "MEME", direction = "episodic", n_sites = n())
  } else {
    meme_summary <- tibble()
  }

  selection_summary <- bind_rows(fel_summary, meme_summary)
} else {
  selection_summary <- tibble()
}

write_csv(selection_summary, file.path(sel_dir, "nav_selection_direction_summary.csv"))