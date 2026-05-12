library(tidyverse)
library(Biostrings)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

# ---- Load Stuff ----

aln_file <- file.path(root, "nav_selection/outputs/alignments/nav_protein_alignment_trimmed.fasta")
rep_meta_file <- file.path(
  root,
  "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_metadata.csv"
)
rnaseq_file <- file.path(root, "data/metadata/MASTER_caddis_rnaseq_samples.csv")
ann_file <- file.path(root, "nav_selection/outputs/alignments/nav_domain_annotation.csv")

out_dir <- file.path(root, "nav_selection/outputs/selection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

subs_file <- file.path(out_dir, "outgroup_polarized_substitutions.csv")
burden_seq_file <- file.path(out_dir, "outgroup_polarized_burden_by_seq.csv")
burden_group_file <- file.path(out_dir, "outgroup_polarized_burden_by_group.csv")
contrast_file <- file.path(out_dir, "pore_window_focal_nonfocal_contrast.csv")
fisher_file <- file.path(out_dir, "pore_window_fisher_test.csv")

# ---- Functions ----

norm_sample_id <- function(value) {
  sample_num <- str_match(str_trim(as.character(value)), "(?i)CF[-_ ]*0*([0-9]+)")[, 2]
  ifelse(is.na(sample_num), NA_character_, sprintf("CF%03d", as.integer(sample_num)))
}

modal_state <- function(values) {
  clean_vals <- values[!values %in% c("-", "X")]
  if (length(clean_vals) == 0) {
    return(NA_character_)
  }
  names(sort(table(clean_vals), decreasing = TRUE))[1]
}

read_sample_ann <- function(rep_meta, rnaseq_path) {
  rnaseq <- read_csv(rnaseq_path, show_col_types = FALSE) %>%
    mutate(
      sample_id = norm_sample_id(SAMPLE),
      species_label = paste(str_trim(`Genbank genus`), str_trim(`Genbank species`)),
      focal = str_trim(`Genbank species`) %in% c("flavastellus", "externus")
    ) %>%
    dplyr::select(sample_id, species_label, focal, Location)

  rep_meta %>%
    dplyr::select(sample_id, seq_id, orf_type) %>%
    left_join(rnaseq, by = "sample_id") %>%
    dplyr::rename(locality = Location)
}

# ---- Load and subset ----

aln <- readAAStringSet(aln_file)
names(aln) <- sub("^(\\S+).*", "\\1", names(aln))
aln_mat <- as.matrix(aln)

rep_meta <- read_csv(rep_meta_file, show_col_types = FALSE)
sample_ann <- read_sample_ann(rep_meta, rnaseq_file)

sample_ids <- sample_ann$seq_id[sample_ann$seq_id %in% rownames(aln_mat)]
non_sample_ids <- setdiff(rownames(aln_mat), sample_ids)

lim_ids <- non_sample_ids[str_detect(non_sample_ids, "^ENSLLST")] # Lim luna
out_ids <- setdiff(non_sample_ids, lim_ids)                        # all other refs

if (length(lim_ids) == 0) {
  stop("No Limnephilus lunatus ref IDs")
}
if (length(out_ids) == 0) {
  stop("No non-Limnephilus ref IDs")
}

ann_tbl <- read_csv(ann_file, show_col_types = FALSE) %>%
  dplyr::select(aln_col, lim_pos, domain, in_ttx_window)

# ---- Polarized substitutions ----

# Use Lim luna as coordinate reference, skip gaps
lim_anchor <- aln_mat[lim_ids[1], ]
lim_pos <- cumsum(lim_anchor != "-")
lim_pos[lim_anchor == "-"] <- NA_integer_

# Consensus amino acid at each column for outgroups
lim_mode <- apply(aln_mat[lim_ids, , drop = FALSE], 2, modal_state)
out_mode <- apply(aln_mat[out_ids, , drop = FALSE], 2, modal_state)

# Reference coordinate table
coord_tbl <- tibble(
  aln_col = seq_len(ncol(aln_mat)),
  lim_pos = lim_pos,
  lim_ref = lim_mode,
  out_ref = out_mode
) %>%
  filter(!is.na(lim_pos), !is.na(lim_ref))

long_tbl <- map_dfr(sample_ids, function(seq_id) {
  aa_vals <- aln_mat[seq_id, ]
  tibble(
    seq_id = seq_id,
    aln_col = seq_along(aa_vals),
    aa = aa_vals
  )
}) %>%
  left_join(sample_ann, by = "seq_id") %>%
  left_join(coord_tbl, by = "aln_col") %>%
  filter(!is.na(lim_pos), !aa %in% c("-", "X"))

# Derived - limnephilids match the outgroup but this sample doesn't
# Ancestral-retained - this sample matches the outgroup instead of limnephilids
# Third-state - everyone is different from each other
# Outgroup unresolved - not enough outgroup coverage to polarize
subs_tbl <- long_tbl %>%
  filter(aa != lim_ref) %>%
  mutate(
    polarization = case_when(
      is.na(out_ref) ~ "Outgroup unresolved",
      aa == out_ref ~ "Ancestral-retained",
      lim_ref == out_ref & aa != out_ref ~ "Derived",
      lim_ref != out_ref & aa != out_ref ~ "Third-state",
      TRUE ~ "Unclassified"
    )
  ) %>%
  dplyr::select(-lim_pos) %>%
  left_join(ann_tbl, by = "aln_col") %>%
  arrange(desc(focal), species_label, lim_pos, seq_id)

write_csv(subs_tbl, subs_file)

# ---- Mutation burden summaries ----

burden_seq <- subs_tbl %>%
  group_by(seq_id, sample_id, species_label, focal) %>%
  summarise(
    n_total_changes = n(),
    n_derived = sum(polarization == "Derived"),
    n_third_state = sum(polarization == "Third-state"),
    n_ancestral_retained = sum(polarization == "Ancestral-retained"),
    n_unresolved = sum(polarization == "Outgroup unresolved"),
    .groups = "drop"
  )

burden_group <- burden_seq %>%
  group_by(focal) %>%
  summarise(
    n_seq = n(),
    mean_derived = mean(n_derived),
    mean_third_state = mean(n_third_state),
    mean_ancestral_retained = mean(n_ancestral_retained),
    mean_unresolved = mean(n_unresolved),
    .groups = "drop"
  )

write_csv(burden_seq, burden_seq_file)
write_csv(burden_group, burden_group_file)

# ---- Pore-window contrast ----

# Sum changes inside vs outside TTX pore window
contrast_tbl <- subs_tbl %>%
  mutate(
    focal_label = if_else(focal, "focal", "non_focal"),
    window_label = if_else(in_ttx_window, "pore_window", "outside_window")
  ) %>%
  group_by(focal_label, window_label, polarization) %>%
  summarise(n_changes = n(), .groups = "drop") %>%
  arrange(focal_label, window_label, polarization)

write_csv(contrast_tbl, contrast_file)

test_tbl <- subs_tbl %>%
  filter(polarization %in% c("Derived", "Ancestral-retained")) %>%
  mutate(
    focal_label = if_else(focal, "focal", "non_focal"),
    window_label = if_else(in_ttx_window, "pore_window", "outside_window")
  )

fisher_rows <- list()

for (group_label in c("focal", "non_focal")) {
  group_tbl <- test_tbl %>% filter(focal_label == group_label)

  count_mat <- matrix(0, nrow = 2, ncol = 2)
  rownames(count_mat) <- c("Derived", "Ancestral-retained")
  colnames(count_mat) <- c("pore_window", "outside_window")

  raw_tab <- table(group_tbl$polarization, group_tbl$window_label)
  count_mat[rownames(raw_tab), colnames(raw_tab)] <- raw_tab
  fisher_res <- fisher.test(count_mat)

  fisher_rows[[group_label]] <- tibble(
    focal_label = group_label,
    odds_ratio = as.numeric(fisher_res$estimate),
    p_value = fisher_res$p.value,
    derived_pore = count_mat["Derived", "pore_window"],
    derived_outside = count_mat["Derived", "outside_window"],
    ancestral_pore = count_mat["Ancestral-retained", "pore_window"],
    ancestral_outside = count_mat["Ancestral-retained", "outside_window"]
  )
}

bind_rows(fisher_rows) %>% write_csv(fisher_file)