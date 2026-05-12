library(tidyverse)
library(Biostrings)
library(DECIPHER)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

# ---- Load Stuff ----

aln_file <- file.path(root, "nav_selection/outputs/alignments/nav_protein_alignment_trimmed.fasta")
orth_file <- file.path(root, "Reference/sodium_channel/allrefspecies_YY_alignment.fasta")
rep_meta_file <- file.path(
  root,
  "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_metadata.csv"
)
rnaseq_file <- file.path(root, "data/metadata/MASTER_caddis_rnaseq_samples.csv")
ann_file <- file.path(root, "nav_selection/outputs/alignments/nav_domain_annotation.csv")

aln_dir <- file.path(root, "nav_selection/outputs/alignments")
out_dir <- file.path(root, "nav_selection/outputs/selection")
dir.create(aln_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_aln_file <- file.path(aln_dir, "nav_protein_with_orthologs_alignment.fasta")
subs_file <- file.path(out_dir, "ortholog_polarized_substitutions.csv")
burden_seq_file <- file.path(out_dir, "ortholog_polarized_burden_by_seq.csv")
burden_group_file <- file.path(out_dir, "ortholog_polarized_burden_by_group.csv")
group_test_file <- file.path(out_dir, "ortholog_polarized_group_test.csv")

# ---- Select groups ----

deep_outgroups <- c(
  "Bmori",
  "Dmel|paraPA",
  "Housefly|X96668.1:146340",
  "Tcast|XM_015982008.1"
)
lim_refs <- c(
  "gAfus_ENSSFIG00005005828.1",
  "gGpel_ENSGPLG00000012611.1",
  "gLmar_ENSLMMG00005001277.1",
  "gLlun_ENSLLSG00015002091.1",
  "gLrho_ENSLRHG00005003296.1"
)
exclude_orth <- c(
  "Annulipalpia_sp",
  "LF_TRINITY_DN18073_c8_g1_i4_12",
  "HO_TRINITY_DN18886_c1_g1_i4_3_DN18886_c1_g2_i10_4"
)
anchor_id <- "gLlun_ENSLLSG00015002091.1"

# ---- Functions ----

norm_sample_id <- function(value) {
  sample_num <- str_match(str_trim(as.character(value)), "(?i)CF[-_ ]*0*([0-9]+)")[, 2]
  ifelse(is.na(sample_num), NA_character_, sprintf("CF%03d", as.integer(sample_num)))
}

# Get most common AA across groups
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

# ---- Profile alignment ----

ing_aln <- readAAStringSet(aln_file)
names(ing_aln) <- sub("^(\\S+).*", "\\1", names(ing_aln))

orth_aln <- readAAStringSet(orth_file)
orth_aln <- orth_aln[!names(orth_aln) %in% exclude_orth]

# Align and preserve gaps
combined_aln <- AlignProfiles(ing_aln, orth_aln, processors = NULL)
writeXStringSet(combined_aln, out_aln_file)
aln_mat <- as.matrix(combined_aln)

if (!anchor_id %in% rownames(aln_mat)) {
  stop("No Lim luna")
}
anchor <- aln_mat[anchor_id, ]
lim_pos <- cumsum(anchor != "-")
lim_pos[anchor == "-"] <- NA_integer_

og_present <- intersect(deep_outgroups, rownames(aln_mat))
if (length(og_present) == 0) {
  stop("No ancestral outgroups")
}

lim_present <- intersect(lim_refs, rownames(aln_mat))
if (length(lim_present) == 0) {
  stop("No limnephilidae groups")
}

og_state <- apply(aln_mat[og_present, , drop = FALSE], 2, modal_state)
lim_state <- apply(aln_mat[lim_present, , drop = FALSE], 2, modal_state)
n_og_resolved <- apply(
  aln_mat[og_present, , drop = FALSE],
  2,
  function(values) sum(!values %in% c("-", "X"))
)

# Reference coordinate table
coord_tbl <- tibble(
  aln_col = seq_len(ncol(aln_mat)),
  lim_pos = lim_pos,
  lim_ref = lim_state,
  out_ref = og_state,
  n_og_resolved = n_og_resolved
) %>%
  filter(!is.na(lim_pos), !is.na(lim_ref))

# ---- Polarized substitutions ----

rep_meta <- read_csv(rep_meta_file, show_col_types = FALSE)
sample_ann <- read_sample_ann(rep_meta, rnaseq_file)
sample_ids <- sample_ann$seq_id[sample_ann$seq_id %in% rownames(aln_mat)]

ann_tbl <- read_csv(ann_file, show_col_types = FALSE) %>%
  dplyr::select(lim_pos, domain, in_ttx_window)

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
  left_join(ann_tbl, by = "lim_pos") %>%
  arrange(desc(focal), species_label, lim_pos, seq_id)

write_csv(subs_tbl, subs_file)

# Count substitutions
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

# Count LF
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

# Wilcoxon test
group_test <- wilcox.test(n_derived ~ focal, data = burden_seq)
group_test_tbl <- tibble(
  test = "wilcoxon_derived_burden_by_focal_status",
  statistic = as.numeric(group_test$statistic),
  p_value = group_test$p.value
)

write_csv(burden_seq, burden_seq_file)
write_csv(burden_group, burden_group_file)
write_csv(group_test_tbl, group_test_file)
