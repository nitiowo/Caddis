library(tidyverse)
library(Biostrings)
library(jsonlite)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

aln_file <- file.path(root, "nav_selection/outputs/alignments/nav_protein_alignment.fasta")
rep_meta_file <- file.path(root, "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_metadata.csv")
site_meta_file <- file.path(root, "nav_selection/config/nav_site_metadata.csv")
hyphy_dir <- file.path(root, "nav_selection/outputs/model_raw/hyphy")

out_dir <- file.path(root, "nav_selection/outputs/selection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

site_context_file <- file.path(out_dir, "nav_site_context.csv")
hyphy_site_file <- file.path(out_dir, "nav_hyphy_site_summary.csv")

# ---- Functions ----

safe_col <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  gsub("^_+|_+$", "", value)
}

species_prefix <- function(seq_id) {
  sub("_CF_.*", "", seq_id)
}

read_hyphy_sites <- function(path, method) {
  if (!file.exists(path)) {
    return(tibble())
  }

  fit <- fromJSON(path, simplifyVector = FALSE)
  mat <- do.call(rbind, fit$MLE$content[["0"]])

  if (method == "fel") {
    alpha <- as.numeric(mat[, 1])
    beta <- as.numeric(mat[, 2])
    p_value <- as.numeric(mat[, 4])

    return(tibble(
      hyphy_site = seq_len(nrow(mat)),
      fel_alpha = alpha,
      fel_beta = beta,
      fel_lrt = as.numeric(mat[, 3]),
      fel_p_value = p_value,
      fel_direction = if_else(beta > alpha, "positive", "negative"),
      fel_significant_p_0_05 = p_value <= 0.05
    ))
  }

  p_value <- as.numeric(mat[, 6])
  tibble(
    hyphy_site = seq_len(nrow(mat)),
    meme_p_value = p_value,
    meme_significant_p_0_05 = p_value <= 0.05
  )
}

# ---- Alignment Sites ----

aln <- readAAStringSet(aln_file)
names(aln) <- sub("^(\\S+).*", "\\1", names(aln))
aln_mat <- as.matrix(aln)

gap_frac <- colMeans(aln_mat == "-")
codon_keep_cols <- which(gap_frac < 0.70)

site_tbl <- tibble(
  aln_col = seq_len(ncol(aln_mat)),
  codon_alignment_site = match(aln_col, codon_keep_cols),
  gap_frac = gap_frac
)

if (file.exists(site_meta_file)) {
  site_meta <- read_csv(site_meta_file, show_col_types = FALSE) %>%
    mutate(aln_col = as.integer(aln_col))

  extra_cols <- setdiff(names(site_meta), names(site_tbl))
  site_tbl <- site_tbl %>%
    left_join(site_meta %>% select(aln_col, all_of(extra_cols)), by = "aln_col")
}

# ---- Representative Amino Acids ----

rep_meta <- read_csv(rep_meta_file, show_col_types = FALSE)

if ("is_representative" %in% names(rep_meta)) {
  rep_meta <- rep_meta %>%
    filter(toupper(as.character(is_representative)) == "TRUE")
}

rep_meta <- rep_meta %>%
  filter(seq_id %in% rownames(aln_mat)) %>%
  arrange(sample_id) %>%
  mutate(
    species_prefix = species_prefix(seq_id),
    aa_col = safe_col(paste0("rep_aa_", sample_id, "_", species_prefix))
  )

for (rep_row in seq_len(nrow(rep_meta))) {
  site_tbl[[rep_meta$aa_col[rep_row]]] <- aln_mat[rep_meta$seq_id[rep_row], site_tbl$aln_col]
}

# ---- HyPhy Site Metrics ----

fel_tbl <- read_hyphy_sites(file.path(hyphy_dir, "fel.json"), "fel")
meme_tbl <- read_hyphy_sites(file.path(hyphy_dir, "meme.json"), "meme")

if (nrow(fel_tbl) > 0 && nrow(meme_tbl) > 0) {
  hyphy_tbl <- full_join(fel_tbl, meme_tbl, by = "hyphy_site")
} else if (nrow(fel_tbl) > 0) {
  hyphy_tbl <- fel_tbl
} else if (nrow(meme_tbl) > 0) {
  hyphy_tbl <- meme_tbl
} else {
  hyphy_tbl <- tibble()
}

if (nrow(hyphy_tbl) > 0) {
  write_csv(hyphy_tbl, hyphy_site_file)
  site_tbl <- site_tbl %>%
    left_join(hyphy_tbl, by = c("codon_alignment_site" = "hyphy_site"))
}

write_csv(site_tbl, site_context_file)