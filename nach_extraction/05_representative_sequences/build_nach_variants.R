library(tidyverse)
library(Biostrings)

aln_file  <- "pipeline/06_variant_analysis/alignments/nach_protein_alignment_trimmed.fasta"
meta_file <- "pipeline/06_variant_analysis/representatives/all_samples_rep_metadata.csv"
rnaseq_f  <- "data/metadata/MASTER_caddis_rnaseq_samples.csv"
var_dir   <- "pipeline/06_variant_analysis/variants"
hap_dir   <- "pipeline/06_variant_analysis/haplotypes"

dir.create(var_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(hap_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load alignment and metadata ----

aln      <- readAAStringSet(aln_file)
names(aln) <- sub("^(\\S+).*", "\\1", names(aln))

rep_meta <- read_csv(meta_file, show_col_types = FALSE)
rnaseq   <- read_csv(rnaseq_f,  show_col_types = FALSE) %>%
  mutate(
    sample_id    = sprintf("CF%03d", as.integer(sub("CF", "", str_trim(SAMPLE)))),
    species_label = paste(str_trim(`Genbank genus`), str_trim(`Genbank species`)),
    focal         = str_trim(`Genbank species`) %in% c("flavastellus", "externus")
  )

# Build sequence annotation table for all entries in the alignment
# Distinguish sample sequences vs reference sequences
is_sample <- names(aln) %in% rep_meta$seq_id

sample_ann <- rep_meta %>%
  left_join(rnaseq %>% dplyr::select(sample_id, species_label, focal, Location),
            by = "sample_id") %>%
  dplyr::select(seq_id, sample_id, species_label, focal, Location, orf_type) %>%
  dplyr::rename(locality = Location)

# Annotate refs
ref_ids  <- names(aln)[!is_sample]
ref_ann  <- data.frame(
  seq_id       = ref_ids,
  sample_id    = NA_character_,
  species_label = "reference",
  focal         = FALSE,
  locality      = NA_character_,
  orf_type      = "reference",
  stringsAsFactors = FALSE
)

ann <- bind_rows(sample_ann, ref_ann)

# ---- Variable site extraction ----

aln_mat <- as.matrix(aln)
n_seq   <- nrow(aln_mat)

# Only consider columns with at least 2 different non-gap characters
is_variable <- apply(aln_mat, 2, function(col) {
  residues <- col[col != "-" & col != "X"]
  length(unique(residues)) > 1
})

var_cols <- which(is_variable)

# ---- Protein variant matri ----

sample_rows <- rownames(aln_mat) %in% rep_meta$seq_id

var_mat <- aln_mat[sample_rows, var_cols, drop = FALSE]
var_df  <- as.data.frame(var_mat)
names(var_df) <- paste0("pos_", var_cols)
var_df$seq_id <- rownames(var_mat)

variant_matrix <- var_df %>%
  left_join(sample_ann, by = "seq_id") %>%
  dplyr::select(seq_id, sample_id, species_label, focal, locality, orf_type,
                everything())

write_csv(variant_matrix, file.path(var_dir, "protein_variant_matrix.csv"))

# # Lifl vs non-Lifl matrix
# lifl_mat <- variant_mat %>%