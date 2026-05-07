library(tidyverse)
library(Biostrings)
library(msa)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

rep_pep_file  <- file.path(root, "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_protein.fasta")
rep_cds_file  <- file.path(root, "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_cds.fasta")
rep_meta_file <- file.path(root, "nach_extraction/05_representative_sequences/outputs/representatives/all_samples_rep_metadata.csv")
ref_pep_file  <- file.path(root, "nach_extraction/04_transdecoder/outputs/transdecoder/Transdecoder/p1_noflag/refs_84.fasta.p1.pep")

aln_dir <- file.path(root, "nav_selection/outputs/alignments")
dir.create(aln_dir, showWarnings = FALSE, recursive = TRUE)

protein_aln_file     <- file.path(aln_dir, "nav_protein_alignment.fasta")
protein_trimmed_file <- file.path(aln_dir, "nav_protein_alignment_trimmed.fasta")
codon_aln_file       <- file.path(aln_dir, "nav_codon_alignment.fasta")
codon_filtered_file  <- file.path(aln_dir, "nav_codon_alignment_filtered.fasta")

gap_limit     <- 0.70
occ_threshold <- 0.60

# ---- Protein alignment ----

rep_pep <- readAAStringSet(rep_pep_file)
ref_pep <- readAAStringSet(ref_pep_file)
names(rep_pep) <- sub("^(\\S+).*", "\\1", names(rep_pep))
names(ref_pep) <- sub("^(\\S+).*", "\\1", names(ref_pep))

combined_pep <- c(rep_pep, ref_pep[!names(ref_pep) %in% names(rep_pep)])
protein_aln  <- as(msa(combined_pep, method = "ClustalOmega", type = "protein", verbose = FALSE), "AAStringSet")
writeXStringSet(protein_aln, protein_aln_file)

# Drop columns with too many gaps before downstream use
protein_mat <- as.matrix(protein_aln)
keep_cols   <- colMeans(protein_mat == "-") < gap_limit
trimmed_mat <- protein_mat[, keep_cols, drop = FALSE]
trimmed_aln <- AAStringSet(apply(trimmed_mat, 1, paste, collapse = ""))
names(trimmed_aln) <- rownames(trimmed_mat)
writeXStringSet(trimmed_aln, protein_trimmed_file)

# ---- Codon alignment ----

rep_meta <- read_csv(rep_meta_file, show_col_types = FALSE)
rep_cds  <- readDNAStringSet(rep_cds_file)
names(rep_cds) <- sub("^(\\S+).*", "\\1", names(rep_cds))

sample_seq_ids <- rep_meta$seq_id[rep_meta$seq_id %in% names(protein_aln)]
seq_to_sample  <- setNames(rep_meta$sample_id, rep_meta$seq_id)

# Check each protein-aligned sample and place its CDS codons next to the matching residues
codon_strings <- character(length(sample_seq_ids))
for (i in seq_along(sample_seq_ids)) {
  seq_id <- sample_seq_ids[i]
  aa_full <- strsplit(as.character(protein_aln[[seq_id]]), "")[[1]]
  dna_seq <- toupper(as.character(rep_cds[[seq_id]]))
  usable_nt <- nchar(dna_seq) - (nchar(dna_seq) %% 3)
  codons <- substring(dna_seq, seq(1, usable_nt, by = 3), seq(3, usable_nt, by = 3))

  codon_idx <- 1
  codon_full <- character(length(aa_full))
  for (j in seq_along(aa_full)) {
    if (aa_full[j] == "-") {
      codon_full[j] <- "---"
    } else if (codon_idx > length(codons)) {
      codon_full[j] <- "NNN"
    } else {
      codon_full[j] <- codons[codon_idx]
      codon_idx <- codon_idx + 1
    }
  }
  codon_strings[i] <- paste(codon_full[keep_cols], collapse = "")
}

codon_aln <- DNAStringSet(codon_strings)
names(codon_aln) <- seq_to_sample[sample_seq_ids]
writeXStringSet(codon_aln, codon_aln_file)

# ---- Codon occupancy filter ----

codon_mat <- as.matrix(codon_aln)
codon_n   <- ncol(codon_mat) / 3

occupancy <- vapply(seq_len(nrow(codon_mat)), function(r) {
  nt <- codon_mat[r, ]
  codon_vec <- vapply(seq_len(codon_n), function(k) {
    paste0(nt[((k - 1) * 3 + 1):((k - 1) * 3 + 3)], collapse = "")
  }, character(1))
  sum(codon_vec != "---") / codon_n
}, numeric(1))

writeXStringSet(codon_aln[occupancy >= occ_threshold], codon_filtered_file)
