library(tidyverse)
library(Biostrings)
library(msa)

rep_pep    <- "pipeline/06_variant_analysis/representatives/all_samples_rep_protein.fasta"
# Reference proteins come from the TransDecoder run on refs_84.fasta, not the nucleotide source
ref_pep    <- "pipeline/04_transdecoder/Transdecoder/p1_noflag/refs_84.fasta.p1.pep"
meta_file  <- "pipeline/06_variant_analysis/representatives/all_samples_rep_metadata.csv"
aln_dir    <- "pipeline/06_variant_analysis/alignments"

dir.create(aln_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load sequences ----

seqs_reps <- readAAStringSet(rep_pep)
seqs_refs <- readAAStringSet(ref_pep)

# Shorten names: keep first token only (already done in build_nach_dataset.R)
names(seqs_reps) <- sub("^(\\S+).*", "\\1", names(seqs_reps))
names(seqs_refs) <- sub("^(\\S+).*", "\\1", names(seqs_refs))

# Combine samples + references, deduplicate by sequence ID
# refs_84.fasta.p1.pep may contain some of the same sample sequences - drop those
ref_only <- seqs_refs[!names(seqs_refs) %in% names(seqs_reps)]
combined <- c(seqs_reps, ref_only)

# ---- ClustalOmega alignment ----

aln <- msa(combined, method = "ClustalOmega", type = "protein", verbose = FALSE)

# Convert to AAStringSet
aln_ss <- as(aln, "AAStringSet")

# Save full alignment
writeXStringSet(aln_ss, file.path(aln_dir, "nach_protein_alignment.fasta"))

# ---- Trim highly gapped ends ----

# Remove columns where >70% of sequences have a gap
aln_mat <- as.matrix(aln_ss)
gap_frac <- colMeans(aln_mat == "-")
keep_cols <- gap_frac < 0.70

aln_trimmed_mat <- aln_mat[, keep_cols, drop = FALSE]
aln_trimmed_ss <- AAStringSet(apply(aln_trimmed_mat, 1, paste, collapse = ""))
names(aln_trimmed_ss) <- rownames(aln_trimmed_mat)

writeXStringSet(aln_trimmed_ss, file.path(aln_dir, "nach_protein_alignment_trimmed.fasta"))