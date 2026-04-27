library(tidyverse)
library(Biostrings)

td_dir    <- "pipeline/04_transdecoder/Transdecoder/p1_noflag"
ref_pep   <- "reference/sodium_channel/refs_84.fasta"
out_dir   <- "pipeline/06_variant_analysis/representatives"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Parse TransDecoder headers ----

# Extract fields from a TransDecoder .p1.pep header line
parse_td_header <- function(header) {
  seq_id   <- sub("^>?(\\S+).*", "\\1", header)
  orf_type <- sub(".*ORF type:(\\S+) .*", "\\1", header)
  score    <- as.numeric(sub(".*score=([0-9.]+).*", "\\1", header))
  len      <- as.integer(sub(".*len:([0-9]+).*", "\\1", header))
  data.frame(seq_id, orf_type, score, len, stringsAsFactors = FALSE)
}

# ---- Read all .p1.pep files ----

pep_files <- list.files(td_dir, pattern = "\\.p1\\.pep$", full.names = TRUE)

# Exclude reference files — user the refs_84_pep.sh script for those
sample_pep_files <- pep_files[!grepl("refs_84|dmel-musdo-lilu", basename(pep_files))]

read_pep_file <- function(f) {
  seqs    <- readAAStringSet(f)
  headers <- names(seqs)
  parsed  <- map_dfr(headers, parse_td_header)
  parsed$source_file <- basename(f)
  # Extract sample_id: CF_3 from "LiNo_CF_3_DN..." -> "CF_3"
  parsed$sample_raw <- sub(".*_(CF_[0-9]+)_.*", "\\1", parsed$seq_id)
  # Unify to 3-digit format: CF_3 -> CF003
  parsed$sample_id <- sprintf("CF%03d", as.integer(sub("CF_", "", parsed$sample_raw)))
  parsed
}

all_meta <- map_dfr(sample_pep_files, read_pep_file)

# ---- Select one representative per sample ----

# Priority: complete > partial > internal
# Within complete: pick highest scoring, ties broken by length
orf_rank <- c("complete" = 1, "5prime_partial" = 2, "3prime_partial" = 3, "internal" = 4)

reps <- all_meta %>%
  mutate(orf_rank = orf_rank[orf_type]) %>%
  group_by(sample_id) %>%
  arrange(orf_rank, desc(score), desc(len), .by_group = TRUE) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(is_representative = TRUE)

# Mark all sequences as representative or not
all_meta <- all_meta %>%
  left_join(reps %>% dplyr::select(seq_id, is_representative), by = "seq_id") %>%
  mutate(is_representative = replace_na(is_representative, FALSE))

# ---- Load sequences for output ----

all_pep_seqs <- do.call(c, lapply(sample_pep_files, readAAStringSet))

# Trim seq_id in seqs to match parsed seq_id
names(all_pep_seqs) <- sub("^(\\S+).*", "\\1", names(all_pep_seqs))

rep_ids    <- reps$seq_id
isoform_ids <- all_meta$seq_id

rep_pep    <- all_pep_seqs[rep_ids]
isoform_pep <- all_pep_seqs[isoform_ids]

# ---- Load CDS sequences ----

cds_files <- sub("\\.pep$", ".cds", sample_pep_files)
cds_files <- cds_files[file.exists(cds_files)]

all_cds_seqs <- do.call(c, lapply(cds_files, readDNAStringSet))
names(all_cds_seqs) <- sub("^(\\S+).*", "\\1", names(all_cds_seqs))

rep_cds     <- all_cds_seqs[intersect(rep_ids, names(all_cds_seqs))]
isoform_cds <- all_cds_seqs[intersect(isoform_ids, names(all_cds_seqs))]

# ---- Write FASTA outputs ----

writeXStringSet(rep_pep,     file.path(out_dir, "all_samples_rep_protein.fasta"))
writeXStringSet(rep_cds,     file.path(out_dir, "all_samples_rep_cds.fasta"))
writeXStringSet(isoform_pep, file.path(out_dir, "all_isoforms_protein.fasta"))
writeXStringSet(isoform_cds, file.path(out_dir, "all_isoforms_cds.fasta"))

# ---- Write metadata table ----

metadata_out <- reps %>%
  dplyr::select(sample_id, sample_raw, seq_id, source_file, orf_type, score, len, is_representative) %>%
  dplyr::rename(
    orig_file    = source_file,
    orf_score    = score,
    seq_len_aa   = len
  )

write_csv(metadata_out, file.path(out_dir, "all_samples_rep_metadata.csv"))