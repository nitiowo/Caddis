library(tidyverse)
library(Biostrings)
library(DECIPHER)
library(ape)
library(ggtree)
library(phangorn)

ingroup_fa  <- "data/COI/coi_sequences.fasta"
outgroup_fa <- "data/COI/outgroups/limnephilid_outgroups.fasta"
seq_meta    <- "09_coi_phylogeography/outputs/tables/coi_sequence_metadata_table.csv"

aln_dir  <- "09_coi_phylogeography/outputs/alignments"
dist_dir <- "09_coi_phylogeography/outputs/distances"
fig_dir  <- "09_coi_phylogeography/outputs/figures"

for (d in c(aln_dir, dist_dir, fig_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

aln_full <- file.path(aln_dir,  "coi_with_outgroups_alignment.fasta")
tr_nwk   <- file.path(dist_dir, "coi_rooted_ml_tree.nwk")
fig_tr   <- file.path(fig_dir,  "coi_rooted_ml_tree.png")

# ---- Load and combine ----

ing      <- readDNAStringSet(ingroup_fa)
out      <- readDNAStringSet(outgroup_fa)
all_seqs <- c(ing, out)

# ---- Realign with DECIPHER on the combined set ----

aln      <- AlignSeqs(all_seqs, processors = NULL, verbose = FALSE)
aln_mat  <- do.call(rbind, strsplit(as.character(aln), ""))
col_cov  <- colMeans(aln_mat != "-")
core_mat <- aln_mat[, col_cov >= 0.8, drop = FALSE]
core_seq <- apply(core_mat, 1, paste0, collapse = "")
writeLines(unlist(map2(names(aln), core_seq, ~ c(paste0(">", .x), .y))), aln_full)

dna <- read.dna(aln_full, format = "fasta")
phy <- as.phyDat(dna)

# ---- ML tree (phangorn) ----

set.seed(42)
dm    <- dist.dna(dna, model = "TN93", pairwise.deletion = TRUE)
nj_tr <- nj(dm)
nj_tr$edge.length[nj_tr$edge.length < 0] <- 0

fit0 <- pml(nj_tr, data = phy, k = 4)
fit  <- optim.pml(fit0, model = "GTR", optGamma = TRUE, optEdge = TRUE, optInv = FALSE,
                  control = pml.control(trace = 0))

ml_tr <- fit$tree

# Root on outgroup MRCA
og_tips   <- names(out)
ml_rooted <- root(ml_tr, outgroup = og_tips, resolve.root = TRUE)

write.tree(ml_rooted, tr_nwk)

# ---- Rooted tree colored by taxon ----

seq_tbl <- read_csv(seq_meta, show_col_types = FALSE)
tip_tbl <- tibble(label = ml_rooted$tip.label) %>%
  left_join(seq_tbl %>% select(label = seq_name, species_label, locality_short),
            by = "label") %>%
  mutate(
    species_label = ifelse(is.na(species_label), "Outgroup (Limnephilidae)", species_label),
    is_og         = species_label == "Outgroup (Limnephilidae)",
    og_label      = ifelse(is_og, gsub("_", " ", label), "")
  )

p <- ggtree(ml_rooted, ladderize = TRUE) %<+% tip_tbl +
  geom_tippoint(aes(color = species_label, shape = is_og), size = 1.6, na.rm = TRUE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  geom_tiplab(aes(label = og_label), size = 2, hjust = -0.05) +
  labs(color = "Taxon", title = "Rooted COI ML tree with limnephilid outgroups") +
  theme_tree2() +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(ncol = 1))

ggsave(fig_tr, p, width = 10, height = 12, dpi = 300)
