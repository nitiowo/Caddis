library(tidyverse)
library(Biostrings)
library(DECIPHER)
library(ape)
library(ggtree)
library(phangorn)

# Inputs (CLEAN)
# meta_file: CSV with columns sample_id, species_label, locality_short, lat, lon
#            (sample_id must match FASTA sequence names exactly)
# ttx_file:  CSV with columns site, lat, lon, ttx_skin_mg
seq_file <- "data/COI/coi_sequences.fasta"
meta_file <- "data/metadata/coi_metadata.csv"
ttx_file  <- "data/ttx/newt_ttx_coords.csv"

aln_dir <- "09_coi_phylogeography/outputs/alignments"
dist_dir <- "09_coi_phylogeography/outputs/distances"
hap_dir  <- "09_coi_phylogeography/outputs/haplotypes"
fig_dir  <- "09_coi_phylogeography/outputs/figures"
tab_dir  <- "09_coi_phylogeography/outputs/tables"

for (d in c(aln_dir, dist_dir, hap_dir, fig_dir, tab_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

core_fa  <- file.path(aln_dir, "coi_core_alignment.fasta")
dm_rds   <- file.path(dist_dir, "coi_pairwise_distance_matrix.rds")
tr_nwk   <- file.path(dist_dir, "coi_nj_tree.nwk")

tab_seq  <- file.path(tab_dir, "coi_sequence_metadata_table.csv")
tab_sp   <- file.path(tab_dir, "coi_species_summary.csv")
tab_foc  <- file.path(tab_dir, "coi_flavastellus_locality_summary.csv")
tab_pair <- file.path(tab_dir, "coi_flavastellus_pairwise_distance_summary.csv")

fig_tree <- file.path(fig_dir, "coi_nj_tree.png")
fig_hap  <- file.path(fig_dir, "coi_flavastellus_haplotype_heatmap.png")

# ---- Load sequences and metadata ----

seqs <- readDNAStringSet(seq_file)
meta <- read_csv(meta_file, show_col_types = FALSE)
ttx  <- read_csv(ttx_file,  show_col_types = FALSE)

seq_tbl <- tibble(
  seq_name = names(seqs),
  seq_len  = width(seqs)
) %>%
  left_join(meta, by = c("seq_name" = "sample_id"))

seqs <- seqs[seq_tbl$seq_name]



# ---- Functions ----

hap_div <- function(h) {
  n <- length(h)
  if (n < 2) return(0)
  p <- table(h) / n
  (n / (n - 1)) * (1 - sum(p^2))
}

mean_pair_dist <- function(dm, idx) {
  if (length(idx) < 2) return(NA_real_)
  x <- as.matrix(dm)[idx, idx, drop = FALSE]
  mean(x[upper.tri(x)], na.rm = TRUE)
}

nearest_ttx <- function(lat, lon, tox) {
  d <- sqrt((tox$lat - lat)^2 + (tox$lon - lon)^2)
  i <- which.min(d)
  tibble(
    nearest_site     = tox$site[i],
    ttx_skin_mg      = tox$ttx_skin_mg[i],
    tox_distance_deg = d[i]
  )
}

# ---- Alignments ----

aln     <- AlignSeqs(seqs, processors = NULL, verbose = FALSE)
aln_mat <- do.call(rbind, strsplit(as.character(aln), ""))
col_cov <- colMeans(aln_mat != "-")
core_mat <- aln_mat[, col_cov >= 0.95, drop = FALSE]
core_seq <- apply(core_mat, 1, paste0, collapse = "")

writeLines(
  unlist(map2(names(aln), core_seq, ~ c(paste0(">", .x), .y))),
  core_fa
)

dna <- read.dna(core_fa, format = "fasta")
dm  <- dist.dna(dna, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
tr  <- midpoint(nj(as.dist(dm)))
tr$edge.length[tr$edge.length < 0] <- 0

saveRDS(dm, dm_rds)
write.tree(tr, tr_nwk)

seq_tbl <- seq_tbl %>%
  mutate(
    core_seq = core_seq,
    hap_id   = paste0("H", match(core_seq, unique(core_seq)))
  )

# ---- Species summaries ----

sp_sum <- seq_tbl %>%
  group_by(species_label) %>%
  summarise(
    n_seq               = n(),
    n_localities        = n_distinct(locality_short),
    n_haplotypes        = n_distinct(hap_id),
    haplotype_diversity = hap_div(hap_id),
    .groups = "drop"
  ) %>%
  mutate(
    nucleotide_diversity = map_dbl(species_label, function(sp) {
      mean_pair_dist(dm, which(seq_tbl$species_label == sp))
    })
  ) %>%
  arrange(desc(n_seq), species_label)

# ---- LF locality summaries  ----

foc <- seq_tbl %>% filter(species_label == "Limnephilus flavastellus")

foc_loc <- foc %>%
  group_by(locality_short) %>%
  summarise(
    lat                 = mean(lat, na.rm = TRUE),
    lon                 = mean(lon, na.rm = TRUE),
    n_seq               = n(),
    n_haplotypes        = n_distinct(hap_id),
    haplotype_diversity = hap_div(hap_id),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(nearest = list(nearest_ttx(lat, lon, ttx))) %>%
  unnest_wider(nearest) %>%
  ungroup()

lf_dm    <- dm[foc$seq_name, foc$seq_name]
pair_idx <- which(upper.tri(lf_dm), arr.ind = TRUE)

foc_within_between <- tibble(
  s1   = rownames(lf_dm)[pair_idx[, 1]],
  s2   = colnames(lf_dm)[pair_idx[, 2]],
  dist = lf_dm[pair_idx]
) %>%
  left_join(foc %>% select(seq_name, locality_short), by = c("s1" = "seq_name")) %>%
  rename(loc1 = locality_short) %>%
  left_join(foc %>% select(seq_name, locality_short), by = c("s2" = "seq_name")) %>%
  rename(loc2 = locality_short) %>%
  mutate(comp = if_else(loc1 == loc2, "Within locality", "Between locality")) %>%
  group_by(comp) %>%
  summarise(
    n_pairs         = n(),
    mean_raw_dist   = mean(dist, na.rm = TRUE),
    median_raw_dist = median(dist, na.rm = TRUE),
    .groups = "drop"
  )

hap_loc <- foc %>%
  count(locality_short, hap_id, name = "n") %>%
  complete(locality_short, hap_id, fill = list(n = 0))

# ---- Figures ----

p_tree <- ggtree(tr, layout = "circular") %<+% seq_tbl +
  geom_tippoint(aes(color = species_label), size = 1.6, alpha = 0.9) +
  labs(color = "Species") +
  theme(legend.position = "right")

ggsave(fig_tree, p_tree, width = 10, height = 10, dpi = 300)

p_hap <- hap_loc %>%
  ggplot(aes(x = locality_short, y = hap_id, fill = n)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(n > 0, as.character(n), "")), size = 3) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
  labs(x = NULL, y = "Flavastellus haplotype", fill = "N") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(fig_hap, p_hap, width = 8, height = 5.5, dpi = 300)

# ---- Write outputs ----

write_csv(seq_tbl, tab_seq)
write_csv(sp_sum,  tab_sp)
write_csv(foc_loc, tab_foc)
write_csv(foc_within_between, tab_pair)