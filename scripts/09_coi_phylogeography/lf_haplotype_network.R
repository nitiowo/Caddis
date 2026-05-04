library(tidyverse)
library(ape)
library(pegas)

aln_fa   <- "09_coi_phylogeography/outputs/alignments/coi_core_alignment.fasta"
seq_meta <- "09_coi_phylogeography/outputs/tables/coi_sequence_metadata_table.csv"

fig_dir  <- "09_coi_phylogeography/outputs/figures"
tab_dir  <- "09_coi_phylogeography/outputs/tables"
hap_dir  <- "09_coi_phylogeography/outputs/haplotypes"

for (d in c(fig_dir, tab_dir, hap_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

fig_net      <- file.path(fig_dir, "lf_haplotype_network.png")
tab_hap_loc  <- file.path(tab_dir, "lf_haplotype_by_locality.csv")
tab_hap_freq <- file.path(tab_dir, "lf_haplotype_frequencies.csv")
tab_stats    <- file.path(tab_dir, "lf_haplotype_network_stats.csv")
hap_rds      <- file.path(hap_dir, "lf_haplotypes.rds")

# ---- Load and subset to LF ----

seq_tbl <- read_csv(seq_meta, show_col_types = FALSE)
dna     <- read.dna(aln_fa, format = "fasta")

lf_ids  <- seq_tbl %>%
  filter(species_label == "Limnephilus flavastellus") %>%
  pull(seq_name)

lf_dna  <- dna[lf_ids, ]
lf_meta <- seq_tbl %>% filter(seq_name %in% lf_ids)

# ---- Haplotypes ----

haps       <- haplotype(lf_dna)
hap_assign <- attr(haps, "index")
names(hap_assign) <- rownames(haps)

ind_hap <- tibble(
  seq_name = unlist(lapply(seq_along(hap_assign),
                           function(i) rownames(lf_dna)[hap_assign[[i]]])),
  hap      = unlist(lapply(seq_along(hap_assign),
                           function(i) rep(names(hap_assign)[i], length(hap_assign[[i]]))))
) %>%
  left_join(lf_meta %>% select(seq_name, locality_short), by = "seq_name")

hap_freq <- ind_hap %>%
  count(hap, name = "n_individuals") %>%
  arrange(desc(n_individuals))

hap_loc <- ind_hap %>%
  count(hap, locality_short, name = "n") %>%
  pivot_wider(names_from = locality_short, values_from = n, values_fill = 0) %>%
  arrange(hap)

write_csv(hap_freq, tab_hap_freq)
write_csv(hap_loc,  tab_hap_loc)
saveRDS(list(haps = haps, ind_hap = ind_hap), hap_rds)

# ---- Network ----

net        <- haploNet(haps)
loc_levels <- sort(unique(ind_hap$locality_short))
loc_pal    <- setNames(scales::hue_pal()(length(loc_levels)), loc_levels)

hap_x_loc <- hap_loc %>%
  column_to_rownames("hap") %>%
  as.matrix()
hap_x_loc <- hap_x_loc[rownames(haps), loc_levels, drop = FALSE]

png(fig_net, width = 2400, height = 2000, res = 300)
plot(
  net,
  size           = attr(net, "freq"),
  pie            = hap_x_loc,
  bg             = loc_pal[loc_levels],
  legend         = FALSE,
  show.mutation  = 1,
  scale.ratio    = 1,
  threshold      = 0
)
legend("topleft", legend = loc_levels, fill = loc_pal[loc_levels], bty = "n", cex = 0.8)
title(main = "Limnephilus flavastellus COI haplotype network", line = 0.5)
dev.off()

# ---- Diversity stats ----

n_ind  <- nrow(lf_dna)
n_hap  <- nrow(haps)
hd     <- (n_ind / (n_ind - 1)) * (1 - sum((hap_freq$n_individuals / n_ind)^2))
pi_val <- nuc.div(lf_dna, pairwise.deletion = TRUE)
td     <- tajima.test(lf_dna)

write_csv(
  tibble(
    metric = c("n_individuals", "n_haplotypes", "haplotype_diversity", "nucleotide_diversity",
               "tajimas_D", "tajimas_D_pval_normal", "tajimas_D_pval_beta"),
    value  = c(n_ind, n_hap, hd, pi_val, td$D, td$Pval.normal, td$Pval.beta)
  ),
  tab_stats
)
