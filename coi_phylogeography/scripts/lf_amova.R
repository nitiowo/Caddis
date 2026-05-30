# YY amova code
library(tidyverse)
library(ape)
library(adegenet)
library(pegas)
library(poppr)

aln_fa   <- "09_coi_phylogeography/outputs/alignments/coi_core_alignment.fasta"
seq_meta <- "09_coi_phylogeography/outputs/tables/coi_sequence_metadata_table.csv"

tab_dir  <- "09_coi_phylogeography/outputs/tables"
dist_dir <- "09_coi_phylogeography/outputs/distances"
fig_dir  <- "09_coi_phylogeography/outputs/figures"

for (d in c(tab_dir, dist_dir, fig_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# Output files
tab_amova <- file.path(tab_dir, "lf_amova_table.csv")
tab_phist <- file.path(tab_dir, "lf_pairwise_phist.csv")
fig_phist <- file.path(fig_dir, "lf_pairwise_phist_heatmap.png")

# ---- Load LF sequences ----

seq_tbl <- read_csv(seq_meta, show_col_types = FALSE)
dna     <- read.dna(aln_fa, format = "fasta")

# Keep only Limnephilus flavastellus and locality info
lf_meta <- seq_tbl %>%
  filter(species_label == "Limnephilus flavastellus") %>%
  select(seq_name, locality_short)

lf_dna <- dna[lf_meta$seq_name, ]

d_lf <- dist.dna(lf_dna, model = "raw", pairwise.deletion = TRUE)

# Convert the alignment into a genind object and label each sequence by locality
gid <- DNAbin2genind(lf_dna)
strata(gid) <- data.frame(locality = lf_meta$locality_short)
setPop(gid) <- ~locality

# ---- Run AMOVA ----

amv      <- poppr.amova(gid, ~locality, dist = d_lf, within = FALSE, quiet = TRUE)
amv_test <- randtest(amv, nrepet = 999)

# Turn the AMOVA output into a plain table
# Sigma is the variance component and % is the share of total variance
amv_tab    <- as_tibble(amv$results, rownames = "source") %>%
  mutate(
    sigma2       = amv$componentsofcovariance$Sigma,
    pct_variance = amv$componentsofcovariance$`%`
  )

# This AMOVA has one grouping level, so there is one global PhiST value
global_phi <- as.numeric(amv$statphi$Phi[1])
amv_tab$phi <- c(global_phi, NA_real_, NA_real_)[seq_len(nrow(amv_tab))]
write_csv(amv_tab, tab_amova)

# ---- Calculate pairwise PhiST ----

# Calculate pairwise PhiST by calculating AMOVA for each locality pair
pair_phist <- function(d, grp, a, b, nperm = 999) {
  idx   <- which(grp %in% c(a, b))
  d_sub <- as.dist(as.matrix(d)[idx, idx])
  g_sub <- droplevels(factor(grp[idx]))

  amv_pair <- pegas::amova(d_sub ~ g_sub, nperm = nperm)
  phi_pair <- pegas::getPhi(amv_pair$varcomp[, "sigma2"])[1, 1]
  pval <- amv_pair$varcomp[1, "P.value"]

  list(phi = unname(phi_pair), p = unname(pval))
}

# Keep only localities with at least three sequences for the pairwise tests
loc_n      <- table(lf_meta$locality_short)
keep_locs  <- names(loc_n)[loc_n >= 3]
pairs      <- combn(keep_locs, 2, simplify = FALSE)

# Run the pairwise PhiST function for every locality pair
phist_tbl <- map_dfr(pairs, function(pp) {
  res <- pair_phist(d_lf, lf_meta$locality_short, pp[1], pp[2])
  tibble(loc1 = pp[1], loc2 = pp[2], phi_st = res$phi, p_value = res$p)
})
write_csv(phist_tbl, tab_phist)

# ---- Plot the PhiST heatmap ----

# Build square matrix
phist_mat <- matrix(NA_real_, length(keep_locs), length(keep_locs),
                    dimnames = list(keep_locs, keep_locs))
for (i in seq_len(nrow(phist_tbl))) {
  phist_mat[phist_tbl$loc1[i], phist_tbl$loc2[i]] <- phist_tbl$phi_st[i]
  phist_mat[phist_tbl$loc2[i], phist_tbl$loc1[i]] <- phist_tbl$phi_st[i]
}

diag(phist_mat) <- 0

# Turn the matrix into long format for ggplot
phist_long <- as_tibble(phist_mat, rownames = "loc1") %>%
  pivot_longer(-loc1, names_to = "loc2", values_to = "phi_st") %>%
  mutate(
    loc1 = factor(loc1, levels = keep_locs),
    loc2 = factor(loc2, levels = keep_locs)
  )

p_heat <- ggplot(phist_long, aes(loc1, loc2, fill = phi_st)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = formatC(phi_st, digits = 2, format = "f")), size = 3) +
  scale_fill_gradient2(low = "#2c7fb8", mid = "white", high = "#d7301f",
                       midpoint = 0, limits = c(-0.1, 1)) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = "PhiST",
       title = "LF pairwise PhiST among localities (n >= 3)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(fig_phist, p_heat, width = 6.5, height = 6, dpi = 300)
