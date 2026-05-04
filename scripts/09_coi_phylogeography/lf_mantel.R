library(tidyverse)
library(ape)
library(vegan)

aln_fa   <- "09_coi_phylogeography/outputs/alignments/coi_core_alignment.fasta"
seq_meta <- "09_coi_phylogeography/outputs/tables/coi_sequence_metadata_table.csv"
loc_sum  <- "09_coi_phylogeography/outputs/tables/coi_flavastellus_locality_summary.csv"

tab_dir  <- "09_coi_phylogeography/outputs/tables"
fig_dir  <- "09_coi_phylogeography/outputs/figures"

for (d in c(tab_dir, fig_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

tab_mantel  <- file.path(tab_dir, "lf_mantel_results.csv")
tab_locdist <- file.path(tab_dir, "lf_locality_distance_matrices.csv")
fig_ibd     <- file.path(fig_dir, "lf_isolation_by_distance.png")
fig_ibt     <- file.path(fig_dir, "lf_isolation_by_toxicity.png")

# ---- Load ----

seq_tbl  <- read_csv(seq_meta, show_col_types = FALSE)
dna      <- read.dna(aln_fa, format = "fasta")
loc_data <- read_csv(loc_sum, show_col_types = FALSE)

lf     <- seq_tbl %>% filter(species_label == "Limnephilus flavastellus")
lf_dna <- dna[lf$seq_name, ]
d_ind  <- as.matrix(dist.dna(lf_dna, model = "raw", pairwise.deletion = TRUE))

loc_n      <- table(lf$locality_short)
keep_locs  <- names(loc_n)[loc_n >= 3]

# ---- Per-locality genetic distance (mean pairwise raw distance) ----

n_loc   <- length(keep_locs)
gen_mat <- matrix(0, n_loc, n_loc, dimnames = list(keep_locs, keep_locs))
for (a in seq_len(n_loc)) {
  for (b in seq_len(n_loc)) {
    i <- which(lf$locality_short == keep_locs[a])
    j <- which(lf$locality_short == keep_locs[b])
    gen_mat[a, b] <- mean(d_ind[i, j])
  }
}

# ---- Geographic distance (great-circlenkm) ----

haversine <- function(lat1, lon1, lat2, lon2) {
  R    <- 6371
  dLat <- (lat2 - lat1) * pi / 180
  dLon <- (lon2 - lon1) * pi / 180
  a    <- sin(dLat / 2)^2 +
    cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dLon / 2)^2
  2 * R * asin(sqrt(a))
}

loc_coords <- lf %>%
  filter(locality_short %in% keep_locs) %>%
  group_by(locality_short) %>%
  summarise(lat = mean(lat), lon = mean(lon), .groups = "drop") %>%
  left_join(loc_data %>% select(locality_short, ttx = ttx_skin_mg),
            by = "locality_short") %>%
  arrange(match(locality_short, keep_locs))

geo_mat <- outer(seq_len(n_loc), seq_len(n_loc), Vectorize(function(i, j) {
  haversine(loc_coords$lat[i], loc_coords$lon[i],
            loc_coords$lat[j], loc_coords$lon[j])
}))
dimnames(geo_mat) <- list(keep_locs, keep_locs)

ttx_mat <- as.matrix(dist(loc_coords$ttx))
dimnames(ttx_mat) <- list(keep_locs, keep_locs)

# ---- Mantel tests ----

set.seed(42)
m_geo     <- mantel(as.dist(gen_mat), as.dist(geo_mat), permutations = 9999, method = "spearman")
m_ttx     <- mantel(as.dist(gen_mat), as.dist(ttx_mat), permutations = 9999, method = "spearman")
m_geo_par <- mantel.partial(as.dist(gen_mat), as.dist(ttx_mat), as.dist(geo_mat),
                            permutations = 9999, method = "spearman")
m_ttx_par <- mantel.partial(as.dist(gen_mat), as.dist(geo_mat), as.dist(ttx_mat),
                            permutations = 9999, method = "spearman")

mantel_tbl <- tribble(
  ~test,                                          ~rho,              ~p_value,           ~n_perm,
  "Genetic vs geographic (IBD)",                  m_geo$statistic,   m_geo$signif,       m_geo$permutations,
  "Genetic vs toxicity (IBT)",                    m_ttx$statistic,   m_ttx$signif,       m_ttx$permutations,
  "Genetic vs toxicity | geography (partial)",    m_geo_par$statistic, m_geo_par$signif, m_geo_par$permutations,
  "Genetic vs geography | toxicity (partial)",    m_ttx_par$statistic, m_ttx_par$signif, m_ttx_par$permutations
)
write_csv(mantel_tbl, tab_mantel)

# ---- Pairwise distance long table ----

idx      <- which(upper.tri(gen_mat), arr.ind = TRUE)
pair_tbl <- tibble(
  loc1     = keep_locs[idx[, 1]],
  loc2     = keep_locs[idx[, 2]],
  gen_dist = gen_mat[idx],
  geo_km   = geo_mat[idx],
  ttx_diff = ttx_mat[idx]
)
write_csv(pair_tbl, tab_locdist)

# ---- Figures ----

p_ibd <- ggplot(pair_tbl, aes(geo_km, gen_dist)) +
  geom_point(size = 2.5, color = "#2c7fb8") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(x = "Geographic distance (km)", y = "Mean pairwise COI distance (raw)",
       title = "LF isolation by distance") +
  theme_minimal(base_size = 11)
ggsave(fig_ibd, p_ibd, width = 5.5, height = 4.5, dpi = 300)

p_ibt <- ggplot(pair_tbl, aes(ttx_diff, gen_dist)) +
  geom_point(size = 2.5, color = "#d7301f") +
  geom_smooth(method = "lm", se = TRUE, color = "grey30") +
  labs(x = "Pairwise difference in nearest TTX (mg/cm2)",
       y = "Mean pairwise COI distance (raw)",
       title = "LF isolation by toxicity") +
  theme_minimal(base_size = 11)
ggsave(fig_ibt, p_ibt, width = 5.5, height = 4.5, dpi = 300)
