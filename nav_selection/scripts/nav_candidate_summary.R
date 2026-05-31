library(tidyverse)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

candidate_file <- file.path(root, "nav_selection/config/candidate_sites.csv")
site_context_file <- file.path(root, "nav_selection/outputs/selection/nav_site_context.csv")

out_dir <- file.path(root, "nav_selection/outputs/selection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_summary <- file.path(out_dir, "nav_candidate_summary.csv")
out_support <- file.path(out_dir, "nav_candidate_support_long.csv")

# ---- Functions ----

valid_aa <- function(value) {
  !is.na(value) & !value %in% c("", "-", "?", "X")
}

modal_aa <- function(value) {
  value <- value[valid_aa(value)]
  if (length(value) == 0) {
    return(NA_character_)
  }
  names(sort(table(value), decreasing = TRUE))[1]
}

as_bool <- function(value) {
  toupper(as.character(value)) == "TRUE"
}

# ---- Load Inputs ----

if (!file.exists(site_context_file)) {
  stop("Run nav_site_context.R before nav_candidate_summary.R")
}

candidate_tbl <- read_csv(candidate_file, show_col_types = FALSE) %>%
  mutate(aln_col = as.integer(aln_col))

if ("candidate" %in% names(candidate_tbl)) {
  candidate_tbl <- candidate_tbl %>%
    filter(as_bool(candidate))
}

site_context <- read_csv(site_context_file, show_col_types = FALSE)
rep_cols <- names(site_context)[str_starts(names(site_context), "rep_aa_")]

if (length(rep_cols) == 0) {
  stop("No representative amino-acid columns found in nav_site_context.csv")
}

candidate_sites <- candidate_tbl %>%
  left_join(site_context, by = "aln_col", suffix = c("", "_site"))

if (!"candidate_aa" %in% names(candidate_sites)) {
  candidate_sites$candidate_aa <- NA_character_
}

inferred_states <- candidate_sites %>%
  select(aln_col, all_of(rep_cols)) %>%
  pivot_longer(all_of(rep_cols), names_to = "rep_col", values_to = "aa") %>%
  filter(str_detect(rep_col, "_LiFl$")) %>%
  group_by(aln_col) %>%
  summarise(inferred_candidate_aa = modal_aa(aa), .groups = "drop")

candidate_sites <- candidate_sites %>%
  left_join(inferred_states, by = "aln_col") %>%
  mutate(candidate_aa = coalesce(candidate_aa, inferred_candidate_aa))

# ---- Sample Support ----

support_long <- candidate_sites %>%
  select(aln_col, candidate_aa, all_of(rep_cols)) %>%
  pivot_longer(all_of(rep_cols), names_to = "rep_col", values_to = "aa") %>%
  mutate(
    sample_label = str_remove(rep_col, "^rep_aa_"),
    sample_id = str_extract(sample_label, "CF[0-9]+"),
    focal = str_detect(rep_col, "_LiFl$"),
    called = valid_aa(aa),
    supports_candidate = called & aa == candidate_aa
  )

support_summary <- support_long %>%
  group_by(aln_col) %>%
  summarise(
    lifl_called = sum(focal & called),
    lifl_support = sum(focal & supports_candidate),
    non_lifl_called = sum(!focal & called),
    non_lifl_support = sum(!focal & supports_candidate),
    .groups = "drop"
  )

summary_cols <- c(
  "aln_col", "candidate_tier", "candidate_label", "musca_aa_number",
  "dmel_aa_number", "limluna_aa_number", "limluna_aa", "candidate_aa",
  "feature_class", "primary_context", "why_candidate", "codon_alignment_site",
  "fel_p_value", "fel_direction", "fel_significant_p_0_05",
  "meme_p_value", "meme_significant_p_0_05"
)

candidate_summary <- candidate_sites %>%
  select(any_of(summary_cols)) %>%
  left_join(support_summary, by = "aln_col")

if (!"fel_significant_p_0_05" %in% names(candidate_summary)) {
  candidate_summary$fel_significant_p_0_05 <- FALSE
}
if (!"meme_significant_p_0_05" %in% names(candidate_summary)) {
  candidate_summary$meme_significant_p_0_05 <- FALSE
}

candidate_summary <- candidate_summary %>%
  mutate(
    lifl_support_ratio = lifl_support / lifl_called,
    non_lifl_support_ratio = non_lifl_support / non_lifl_called,
    any_hyphy_p_0_05 = coalesce(fel_significant_p_0_05, FALSE) |
      coalesce(meme_significant_p_0_05, FALSE)
  ) %>%
  arrange(aln_col)

write_csv(candidate_summary, out_summary)
write_csv(support_long, out_support)