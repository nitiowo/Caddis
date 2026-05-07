library(tidyverse)
library(jsonlite)

root <- Sys.getenv("CADDISRNA_ROOT", unset = getwd())

raw_dir <- file.path(root, "nav_selection/outputs/model_raw")
out_dir <- file.path(root, "nav_selection/outputs/model_selection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- IQ-TREE ----

iq_lines <- read_lines(file.path(raw_dir, "iqtree/nav_codon.iqtree"))
contree  <- read_lines(file.path(raw_dir, "iqtree/nav_codon.contree"))

grab <- function(pat) {
  hit <- str_match(iq_lines, pat)[, 2]
  hit[!is.na(hit)][1]
}

support_vals <- as.numeric(str_match_all(contree, "\\)([0-9.]+):")[[1]][, 2])

best_model <- grab("Best-fit model according to BIC: ([^\\r\\n]+)")
log_lik    <- as.numeric(grab("Log-likelihood of the tree: (-?[0-9.]+)"))
omega_iq   <- as.numeric(grab("Nonsynonymous/synonymous ratio \\(omega\\): ([0-9.]+)"))
n_bs95     <- sum(support_vals >= 95)

# ---- HyPhy ----

busted <- fromJSON(file.path(raw_dir, "hyphy/busted.json"), simplifyVector = FALSE)
fel    <- fromJSON(file.path(raw_dir, "hyphy/fel.json"),    simplifyVector = FALSE)
meme   <- fromJSON(file.path(raw_dir, "hyphy/meme.json"),   simplifyVector = FALSE)

# FEL site columns: alpha, beta, alpha=beta LRT, p-value
fel_mat  <- do.call(rbind, fel$MLE$content[["0"]])
fel_p    <- as.numeric(fel_mat[, 4])
fel_dir  <- ifelse(as.numeric(fel_mat[, 2]) > as.numeric(fel_mat[, 1]), "pos", "neg")
fel_pos  <- sum(fel_p <= 0.05 & fel_dir == "pos", na.rm = TRUE)
fel_neg  <- sum(fel_p <= 0.05 & fel_dir == "neg", na.rm = TRUE)

# MEME site column 6 is p-value
meme_mat <- do.call(rbind, meme$MLE$content[["0"]])
meme_p   <- as.numeric(meme_mat[, 6])
meme_pos <- sum(meme_p <= 0.05, na.rm = TRUE)

busted_lrt <- as.numeric(busted$`test results`$LRT)
busted_p   <- as.numeric(busted$`test results`$`p-value`)

# ---- codeml branch model ----

alt_lines  <- read_lines(file.path(raw_dir, "codeml_branch/codeml_alt.out"))
null_lines <- read_lines(file.path(raw_dir, "codeml_branch/codeml_null.out"))

grab_lnl <- function(lines) {
  as.numeric(na.omit(str_match(lines, "lnL\\(.*\\):\\s+(-?[0-9.]+)")[, 2])[1])
}
grab_w <- function(lines) {
  raw <- na.omit(str_match(lines, "w \\(dN/dS\\) for branches:\\s+([^\\r\\n]+)")[, 2])[1]
  as.numeric(str_split(str_squish(raw), " ")[[1]])
}

alt_lnl  <- grab_lnl(alt_lines)
null_lnl <- grab_lnl(null_lines)
alt_w    <- grab_w(alt_lines)

lrt   <- 2 * (alt_lnl - null_lnl)
p_val <- pchisq(lrt, df = 1, lower.tail = FALSE)

# ---- Combined summary ----

model_summary <- tribble(
  ~metric,                        ~value,
  "iqtree_best_model_bic",        best_model,
  "iqtree_log_likelihood",        format(log_lik),
  "iqtree_omega",                 format(omega_iq),
  "iqtree_bootstrap_nodes_ge_95", as.character(n_bs95),
  "busted_lrt",                   format(busted_lrt),
  "busted_p_value",               format(busted_p),
  "fel_positive_sites_p_0_05",    as.character(fel_pos),
  "fel_negative_sites_p_0_05",    as.character(fel_neg),
  "meme_episodic_sites_p_0_05",   as.character(meme_pos),
  "codeml_branch_lrt",            format(lrt),
  "codeml_branch_p_value",        format(p_val),
  "codeml_background_omega",      format(alt_w[1]),
  "codeml_foreground_omega",      format(alt_w[2])
)

write_csv(model_summary, file.path(out_dir, "model_selection_summary.csv"))
