# COI Phylogeography

## Inputs

- COI barcode FASTA with all samples deduplicated and trimmed (from `data/COI/`)
- GenBank COI barcoding metadata table with species assignments and locality information (from `data/metadata/`)
- Limnephilidae outgroup sequences for rooting the ML tree (from `data/COI/outgroups/`)
- TTX locality data (from `data/ttx/`)

## Scripts

- `build_coi_dataset.R` - assembles the final aligned COI dataset from the master FASTA and metadata, produces alignment and species/locality summary tables
- `coi_rooted_tree.R` - builds an ML tree in IQ-TREE and roots it using the limnephilid outgroups
- `lf_haplotype_network.R` - builds a TCS haplotype network for *L. flavastellus* and summarizes haplotype diversity by locality
- `lf_mantel.R` - tests isolation by geographic distance and isolation by TTX toxicity using Mantel and partial Mantel tests

## Outputs

All outputs are written to subfolders of `coi_phylogeography/outputs/`:

- Aligned COI FASTA and species/locality summary tables (in `outputs/alignments/` and `outputs/tables/`)
- ML tree files (in `outputs/trees/`)
- Haplotype network figure and diversity summary (in `outputs/figures/` and `outputs/haplotypes/`)
- Mantel and partial Mantel test results (in `outputs/tables/`)