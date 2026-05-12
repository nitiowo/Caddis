# CaddisRNA

Caddisfly sodium-channel RNA-seq and COI phylogeography project.

Investigating *Limnephilus flavastellus* sodium channel evolution relative to related caddisflies sampled across the PNW and ancestral outgroups to identify candidate large-effect mutations for TTX resistance.

## Repository Layout

- `nach_extraction/`: transcriptome assembly, sodium-channel candidate recovery, manual filtering, ORF prediction, and representative Nav sequence selection
- `nav_selection/`: sodium-channel selection analyses from representative Nav sequences.
- `coi_phylogeography/`: COI barcode dataset construction and phylogeographic analyses.

## Analysis Flow

1. Generate representative Nav sequences from transcriptome files with `nach_extraction/`
2. Test sodium-channel selection patterns with `nav_selection/`
3. Analyze COI lineage structure and geographic context with `coi_phylogeography/`

## Software

CRAN and Bioconductor packages including `Biostrings`, `msa`, `DECIPHER`, `ape`, `adegenet`, `poppr`, `vegan`.

Model-based selection using IQ-TREE, HyPhy, and codeml/PAML uses the conda environment in `nav_selection/env/nav-models.yml`.
