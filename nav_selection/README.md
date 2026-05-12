# Nav Selection

Codon-based multiple sequence alignment, substitution model selection and phylogenetic inference in IQ-TREE, and test for positive selection using HyPhy (BUSTED, FEL, MEME) and PAML codeml on representative NaV sequences. The codeml branch model compares *L. flavastellus* against all other tips.

## Inputs

All inputs come from earlier steps in `nach_extraction/`:

- Representative Nav protein sequences for all samples (from `05_representative_sequences/`)
- Representative Nav CDS sequences for all samples (from `05_representative_sequences/`)
- Representative sequence metadata table with sample IDs and ORF classifications (from `05_representative_sequences/`)
- Reference outgroup NaV protein sequences from the TransDecoder reference set (from `04_transdecoder/`)

## Scripts

- `nav_selection_inputs.R` - builds the protein alignment, back-translates to a codon alignment, and applies an occupancy filter
- `model_inference.sh` - runs IQ-TREE, HyPhy BUSTED/FEL/MEME, and codeml using settings in `config/model_config.sh`
- `model_selection_summary.R` - extracts statistics from IQ-TREE, HyPhy, and codeml outputs into a single CSV
- `pore_window_selection.R` - polarizes Nav substitutions against an outgroup reference, sums derived changes per sample, and tests for enrichment in the TTX pore-window region
- `deep_ortholog_selection.R` - polarization against deep insect orthologs (Musca, Drosophila, Bombyx, Tribolium) and compares derived substitution burden between LF and non-LF cddis species
- `clean_hyphy_alignment.py` - masks stop codons in the codon alignment before passing it to HyPhy
- `prepare_codeml_inputs.py` - cleans and formats the codon alignment for codeml (PHYLIP format)

## Outputs

All outputs are written to `nav_selection/outputs/`.

- Filtered codon alignment used as input to model inference (in `outputs/alignments/`)
- Raw IQ-TREE, HyPhy, and codeml result files (in `outputs/model_raw/`)
- `model_selection_summary.csv` - best-fit model, log-likelihood, omega estimates, and selection test statistics in one table (in `outputs/model_selection/`)
- `ortholog_polarized_substitutions.csv` - per-position substitution table with polarization categories relative to deep insect orthologs (in `outputs/selection/`)
- `outgroup_polarized_substitutions.csv` - per-position substitution table with polarization categories relative to  non-LF caddis outgroup (in `outputs/selection/`)
- Derived substitution burden summaries and group comparison tables (in `outputs/selection/`)
- Pore-window LF vs outgroup contrast table and Fisher test results (in `outputs/selection/`)
