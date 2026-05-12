# Nav Extraction

Assembly, sodium-channel candidate recovery, manual filtering, ORF prediction, clustering, and representative-sequence selection.

## Layout

- `01_assembly/`: Trinity assembly submission scripts
- `02_candidate_search/`: BLAST database and sodium-channel candidate search scripts
- `03_manual_filtering/`: FASTA cleanup and formatting scripts for manually reviewed candidate hits
- `04_transdecoder/`: TransDecoder and CD-HIT scripts for ORF prediction and clustering
- `05_representative_sequences/`: R scripts that choose one representative Nav protein/CDS per sample

## Scripts

- `01_assembly/submit_trinity_jobs.sh` - submits per-sample Trinity assemblies as an array job on the HPC
- `02_candidate_search/run_makeblastdb.sh` - builds a BLAST database from the sodium-channel reference sequences
- `02_candidate_search/run_blastn.sh` - nucleotide BLAST of assembled contigs against the reference database
- `02_candidate_search/run_tblastx.sh` / `run_tblastn.sh` - translated BLAST searches for more distant hits
- `02_candidate_search/arrayblast.job` - UGE array job script for the BLAST searches
- `03_manual_filtering/spec_sampID_rename.sh` - standardizes sample ID prefixes in FASTA headers after manual review
- `03_manual_filtering/fasta_sortandclean.py` - sorts and cleans FASTA files after manual filtering
- `04_transdecoder/transdecoder_array.job` - runs TransDecoder ORF prediction as a HPC array job
- `04_transdecoder/cdhit.sh` - clusters TransDecoder peptide predictions with CD-HIT to reduce redundancy
- `05_representative_sequences/build_nach_dataset.R` - selects one representative Nav protein and CDS per sample based on ORF quality and coverage
- `05_representative_sequences/build_nach_alignment.R` - aligns the representative Nav proteins
- `05_representative_sequences/build_nach_variants.R` - builds a protein variant matrix from the trimmed alignment

## Outputs for downstream analyses

The key outputs used by `nav_selection/` are written by the `05_representative_sequences/` scripts:

- Representative Nav protein sequences for all samples
- Representative Nav CDS sequences for all samples
- Metadata table linking each representative sequence to its sample ID and ORF classification