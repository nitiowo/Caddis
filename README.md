# CaddisRNA

- Trinity scripts and job submission for per-sample transcriptome assembly
- BLASTdb scripts for sodium channel transcript extraction
- TransDecoder array job for ORF prediction on cleaned candidate FASTAs
- R scripts that build representative sodium-channel protein/CDS datasets, align representatives, and summarize amino-acid variation

- You still have to use geneious to filter and clean the FASTAs after extraction

## Script Layout

`scripts/assembly/`

- `run_trinity.sh`
- `submit_trinity_jobs.sh`

`scripts/candidate_search/`

- `run_makeblastdb.sh`
- `run_blastn.sh`
- `run_tblastx.sh`
- `run_tblastn.sh`
- `arrayblast.job`

`scripts/transdecoder/`

- `transdecoder_array.job`

`scripts/variants/`

- `build_nach_dataset.R`
- `build_nach_alignment.R`
- `build_nach_variants.R`
