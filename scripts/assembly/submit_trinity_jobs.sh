#!/bin/bash
# submit_jobs.sh
# Nitin Vincent 2025
# This script creates and submits a job for each sample. This is necessary when working with RNA reads from different species in each sample.
# If your samples are different tissue types/samples that can be mapped to each other for assembly, include all samples in the same Trinity run - don't use this script.
# See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity for more info
# Change filename formatting as required.
# Use trimmed treads

cd "/temp180/mpfrende/nvincen2/Caddis"  # Working directory should be above the directory where the trinity image is

# Define directories and parameters use absolute paths
INPUT_DIR="/temp180/mpfrende/nvincen2/Caddis/CaddisRNA/trimmed"           # Directory containing your trimmed reads
OUTPUT_BASE="/temp180/mpfrende/nvincen2/Caddis/trinity_indiv_out"     # Base directory for Trinity outputs
WRAPPER_SCRIPT="/temp180/mpfrende/nvincen2/Caddis/trinity_workspace/run_trinity.sh"  # Full path to run_trinity.sh
MEMORY="250G"                        # Memory per job
CPU=8                               # Number of CPU cores per job (specifies to both Trinity and jobscript header)
EMAIL="nvincen2@nd.edu"        # Email for job notifications
QUEUE="largemem"                     # SGE queue to use

# Create the output base directory if needed
mkdir -p "${OUTPUT_BASE}"

# Loop through each *_R1.fq.gz file in the input directory
for left_file in "${INPUT_DIR}"/*_R1_trimmed.fastq.gz; do
  # Derive the sample name from the filename (adjust pattern if needed). $left_file and $right_file are the sample names passed to the wrapper script.
  sample=$(basename "${left_file}" _R1_trimmed.fastq.gz)
  right_file="${INPUT_DIR}/${sample}_R2_trimmed.fastq.gz"
  
  # Verify that the right file exists
  if [[ ! -f "${right_file}" ]]; then
    echo "Paired file for sample ${sample} not found, skipping..."
    continue
  fi

  # Set output directory for this sample
  outdir="${OUTPUT_BASE}/${sample}_trinity"
  mkdir -p "${outdir}"

  # Create a SGE job script for this sample. Everything between the two 'EOF's get put into a job_script.
  job_script="${OUTPUT_BASE}/${sample}_trinity_job.sh"
  cat > "${job_script}" <<EOF
#!/bin/bash
#$ -N Trinity_${sample}           # Job name
#$ -q ${QUEUE}                    # Queue selection
#$ -pe smp ${CPU}   # Parallel environment and CPU cores
#$ -M ${EMAIL}                    # Email for job notifications
#$ -m abe                         # Email on (a)bort, (b)egin, and (e)nd
#$ -o ${outdir}/${sample}_trinity.out  # Standard output log
#$ -e ${outdir}/${sample}_trinity.err  # Standard error log

# Suggested use is by downloading the Trinity Singularity image
# If not using Singularity, load Trinity module if required
# module load trinity

# Run the Trinity wrapper script with the current sample's files and output directory
bash ${WRAPPER_SCRIPT} ${left_file} ${right_file} ${outdir} ${CPU} ${MEMORY}
EOF

  # Make the job script executable (optional)
  # chmod +x "${job_script}"

  # Submit the job script to SGE
  qsub "${job_script}"
  echo "Submitted Trinity job for sample ${sample}"
done