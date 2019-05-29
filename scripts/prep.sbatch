#!/usr/bin/env bash
#SBATCH --job-name prep_quick_rna-seq_qc
#SBATCH --output logs/%x_%j.log
#SBATCH --time 00:10:00
#SBATCH --mem 500M
#SBATCH --cpus-per-task 2

module load \
  Bowtie2/2.2.9-foss-2016b \
  pigz/2.3.3-foss-2016b

DIRS=${@}

# Decompress the FASTA files
find ${DIRS} -name "*.fasta.gz" -exec pigz --decompress --force --keep --processes ${SLURM_CPUS_PER_TASK:-1} {} +

# Bowtie2 Index the FASTA files
find ${DIRS} -name "*.fasta" -print0 \
  | xargs -0 -I"{}" bowtie2-build --threads ${SLURM_CPUS_PER_TASK:-1} "{}" "{}"
