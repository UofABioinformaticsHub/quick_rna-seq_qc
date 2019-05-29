#!/usr/bin/env bash
#SBATCH --job-name quick_rna-seq_qc
#SBATCH --output logs/%x_%j.log
#SBATCH --time 00:20:00
#SBATCH --mem 400M
#SBATCH --cpus-per-task 10

module load \
  Bowtie2/2.2.9-foss-2016b \
  SAMtools/1.9-foss-2016b \
  pigz/2.3.3-foss-2016b

READS=${@}

N_READS=100000
SUBSAMPLE_CMD='head'

quick_rna-seq_qc/quick_rna-seq_qc.sh \
  -d ./output_directory \
  -m ${SLURM_CPUS_PER_TASK:-1} \
  -s ${SLURM_CPUS_PER_TASK:-1} \
  -n ${N_READS} \
  -c ${SUBSAMPLE_CMD} \
  -r quick_rna-seq_qc/databases/wheat/wheat_rRNA.fasta \
  -o quick_rna-seq_qc/databases/plantago/plantago_organelles.fasta \
  ${READS}
