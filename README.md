# Obtain the code

```bash
git clone --recursive https://github.com/UofABioinformaticsHub/quick_rna-seq_qc
cd quick_rna-seq_qc
```

# Running the Quick RNA-Seq QC Analysis

By default the supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script will take the first `100,000` reads from each FASTQ file
and map those reads to the wheat rRNA sequences and the plantago chloroplast genome. The output generated will provide a summary of the
number and percentage of reads found to map to each of these sequences for each file processed.

```bash
# Submit the database preparation script to Slurm
JOBID=$( sbatch scripts/prep.sbatch quick_rna-seq_qc/databases | cut -f4 -d " " )

# Submit slurm jobs for mapping raw reads against the wheat rRNA sequences
# and the plantago chloroplast genome
sbatch --dependency afterok:${JOBID} scripts/quick_qc.sbatch \
  *.fastq.gz
```

# Modifying the Quick RNA-Seq QC Tests

## Random Sampling of Reads

The supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script will take the first `100,000` reads from each FASTQ file. To randomly
select reads from throughout the entire file, simply modify the `SUBSAMPLE_CMD='head'` line to read `SUBSAMPLE_CMD='shuf'`.

Since `shuf` requires reading an entire FASTQ file from start to end, it's use can significantly increase runtimes compared with `head`.

## Subsample More/Less Reads

The supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script will take the first `100,000` reads from each FASTQ file. To modify the
number of reads subsampled, simply modify the `N_READS=100000` line to a different number e.g. `N_READS=1000` or `N_READS=500000`.

## Test Proportion of Reads Derived from ERCC92 Spike-in Sequences

Simply modify the call to `quick_rna-seq_qc/quick_rna-seq_qc.sh` found in the supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script
to include the `-e` argument and specify the location of the [`ERCC92.fasta`](databases/common/ERCC92.fasta) file. e.g. `-e databases/common/ERCC92.fasta`.

## Test Proportion of Reads Derived from Coding Sequences

Simply modify the call to `quick_rna-seq_qc/quick_rna-seq_qc.sh` found in the supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script
to include the `-g` argument and specify the location of the FASTA file containing coding sequences (CDS).
e.g. `-g databases/wheat/ta_IWGSC_MIPSv2.1_HCS_REPR_CDS_2013Nov28.fasta`.

## Report Total Number of Reads

Simply modify the call to `quick_rna-seq_qc/quick_rna-seq_qc.sh` found in the supplied [`scripts/quick_qc.sbatch`](scripts/quick_qc.sbatch) script
to include the `-t` option.

# Robust Quick RNA-Seq QC

For a robust Quick RNA-Seq QC analysis, a command similar to the following could be used:

```bash
quick_rna-seq_qc.sh \
  -d ./output_directory \
  -m 2 \
  -s 2 \
  -n 100000 \
  -c shuf \
  -r databases/wheat/wheat_rRNA.fasta \
  -o databases/plantago/plantago_organelles.fasta \
  -e databases/common/ERCC92.fasta \
  -g databases/wheat/ta_IWGSC_MIPSv2.1_HCS_REPR_CDS_2013Nov28.fasta \
  -t \
  *.fastq.gz \
> report.tab
```
