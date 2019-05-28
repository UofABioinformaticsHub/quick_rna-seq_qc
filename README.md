# Obtain the code

```bash
git clone --recursive https://github.com/UofABioinformaticsHub/quick_rna-seq_qc
```

# Setup the Environment

For the script to run, it is expected that the following tools be available on the `PATH`.

## Modules on Phenix

```bash
module load \
  Bowtie2/2.2.9-foss-2016b \
  SAMtools/1.9-foss-2016b \
  pigz/2.3.3-foss-2016b \
  parallel/20180922-foss-2016b
```

## Create Bowtie2 Indexes

```bash
# Decompress the FASTA files
find quick_rna-seq_qc/databases -name "*.fasta.gz" -exec pigz --decompress --keep --processes 2 {} +

# Bowtie2 Index the FASTA files
find quick_rna-seq_qc/databases -name "*.fasta" -print0 \
  | xargs -0 -I"{}" bowtie2-build --threads 4 "{}" "{}"
```

# Test for Extent of rRNA Contamination

For a faster run, we subsample 10,000 reads (`-n 10000`) from the file(s) by just looking at the head (`-c head`) of the file(s). We specify the
bowtie2 index base of just rRNA sequences, other tests will be skipped:

```bash
quick_rna-seq_qc/quick_rna-seq_qc.sh \
  -d ./output_directory \
  -m 1 \
  -s 1 \
  -n 10000 \
  -c head \
  -r quick_rna-seq_qc/databases/wheat_rRNA.fasta \
  my.fastq.gz
```

For a more robust analysis we will subsample 100,000 reads (`-n 100000`) randomly (`-c shuf`) from the file(s) and also include a test
of mapping against organelle genome sequences:

```bash
quick_rna-seq_qc/quick_rna-seq_qc.sh \
  -d ./output_directory \
  -m 1 \
  -s 1 \
  -n 100000 \
  -c shuf \
  -r quick_rna-seq_qc/databases/wheat_rRNA.fasta \
  -o quick_rna-seq_qc/databases/plantago/plantago_organelles.fasta \
  my.fastq.gz
```
