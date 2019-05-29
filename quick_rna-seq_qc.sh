#!/bin/bash
#####
# Get some quick QC stats using a subset of reads from a bunch of gzipped FASTQ files
#   The subset of reads is generated randomly from the file using shuf
# WARNING: If a file contains reads from multiple lanes, we might miss important QC
#          stats/issues. Calculate number of reads/sample/lane in a file (assumes
#          reads from the same lane are consecutive):
#            pigz -dcp2 in.fastq.gz | paste - - - - | cut -f 1 | cut -f 4 -d':' | uniq -c
# It assumes reads are single-ends
#####

# Usage: quick_rna_seq_qc.sh [output directory] <file1.fastq.gz> [...]
#        mkdir output
#        find ./ -name "*.fastq.gz" -exec /home/nhaigh/git/rna-seq_qc/quick_rna-seq_qc.sh ./output {} \; > quick_rna-seq_qc.out

n_mapping_threads=20
n_sorting_threads=10

# pick a subset of about 10% of the number of reads in a file
n_read_subsample_mapping=100000
n_read_subsample_adapter=100000
n_read_subsample_homopolymer=100000
subset_cmd="shuf"

outdir='./'

rRNA_reference_index_base=''
organelle_reference_index_base=''
cds_reference_index_base=''
ercc_reference_index_base=''
adapter_sequences=()

################################
# Parse command line arguments #
################################
usage="USAGE: $(basename "$0") [-h] [-d output_dir] [-m mapping_threads] [-s sorting_threads] [-n num_reads] [[-a adapter seq] ...] [-r rRNA] [-o organelle] [-e ERCC] <file1.fastq.gz> [[<file2.fastq.gz>]...]

where:
  -h Show this helpful help text
  -d Output directory [Default: ./]
  -m Number of bowtie2 mapping threads [Default: 20]
  -s Number of samtools sorting threads [Default: 10]
  -n Number of reads to randomly subsample for generating metrics [Default: 100,000]
  -c Command for getting subset of reads [default: shuf]
  -a Adapter sequence that may be present in the reads [default: AATGATACGGCGACCACCGAGA TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ATCTCGTATGCCGTCTTCTGCTTG]
  -r Perform rRNA contamination quantification using this rRNA FASTA index file
  -o Perform organelle contamination quantification using this Organelle FASTA index file
  -e Perform ERCC spike-in quantification using this ERCC FASTA index file
  -g Perform gene quantification using this CDS FASTA index file
  -t Display total reads column
"

while getopts ":hd:m:s:n:c:a:r:o:e:g:t" opt; do
  case $opt in
    h) echo "$usage"
       exit;;
    d) outdir=$OPTARG;;
    m) n_mapping_threads=$OPTARG;;
    s) n_sorting_threads=$OPTARG;;
    n) n_read_subsample_mapping=$OPTARG; n_read_subsample_adapter=$OPTARG; n_read_subsample_homopolymer=$OPTARG;;
    c) subset_cmd=$OPTARG;;
    a) adapter_sequences+=($OPTARG);;
    r) rRNA_reference_index_base=$OPTARG;;
    o) organelle_reference_index_base=$OPTARG;;
    e) ercc_reference_index_base=$OPTARG;;
    g) cds_reference_index_base=$OPTARG;;
    t) get_total_reads=1;;
    ?) printf "Illegal option: '-%s'\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1;;
  esac
done
# shift all processed options away with
shift $((OPTIND-1))

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#echo "DIR: ${DIR}"


# Set default list of sequences to search for in the event the user didn't specify any
if [[ ${#adapter_sequences[@]} == 0 ]]; then
  # sequences from https://wikis.utexas.edu/display/GSAF/Illumina+-+all+flavors
  #  Single Index adapter design on a standard Illumina HiSeq or MiSeq
  adapter_sequences=(
    'AATGATACGGCGACCACCGAGA'                # P5 primer/capture site
    'TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # TruSeq Read1 primer site - is complementary to Read2 primer site
    'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'    # Read2 primer site
    'ATCTCGTATGCCGTCTTCTGCTTG'              # P7 primer/capture site
  )
fi


# Output files will go into the cwd unless the first argument is a directory
# TODO change all this behaviors to using getopts
#if [[ "$1" ]]; then
#  outdir="$1"
#  shift
#fi
# Remove any trailing slashes
# BUG: is the root dir is given, it will also be removed
outdir=${outdir%/}
#echo "outdir: ${outdir}"

echo -e "Filename\tAdapter Containing Reads (#)\tAdapter Containing Reads (%)\trRNA Mapping Reads (#)\trRNA Mapping Reads (%)\tOrganelle Mapping Reads (#)\tOrganelle Mapping Reads (%)\tCDS Mapping Reads (#)\tCDS Mapping Reads (%)\tERCC92 Mapping Reads (#)\tERCC92 Mapping Reads (%)\tTotal Reads (#)"

for file in "$@"; do
  #echo >&2 "Processing: ${file}"
  echo -n "${file}"
  
  # Setup output dir for each input file. The dir structure will be the same as the raw input file to avoid clobbering of similarly named files
  file_basename=${file##*/}
  file_path=${file%/*};
  if [[ ! -d "${outdir}/${file_path##*.}" ]]; then
    mkdir -p "${outdir}/${file_path##*.}"
  fi
  
  ##########
  # 1 Adapter-containing reads
  #   requires:
  #     GNU parallel
  #     tre-agrep if ED>0
  ##########
  ED=0
  echo -ne "\t"
  n_adapter_containing=$(pigz -dcp 2 "${file}" \
    | paste - - - - \
    | ${subset_cmd} -n ${n_read_subsample_adapter} \
    | ${DIR}/adapter_filters/fuzzy_sequence_filter.sh -e ${ED} -c '50%' -b '1M' -i $(printf " -i %s" "${adapter_sequences[@]}") \
    | wc -l)
  echo -ne "${n_adapter_containing}\t$(bc <<< "scale=2; ${n_adapter_containing} * 100 / ${n_read_subsample_adapter}")%"
  
  ##########
  # 2 Homopolymer runs
  #   requires:
  #     GNU parallel
  #     tre-agrep if ED>0
  ##########
#  ED=2
#  echo -ne "\t"
#  n_homopolymer_run_containing=$(pigz -dcp 2 "${file}" \
#    | paste - - - - \
#    | ${subset_cmd} -n ${n_read_subsample_homopolymer} \
#    | ${DIR}/homopolymer_filters/fuzzy_homopolymer_filter.sh -l 30 -e ${ED} -c '50%' -b '10M' -i \
#    | wc -l)
#  echo -ne "${n_homopolymer_run_containing}\t$(bc <<< "scale=2; ${n_homopolymer_run_containing} * 100 / ${n_read_subsample_homopolymer}")%"
  
  ##########
  # 3 rRNA-mapping
  ##########
  echo -ne "\t"
  if [ ! -z "${rRNA_reference_index_base}" ]; then
    ${DIR}/read_alignment_filters/filterMappingReads.sh \
      -i <(pigz -dcp 2 "${file}" | paste - - - - | ${subset_cmd} -n ${n_read_subsample_mapping} | tr '\t' '\n') \
      -f FASTQ \
      -t ${n_mapping_threads} \
      -r ${rRNA_reference_index_base} \
      -u "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.rRNA_unmapped.fastq" \
      -m "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.rRNA_mapped.fastq" &> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.rRNA_mapping.log"
    n_rRNA_mapped=$(paste - - - - < "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.rRNA_mapped.fastq" | wc -l)
    echo -ne "${n_rRNA_mapped}\t$(bc <<< "scale=2; ${n_rRNA_mapped} * 100 / ${n_read_subsample_mapping}")%"
  else
    echo -ne "NA\tNA%"
  fi

  ##########
  # 4 Organelle-mapping
  ##########
  echo -ne "\t"
  if [ ! -z "${organelle_reference_index_base}" ]; then
    ${DIR}/read_alignment_filters/filterMappingReads.sh \
      -i <(pigz -dcp 2 "${file}" | paste - - - - | ${subset_cmd} -n ${n_read_subsample_mapping} | tr '\t' '\n') \
      -f FASTQ \
      -t ${n_mapping_threads} \
      -r ${organelle_reference_index_base} \
      -u "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.organelle_unmapped.fastq" \
      -m "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.organelle_mapped.fastq" &> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.organelle_mapping.log"
    n_organelle_mapped=$(paste - - - - < "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.organelle_mapped.fastq" | wc -l)
    echo -ne "${n_organelle_mapped}\t$(bc <<< "scale=2; ${n_organelle_mapped} * 100 / ${n_read_subsample_mapping}")%"
  else
    echo -ne "NA\tNA%"
  fi

  ##########
  # 5 CDS-mapping
  ##########
  echo -ne "\t"
  if [ ! -z "${cds_reference_index_base}" ]; then
    ED=3
    bowtie2 --time --threads ${n_mapping_threads} -N 1 -L 10 -i S,1,0.25 -D 30 -R 4 --ma 0 --mp 2 --np 2 --score-min "L,-$((ED*2)),0" --rdg 999,999 --rfg 999,999 -x ${cds_reference_index_base} -U <(pigz -dcp 2 "${file}" | paste - - - - | ${subset_cmd} -n ${n_read_subsample_mapping} | tr '\t' '\n') 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.CDS_mapping.log" \
      | samtools calmd -S - ${cds_reference_index_base} 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.CDS_calmd.log" \
      | samtools view -buS - 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.CDS_view.log" \
      | samtools sort -@ ${n_sorting_threads} -m 25G -o "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_CDS.bam" - &> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.CDS_sort.log"
    samtools index -c "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_CDS.bam"
    n_cds_mapped=$(samtools view -F 4 "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_CDS.bam" 2> /dev/stderr | cut -f 1 | sort -u --parallel=${n_sorting_threads} -S 25G | wc -l)
    echo -ne "${n_cds_mapped}\t$(bc <<< "scale=2; ${n_cds_mapped} * 100 / ${n_read_subsample_mapping}")%"
  else
    echo -ne "NA\tNA%"
  fi

  ##########
  # 6 ERCC Mix 1 spike-in mapping
  ##########
  echo -ne "\t"
  if [ ! -z "${ercc_reference_index_base}" ]; then
    ED=2
    #${bowtie2_dir}/bowtie2-build --offrate 1 ${reference} ${reference}
    bowtie2 --time --threads ${n_mapping_threads} -N 1 -L 10 -i S,1,0.25 -D 30 -R 4 --ma 0 --mp 2 --np 2 --score-min "L,-$((ED*2)),0" --rdg 999,999 --rfg 999,999 -x ${ercc_reference_index_base} -U <(pigz -dcp 2 "${file}" | paste - - - - | ${subset_cmd} -n ${n_read_subsample_mapping} | tr '\t' '\n') 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.ERCC92_mapping.log" \
      | samtools calmd -S - ${ercc_reference_index_base} 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.ERCC92_calmd.log" \
      | samtools view -buS - 2> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.ERCC92_view.log" \
      | samtools sort -@ ${n_sorting_threads} -m 25G -o "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_ERCC92.bam" - &> "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}.ERCC92_sort.log"
    samtools index -c "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_ERCC92.bam"
    n_ercc_mapped=$(samtools view -F 4 "${outdir}/${file_path##*.}/${file_basename%%.fastq.gz}_vs_ERCC92.bam" 2> /dev/stderr | cut -f 1 | sort -u --parallel=${n_sorting_threads} -S 25G | wc -l)
    echo -ne "${n_ercc_mapped}\t$(bc <<< "scale=2; ${n_ercc_mapped} * 100 / ${n_read_subsample_mapping}")%"
  else
    echo -ne "NA\tNA%"
  fi

  ##########
  # 7 Number of reads
  ##########
  if [ ! -z "${get_total_reads}" ]; then
    n_reads=$(pigz -dcp 2 "${file}" | paste - - - - | wc -l)
    echo -ne "\t${n_reads}"
  else
    echo -ne "\tNA"
  fi
  
  echo -ne "\n"
done
