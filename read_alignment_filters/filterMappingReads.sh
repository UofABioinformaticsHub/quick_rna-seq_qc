#!/bin/bash
# Example: filterMappingReads.sh \
# -f FASTA \
# -r /mnt/storage/storage1/Nathan/BuJunsWheatMiRNA/fasta/E.fas \
# -l /mnt/storage/storage1/wheat_resources/RNA-seq_contaminants/wheat_RNA.fasta \
# -m E.RNA_mapped.fastq \
# -u E.RNA_unmapped.fastq \
# -b /usr/local/Programs/samtools/samtools-0.1.19 -b /usr/local/Programs/bowtie2/bowtie2-2.0.6 -b /home/nhaigh/bioinf/repos/common-code/converters

#####
# Set default command line options
se_read_files=''
pe1_read_files=''
pe2_read_files=''
read_format='FASTQ'
index=''
mapped_reads_file='/dev/null'
unmapped_reads_file='/dev/stdout'
unmapped_orphaned_reads_file='/dev/null'
bin_dirs=()
max_edit_distance=2
threads=10

# Bowtie2 parameters
# TODO expose more of the bowtie2 parameters to the user
read_gap_penalties='999,999'
ref_gap_penalties='999,999'

####
# Parse command line options
#####
usage="USAGE: $(basename $0) [-h] [-f <read_format>] {-i <SE_read_files> | -1 <PE1_read_files> -2 <PE2_read_files> } -r <reference_sequence_index> [-e <max_alignment_edit_distance>] [-t <number_threads>] [-m <mapped_reads_file>] [-u <unmapped_reads_file>] [-p <dir> ... ]
Sorts reads into different files depending on whether they map to the reference sequences or not.

  where:
  -h Show this help text
  -f Format of reads [FASTQ|FASTA] (default: FASTQ)
  -i Single end input read filename(s) - comma-separated list for multiple files
  -1 #1 mates, paired with files in -2
  -2 #2 mates, paired with files in -1
  -r Reference sequence index file against which reads will be screened
  -e Maximum edit distance allowed in the Bowtie2 alignments (default: 2)
  -t Number of threads to use in Bowtie alignment (default: 10)
  -m Output filename for mapped reads (default: /dev/null)
  -u Output filename for unmapped reads (default: /dev/stdout)
  -o Output filename for unmapped orphaned paired reads (default: /dev/null)
  -p Additional directory to add to the PATH for the executables: bowtie2, samtools (v1.0 or later)"

while getopts ":hf:i:1:2:r:e:t:m:u:o:p:" opt; do
  case $opt in
    h) echo "$usage"
       exit
      ;;
    f) read_format=$OPTARG
      ;;
    i) se_read_files=$OPTARG
      ;;
    1) pe1_read_files=$OPTARG
      ;;
    2) pe2_read_files=$OPTARG
      ;;
    r) index=$OPTARG
      ;;
    e) max_edit_distance=$OPTARG
      ;;
    t) threads=$OPTARG
      ;;
    m) mapped_reads_file=$OPTARG
      ;;
    u) unmapped_reads_file=$OPTARG
      ;;
    o) unmapped_orphaned_reads_file=$OPTARG
      ;;
    p) bin_dirs+=($OPTARG)
      ;;
    ?) printf "Illegal option: '-%s'\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      echo "$usage" >&2
      exit 1
      ;;
  esac
done

if [[ "${index}" == "" ]]; then
  echo -e "No reference provided using -r\n"
  echo "$usage"
  exit
fi
if [[ "${se_read_files}" != "" && "${pe1_read_files}" != "" && "${pe2_read_files}" != "" ]]; then
  echo -e "Can't supply both single-end and paired-ends at the same time\n"
  echo "$usage"
  exit
fi
if [[ "${se_read_files}" != "" &&  "${unmapped_orphaned_reads_file}" != "/dev/null" ]]; then
  echo -e "Can't use -o with -i since unmapped single-end reads will go into the file specified by -u\n"
  echo "$usage"
  exit
fi
if [[ "${unmapped_orphaned_reads_file}" != "/dev/null" && ("${unmapped_orphaned_reads_file}" == "${unmapped_reads_file}" || "${unmapped_orphaned_reads_file}" == "${mapped_reads_file}") ]]; then
  echo -e "Output specified using -o must be different to that of -u and -m\n"
  echo "$usage"
  exit
fi

case $read_format in
  FASTQ) read_format='-q'
    ;;
  FASTA) read_format='-f'
    ;;
  *) echo "Unknown format: ${read_format}"
    exit 1
    ;;
esac

# Add bin_dirs to the front of the path
PATH=$(echo ${bin_dirs[*]} | tr ' ' ':')":$PATH"

# Use bowtie2 for mapping
#CMD="bowtie2 --threads ${threads} -N 1 -L 10 -i S,1,0.25 -D 30 -R 4 ${read_format} --mp 2 --np 2 --ma 0 --score-min 'L,-'$[max_edit_distance*2]',0' --rdg ${read_gap_penalties} --rfg ${ref_gap_penalties} -x ${index}"
CMD="bowtie2 --threads ${threads} -N 1 -L 10 -i S,1,0.25 -D 30 -R 4 ${read_format} -x ${index}"

if [[ "${se_read_files}" != "" ]]; then
  CMD="${CMD} -U ${se_read_files}"
  #echo "SE CMD: ${CMD}"
  # Write single ends to the mapped and unmapped files
  {
    ${CMD} \
      | samtools calmd -S - ${index} \
      | tee \
        >(samtools view -hSF 260 - | samtools bam2fq -s ${mapped_reads_file} -0 ${mapped_reads_file} -) \
        | samtools view -hSf 4 - | samtools bam2fq -s ${unmapped_reads_file} -0 ${unmapped_reads_file} -
  } &> /dev/stderr

fi
if [[ "${pe1_read_files}" != "" && "${pe2_read_files}" != "" ]]; then
  CMD="${CMD} -1 ${pe1_read_files} -2 ${pe2_read_files}"
  echo "PE CMD: ${CMD}"
  # Write paired ends to the mapped and unmapped files and write unmapped orphans to the SE file
  {
    ${CMD} \
      | samtools calmd -S - ${index} \
      | tee \
          >(samtools view -hSF 260 - | samtools bam2fq -s /dev/null - > ${mapped_reads_file}) \
          | samtools view -hSf 4 - | samtools bam2fq -s ${unmapped_orphaned_reads_file} - > ${unmapped_reads_file}
  } 2> /dev/null
fi

