#!/bin/bash

#####
# Set default command line options
kmer=15
kmer_offset=3
edit_distance=0
rev_comp=1
block_size='10M'
max_cores='50%'
invert_match=0
verbose=0
contaminants=()

#####
# Parse command line options
#####
usage="USAGE: $(basename $0) [-h] [-s <sequence>] [-k <k-mer length>] [-o <k-mer offset>] [-e <edit distance>] [-c <max cores>] [-b <block size>] [-R] [-i] [-v]
Filter out reads containing fuzzy matching of k-mers generated using a sliding window over a set of sequences (e.g. adapters).
Sequences are reverse complimented before k-mers are generated. Input and Output in the same form as paste - - - -

  where:
    -h Show this help text
    -s Sequence that may be present in the reads (default: AATGATACGGCGACCACCGAGA TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ATCTCGTATGCCGTCTTCTGCTTG)
    -k Generate k-mers of this length (default: 15)
    -o Sliding window offset (default: 3)
    -e Edit distance for fuzzy matching of k-mers to reads (default: 0)
    -c The number of cores to utilise in the screening. Specified as an absolute number or a percentage (default: 50%)
    -b The size of input data blocks to pass to each core (default: 10M)
    -R Switch to disable reverse compliment of the sequences specified by -s
    -i Invert output by reporting reads containing fuzzy matches to the k-mers
    -v Verbose output"

# parse any command line options to change default values
while getopts ":hs:k:o:e:c:b:Riv" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    s) contaminants+=($OPTARG)
       ;;
    k) kmer=$OPTARG
       ;;
    o) kmer_offset=$OPTARG
       ;;
    e) edit_distance=$OPTARG
       ;;
    c) max_cores=$OPTARG
       ;;
    b) block_size=$OPTARG
       ;;
    R) rev_comp=0
       ;;
    i) invert_match=1
       ;;
    v) verbose=1
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

# Set default list of sequences to search for in the event the user didn't specify any
if [[ ${#contaminants[@]} == 0 ]]; then
	# sequences from https://wikis.utexas.edu/display/GSAF/Illumina+-+all+flavors
	#  Single Index adapter design on a standard Illumina HiSeq or MiSeq
	contaminants=(
	  'AATGATACGGCGACCACCGAGA'                # P5 primer/capture site
	  'TCTACACTCTTTCCCTACACGACGCTCTTCCGATCT'  # TruSeq Read1 primer site - is complementary to Read2 primer site
	  'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'    # Read2 primer site
	  'ATCTCGTATGCCGTCTTCTGCTTG'              # P7 primer/capture site
	)
fi

#####

#####
# Now onto the program proper

if [[ $rev_comp == 1 ]]; then
  # Append reverse complemented sequences to the list
  for contaminant in "${contaminants[@]}"; do
    contaminant_rev_comp=`echo $contaminant | rev | tr 'ATGC' 'TACG'`
    contaminants+=($contaminant_rev_comp)
  done
  #for contaminant in "${contaminants[@]}"; do
  #  echo $contaminant >&2
  #done
fi

kmers=()
for contaminant in "${contaminants[@]}"; do
  if [[ $verbose != 0 ]]; then echo "Contaminant: $contaminant" >&2; fi
  for i in $(seq 0 $kmer_offset $[${#contaminant}-$kmer]); do
    #echo "  $i ${contaminant:$i:$kmer}"
    if [[ $verbose != 0 ]]; then echo -n "             " >&2; fi
    if [[ $verbose != 0 ]]; then
      for j in $(seq 1 $i);do echo -n ' ' >&2; done;
      echo "${contaminant:$i:$kmer}" >&2
    fi
    kmers+=(${contaminant:$i:$kmer})
  done
done
# remove redundant k-mers
kmers=($(for each in ${kmers[@]}; do echo $each; done | sort -u))
# stringify the k-mer string to build a regex

if [[ "${edit_distance}" == 0 ]]; then
  # We'll use fgrep for speed
  CMD="fgrep -f <(echo ${kmers[*]} | tr ' ' '\n')"
else
  # we'll use tre-agrep for fuzzy matching
if [[ $verbose != 0 ]]; then echo "NR kmers: $kmers" >&2; fi
# build the tre-agrep regex
# TODO A bug in tre-agrep (https://github.com/laurikari/tre/issues/20) means we can't negate tab characters correctly. We'll use .+? instead as an alternative
  regex="^@.+?\t+[ATGCN]+?("$(IFS='|'; echo "${kmers[*]}")"){#${edit_distance}}[ATGCN]+?\t\+"
  CMD="tre-agrep -e ${regex}"
fi

if [[ $verbose != 0 ]]; then echo "CMD: ${CMD}" >&2; fi
  
#####
#  Correctly invert the matches if requested
#####
if [[ $invert_match == 0 ]]; then
  treagrep_v=''
  CMD+=" -v"
fi

#tre-agrep ${treagrep_v} -${edit_distance} -e "${regex}"
parallel --no-notice --gnu --pipe --blocks $block_size --max-procs $max_cores "${CMD}"

