#!/bin/bash

#####
# Set default command line options
homopolymer_min_length=30
edit_distance=2

block_size='10M'
max_cores='50%'
invert_match=0

#####
# Parse command line options
#####
usage="USAGE: $(basename $0) [-h] [-l <homopolymer length>] [-e <edit distance>] [-c <max cores>] [-b <block size>] [-i]
Filter reads by fuzzy matching of homopolymers runs of ATGCN of a given length.
Input and Output in the same form as paste - - - -

  where:
    -h Show this help text
    -l Homopolymer length (default: 30)
    -e Edit distance for fuzzy matching of k-mers to reads (default: 2)
    -c The number of cores to utilise in the screening. Specified as an absolute number or a percentage (default: 50%)
    -b The size of input data blocks to pass to each core (default: 10M)
    -i Invert output by reporting reads containing fuzzy matches to the homopolymers"

# parse any command line options to change default values
while getopts ":hl:e:c:b:i" opt; do
  case $opt in
    h) echo "$usage"
       exit
       ;;
    l) homopolymer_min_length=$OPTARG
       ;;
    e) edit_distance=$OPTARG
       ;;
    c) max_cores=$OPTARG
       ;;
    b) block_size=$OPTARG
       ;;
    i) invert_match=1
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

#####
# Now onto the program proper

# build the tre-agrep regex
regex="^@[^\t]+\t+[ATGCN]+?(A{${homopolymer_min_length},}|T{${homopolymer_min_length},}|G{${homopolymer_min_length},}|C{${homopolymer_min_length},}|N{${homopolymer_min_length},}){#${edit_distance}}[ATGCN]+?\t\+"

#####
# Remove the reads matching the regex
#####
treagrep_v='-v'
if [[ $invert_match == 1 ]]; then
  treagrep_v=''
fi

#tre-agrep ${treagrep_v} -${edit_distance} -e "${regex}"
#cmd="parallel --gnu --pipe --blocks ${block_size} --max-procs ${max_cores} 'tre-agrep ${treagrep_v} -${edit_distance} -e \"${regex}\"'"
#echo $cmd >&2

parallel --no-notice --gnu --pipe --blocks $block_size --max-procs $max_cores 'tre-agrep '$treagrep_v' -e "'$regex'"'
