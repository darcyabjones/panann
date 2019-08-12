#!/usr/bin/env bash

set -euo pipefail

# Default temp directory
TMPDIR="tmp$$"
NCPU="$(grep -c '^processor' /proc/cpuinfo)"
OUTFILE="/dev/stdout"

usage() {
  echo 'Runs spaln in parallel within regions specified by a BED file.

USAGE: exonerate_parallel.sh -g GENOME -p PROTEINS -b BED -n NCPUS -t TMPDIR > results.gff

Arguments:
-g -- The genome in fasta format [REQUIRED].
-p -- The proteins in fasta format [REQUIRED if not -q].
-q -- The proteins in tsv format [REQUIRED if not -p]
-b -- A bed file with regions to search in the first 3 columns
      and a comma separated list of protein ids to search against that
      region in the 4th column [REQUIRED].
-n -- How many concurrent tasks to run [DEFAULT: All detected cpus].
-t -- Where to store temporary files [DEFAULT: "tmp$$"].
-o -- Where to write the output gff to [DEFAULT: stdout].


Requires the following on PATH:
- spaln_region.sh
- tsv_to_fasta.sh fasta_to_tsv.sh
- spaln
- standard linux utilities: awk, xargs, sed, tr, grep, bash


Example:
```
spaln_parallel.sh -g genome.fasta -p proteins.fasta -b merged.bed -n 16 -o results.gff
```
'
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":hg:p:q:b:n:t:o:" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    g ) GENOME="${OPTARG}" ;;
    p ) PROTEINS="${OPTARG}" ;;
    q ) PROTEINS_TSV="${OPTARG}" ;;
    b ) BED="${OPTARG}" ;;
    n ) NCPU="${OPTARG}" ;;
    t ) TMPDIR="${OPTARG}" ;;
    o ) OUTFILE="${OPTARG}" ;;
    : )
      usage
      echo "ERROR: Invalid option: -${OPTARG} requires an argument" 1>&2
      exit 1
      ;;
    \? )
      usage
      echo "ERROR: Invalid option: -${OPTARG}" 1>&2
      exit 1
      ;;
  esac
done


if   [ -z "${GENOME-}" ] || [ -z "${BED-}" ]
then
  usage
  echo "ERROR: One or more required arguments have not been provided"
  exit 1
fi

if [ -z "${PROTEINS-}" ] && [ -z "${PROTEINS_TSV-}" ]
then
  usage
  echo "ERROR: Either -p or -q must be specified."
  exit 1
fi

if [ "${OUTFILE}" = "-" ]
then
  OUTFILE="/dev/stdout"
fi

mkdir -p "${TMPDIR}"

# Convert proteins to tsv
if [ -z "${PROTEINS_TSV-}" ]
then
  PROTEINS_TSV="${TMPDIR}/proteins.tsv"
  fasta_to_tsv.sh "${PROTEINS}" > "${PROTEINS_TSV}"
fi

GENOME_TSV="${TMPDIR}/genome.tsv"
fasta_to_tsv.sh "${GENOME}" > "${GENOME_TSV}"


# Needed for child processes
export GENOME_TSV
export PROTEINS_TSV
export TMPDIR

# The exit 255 guard is necessary because xargs ignores all other error-codes.
grep -v "^#" "${BED}" \
| cut -f1-4 \
| tr '\n' '\t' \
| xargs -d '\t' -n 4 -P "${NCPU}" \
  bash -eu -c '
    spaln_region.sh \
      -g "${GENOME_TSV}" \
      -p "${PROTEINS_TSV}" \
      -c "$0" \
      -s "$1" \
      -e "$2" \
      -i "$3" \
      -o ${TMPDIR}/tmp$$ \
    || exit 255
  '

sortgrcd -C 60 -E 2 -F 2 -H 4 -S a -O 0 ${TMPDIR}/tmp*.grd > "${OUTFILE}"