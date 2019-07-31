#!/usr/bin/env bash

set -euo pipefail

# Default temp directory
TMPDIR="tmp$$"
NCPU="$(grep -c '^processor' /proc/cpuinfo)"
OUTFILE="/dev/stdout"

usage() {
  echo 'Runs exonerate in parallel within regions specified by a BED file.

USAGE: exonerate_parallel.sh -g GENOME -p PROTEINS -b BED -n NCPUS -t TMPDIR > results.gff

Arguments:
-g -- The genome in fasta format [REQUIRED].
-p -- The proteins in fasta format [REQUIRED].
-b -- A bed file with regions to search in the first 3 columns
      and a comma separated list of protein ids to search against that
      region in the 4th column [REQUIRED].
-n -- How many concurrent tasks to run [DEFAULT: All detected cpus].
-t -- Where to store temporary files [DEFAULT: "tmp$$"].
-o -- Where to write the output gff to [DEFAULT: stdout].


Requires the following on PATH:
- samtools
- exonerate_region.sh
- exonerate
- standard linux utilities: awk, xargs, sed, tr, grep, bash


Example:
```
exonerate_parallel.sh -g genome.fasta -p proteins.fasta -b merged.bed -n 16 -o results.gff
```
'
}


fasta_to_tsv() {
  awk '
    /^>/ {
      b=gensub(/^>\s*(\S+).*$/, "\\1", "g", $0);
      printf("%s%s\t", (N>0?"\n":""), b);
      N++;
      next;
    }
    {
      printf("%s", $0)
    }
    END {
      printf("\n");
    }
  ' < "$1"
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":hg:p:b:n:t:" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    g ) GENOME="${OPTARG}" ;;
    p ) PROTEINS="${OPTARG}" ;;
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


if   [ -z "${GENOME-}" ] \
  || [ -z "${PROTEINS-}" ] \
  || [ -z "${BED-}" ]
then
  usage
  echo "ERROR: One or more required arguments have not been provided"
  exit 1
fi

if [ "${OUTFILE}" = "-"]
then
  OUTFILE="/dev/stdout"
fi

mkdir -p "${TMPDIR}"

# Needed for child processes
export GENOME
export TMPDIR

# Convert proteins to tsv
fasta_to_tsv "${PROTEINS}" > "${TMPDIR}/proteins.tsv"

# Get .fai to avoid race if it doesn't exist already.
samtools faidx "${GENOME}"

grep -v "^#" "${BED}" \
| cut -f1-4 \
| tr '\n' '\t' \
| xargs -d '\t' -n 4 -P "${NCPU}" \
  bash -eu -c '
    exonerate_region.sh \
      -g "${GENOME}" \
      -p "${TMPDIR}/proteins.tsv" \
      -c "$0" \
      -s "$1" \
      -e "$2" \
      -i "$3" \
      -o ${TMPDIR}/tmp$$
  '

cat ${TMPDIR}/tmp*.gff > "${OUTFILE}"
