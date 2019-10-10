#!/usr/bin/env bash

set -euo pipefail

# Default temp directory
TMPDIR="tmp$$"
NCPU="$(grep -c '^processor' /proc/cpuinfo)"
OUTDIR="out"
MIN_INTRON=5
UTR=false
SPLICE_SITES=""

usage() {
  echo 'Runs augustus in parallel within regions specified by a BED file.

USAGE: augustus_region_parallel.sh -f GENOME -b BED -g HINTS -s SPECIES -n NCPUS -t TMPDIR > results.gff3

Arguments:
-f -- The genome in fasta format [REQUIRED].
-b -- A bed 6 file with regions to search in the first 3 columns
      and the strand (+/-) 6th column [REQUIRED].
-g -- The hints GFF to use [REQUIRED].
-s -- The augustus species to use [REQUIRED].
-c -- The hint weights config file [REQUIRED].
-a -- The path to the augustus config directory [DEFAULT: taken from environment variable AUGUSTUS_CONFIG_PATH]
-p -- A comma separated list of valid splice sites for augustus to predict with evidence.
-u -- Predict genes with UTR models.;
-n -- How many concurrent tasks to run [DEFAULT: All detected cpus].
-t -- Where to store temporary files [DEFAULT: "tmp$$"].
-m -- The minimum intron length in bp (default 20).
-r -- The ncbi genetic code number to use (default 1).
-o -- Directory to write the output gffs to [DEFAULT: out].


Requires the following on PATH:
- samtools
- standard linux utilities: awk, xargs, sed, tr, grep, bash


Example:
```
exonerate_parallel.sh -g genome.fasta -p proteins.fasta -b merged.bed -n 16 -o results.gff
```
'
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":hf:b:g:us:c:a:p:n:t:m:r:o:" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    f ) GENOME="${OPTARG}" ;;
    b ) BED="${OPTARG}" ;;
    g ) HINTS="${OPTARG}" ;;
    u ) UTR=true ;;
    s ) SPECIES="${OPTARG}" ;;
    c ) CONFIG="${OPTARG}" ;;
    a ) AUGUSTUS_CONFIG_PATH="${OPTARG}" ;;
    p ) SPLICE_SITES="${OPTARG}" ;;
    n ) NCPU="${OPTARG}" ;;
    t ) TMPDIR="${OPTARG}" ;;
    m ) MIN_INTRON="${OPTARG}" ;;
    o ) OUTDIR="${OPTARG}" ;;
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

if   [ -z "${GENOME-}" ] || [ -z "${BED-}" ] || [ -z "${HINTS-}" ] || [ -z "${SPECIES-}" ] || [ -z "${CONFIG}" ]
then
  usage
  echo "ERROR: One or more required arguments have not been provided"
  exit 1
fi

mkdir -p "${TMPDIR}"
mkdir -p "${OUTDIR}"

# Because samtools places the fai next to the genome
# it can mess with checkpointing in nextflow.
# Copying avoids that issue.
cp "${GENOME}" "${TMPDIR}/genome.fasta"

SCAFFOLDS=$(cut -f1 "${BED}" | sort -u)
for scaf in ${SCAFFOLDS}
do
  samtools faidx "${GENOME}" "${scaf}" > "${TMPDIR}/${scaf}.fasta"
  getLinesMatching.pl <(echo ${scaf}) 1 < "${HINTS}" > "${TMPDIR}/${scaf}.gff3"
done


# Needed for child processes
export BED
export UTR
export SPECIES
export CONFIG
export AUGUSTUS_CONFIG_PATH
export SPLICE_SITES
export MIN_INTRON
export TMPDIR


run_augustus() {

    FASTA="$1"
    HINTS="$2"
    SCAFFOLD="$3"
    START="$4"
    END="$5"
    STRAND="$6"
    MIN_INTRON="$7"
    SPLICE_SITES="$8"
    UTR="$9"
    SPECIES="${10}"
    AUGUSTUS_CONFIG_PATH="${11}"
    CONFIG="${12}"
    TMPDIR="${13}"

    STRAND_FLAG=$( [ "${STRAND}" == "-" ] && echo "--strand=backward" || echo "--strand=forward" )
    UTR_FLAG=$( [ "${UTR}" == "true" ] && echo "--UTR" || echo "" )

    OUTPREFIX="${TMPDIR}/${SCAFFOLD}_${START}_${END}_${STRAND}_preds"

    augustus \
      --species="${SPECIES}" \
      --AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}" \
      --genemodel=complete \
      --extrinsicCfgFile="${CONFIG}" \
      --hintsfile="${HINTS}" \
      ${STRAND_FLAG} \
      ${UTR_FLAG} \
      --allow_hinted_splicesites="${SPLICE_SITES}" \
      --softmasking=on \
      --alternatives-from-evidence=true \
      --min_intron_len="${MIN_INTRON}" \
      --predictionStart="${START}" \
      --predictionEnd="${END}" \
      --noInFrameStop=true \
      --start=off \
      --stop=off \
      --introns=off \
      --gff3=on \
      --outfile="${OUTPREFIX}.gff3" \
      --errfile="${OUTPREFIX}.err" \
      "${FASTA}"
}

export -f run_augustus

# The exit 255 guard is necessary because xargs ignores all other error-codes.
grep -v "^#" "${BED}" \
| tr '\n' '\t' \
| xargs -d '\t' -n 4 -P "${NCPU}" \
  bash -eu -c '
    run_augustus \
      "${TMPDIR}/${0}.fasta" \
      "${TMPDIR}/${0}.gff3" \
      "${0}" \
      "${1}" \
      "${2}" \
      "${3}" \
      "${MIN_INTRON}" \
      "${SPLICE_SITES}" \
      "${UTR}" \
      "${SPECIES}" \
      "${AUGUSTUS_CONFIG_PATH}" \
      "${CONFIG}" \
      "${TMPDIR}" \
    || (echo "Failed $0, $1, $2, $3"; exit 255)
  '

mkdir -p "${OUTDIR}"
mv ${TMPDIR}/*_preds.gff3 "${OUTDIR}"
mv ${TMPDIR}/*_preds.err "${OUTDIR}"
rm -rf -- "${TMPDIR}"
