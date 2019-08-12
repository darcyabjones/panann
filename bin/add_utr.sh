#!/usr/bin/env bash

set -euo pipefail

# Default temp directory
TMPDIR="tmp$$"
NCPU="$(grep -c '^processor' /proc/cpuinfo)"
OUTFILE="/dev/stdout"

# Set some defaults
TIDY=0
GFF3_OUT=0

usage() {
  echo 'Runs

USAGE: exonerate_parallel.sh -g GENOME -p PROTEINS -b BED -n NCPUS -t TMPDIR > results.gff

Arguments:
-g -- The gene annotations in GFF3 format [REQUIRED].
-e -- The est matches in GFF3 gene format [REQUIRED].
-i -- Tidy the input and add introns. Required genometools installed on your PATH.
-f -- Write output as GFF3 instead of GTF. Requires genometools installed on your PATH.
-n -- How many concurrent tasks to run [DEFAULT: All detected cpus].
-t -- Where to store temporary files [DEFAULT: "tmp$$"].
-o -- Where to write the output gff to [DEFAULT: stdout].


Requires the following on PATH:
- bedtools
- genometools [Optional]
- standard linux utilities: awk, sed, tr, grep, bash


Example:
```
```
'
}


select_gff_type() {
  TYPE="$1"
  FILE="$2"
  awk -v type="${TYPE}" '$3 == type' "${FILE}"
}

get_overlapping_features() {
  TARGET="$1"
  QUERY="$2"

  bedtools intersect \
    -a "${TARGET}" \
    -b "${QUERY}" \
    -wa \
    -s \
    -F 1.0 \
  | sort -k1,1 -k4,4n \
  | uniq \
  | bedtools intersect \
    -a - \
    -b "${QUERY}" \
    -wa \
    -s \
    -F 0.5 \
  | sort -k9 \
  | uniq -u
}

get_gff3_ids() {
  GFF="$1"

  awk '{b=gensub(/.*ID=([^;]+).*/, "\\1", "g", $9); print b;}' \
  < "${OUTDIR}/transcript_overlap.gff3" \
  | sort -u
}

get_overlapping_features_with_ids() {
  TARGET="$1"
  QUERY="$2"
  IDS="$3"

  bedtools intersect \
    -a "${TARGET}" \
    -b "${QUERY}" \
    -s \
  | grep -F -f "${IDS}"
}

get_exact_matching_features() {
  TARGET="$1"
  QUERY="$2"

  bedtools intersect \
    -a "${TARGET}" \
    -b "${QUERY}" \
    -s \
    -f 1.0 \
    -F 1.0 \
    -wa \
    -wb \
    -loj
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":hg:e:n:t:o:if" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    g ) GENE="${OPTARG}" ;;
    e ) EST="${OPTARG}" ;;
    i ) TIDY=1 ;;
    f ) GFF3_OUT=1 ;;
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


if [ -z "${GENE}" ] || [ -z "${EST}" ]
then
  usage
  echo "ERROR: One or more required arguments have not been provided."
  exit 1
fi

if [ ${TIDY} -eq 0 ]
then
  NGENE_INTRONS=$(awk '$3 == "intron"' < "${GENE}" | wc -l)
  NEST_INTRONS=$(awk '$3 == "intron"' < "${EST}" | wc -l)
else
  NGENE_INTRONS=1
  NEST_INTRONS=1
fi

if [ ${NGENE_INTRONS} -eq 0 ] || [ ${NEST_INTRONS} -eq 0 ]
then
  echo "ERROR: one of the input GFFs doesn't have any introns."
  echo "Please add these features or use the '-i' flag."
  exit 2
fi

mkdir -p "${TMPDIR}"


# We need the introns to check that the transcript matches the gene completely.
# Calculate them if user tells us to.
if [ ${TIDY} -eq 1 ]
then
  gt gff3 \
    -tidy \
    -retainids \
    --addintrons \
    -sort <(grep -v "^#" "${GENE}") \
  > "${TMPDIR}/genes.gff3"
  GENE="${TMPDIR}/genes.gff3"

  gt gff3 \
    -tidy \
    -retainids \
    --addintrons \
    -sort <(grep -v "^#" "${EST}") \
  > "${TMPDIR}/ests.gff3"
  EST="${TMPDIR}/ests.gff3"
fi


# Select just the mRNA bits from gffs.
GENE_MRNA="${TMPDIR}/genes_mrna.gff3"
select_gff_type "mRNA" "${GENE}" > "${GENE_MRNA}"

EST_MRNA="${TMPDIR}/est_mrna.gff3"
select_gff_type "mRNA" "${EST}" > "${EST_MRNA}"


# Select just the intron bits from gffs.
GENE_INTRON="${TMPDIR}/genes_intron.gff3"
select_gff_type "intron" "${GENE}" > "${GENE_INTRON}"

EST_INTRON="${TMPDIR}/est_intron.gff3"
select_gff_type "intron" "${EST}" > "${EST_INTRON}"


## STEP 1. find overlapping ests/mRNAs

# Identify which transcript alignments actually hit a known gene.
EST_MRNA_OVERLAP="${TMPDIR}/est_mrna_overlap.gff3"
get_overlapping_features "${EST_MRNA}" "${GENE_MRNA}" > "${EST_MRNA_OVERLAP}"

EST_MRNA_OVERLAP_IDS="${TMPDIR}/est_mrna_overlap_ids.txt"
get_gff3_ids "${EST_MRNA_OVERLAP}" > "${EST_MRNA_OVERLAP_IDS}"


# Identify introns from the overlapping mRNAs that are within target mRNAs.
# This is necessary to exclude UTR introns.
EST_INTRON_OVERLAP="${TMPDIR}/est_intron_overlap.gff3
get_overlapping_features_with_ids \
  "${EST_INTRON}" \
  "${GENE_MRNA}" \
  "${EST_MRNA_OVERLAP_IDS}" \
> "${EST_INTRON_OVERLAP}"


## STEP 2. find introns that are common to ests and mRNAs and ones that are
## novel to one or the other.

# Find all est introns in the overlapping set that are the same as in the gene set.
# Output will have 18 columns, with est introns in first 9 cols,
# and joined gene introns in last 9
EST_SHARED_INTRONS="${TMPDIR}/est_shared_introns.gffish"
get_exact_matching_features \
  "${EST_INTRON_OVERLAP}" \
  "${GENE_INTRON}" \
> "${EST_SHARED_INTRONS}"

# Find all gene introns that are the same as in the overlapping est set.
# Like est shared introns but gene introns are in first 9 columns.
GENE_SHARED_INTRONS="${TMPDIR}/gene_shared_introns.gffish"
get_exact_matching_features \
  "${GENE_INTRON}" \
  "${EST_INTRON_OVERLAP}" \
> "${GENE_SHARED_INTRONS}"



