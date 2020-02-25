#!/usr/bin/env bash

set -euo pipefail

OUTFILE="/dev/stdout"
TMPDIR="./get_hint_coverage$$"
TYPE="CDS"


usage() {
  echo 'USAGE: get_hint_coverage.sh -o out.gff3 ref.gff3 hints*.bed

Arguments:
-o -- Where to write the output to [DEFAULT: stdout].
-t -- The type of features to annotate with hints. Should match column 3 in gff of REF.
-m -- Where to store temporary files.

Requires standard linux utilities GNU awk and bash, and bedtools

The hint beds should have the structure (tab separated):
seqid start end name score strand source

The output GFF will be annotated with hint matches in the attributes column.
The hint= attribute has the structure 'source name start end score num_intersect'.
Where source, name, start, end and score are columns 7, 4, 2, 3, and 5.
The seqid and strand will always be the same as the attached feature.
num_intersect is the number of bases that the hint matches the feature by.
A perfect match will have the same start and end coordinates.
Coordinates are 1-based inclusive (like gff3).
'
}

if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":ho:t:m:" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    o ) OUTFILE="${OPTARG}" ;;
    t ) TYPE="${OPTARG}" ;;
    m ) TMPDIR="${OPTARG}" ;;
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

shift $((OPTIND - 1))

if [ $# -lt 2 ]
then
  usage
  echo "ERROR: Require at least 2 positional arguments." 1>&2
  exit 1
else
  REF="${1}"
  HINTS="${@:2}"
fi

# Bedtools adds an extra column more than 1 -b file provided to intersect.
# We need to offset some of downstream analyses because of this.
if [ $# -gt 2 ]
then
  EXTRA_COL=1
else
  EXTRA_COL=0
fi


mkdir -p "${TMPDIR}"

awk -F '\t' '
  BEGIN {OFS="\t"}
  !/^#/ {
    if ($4 > $5) {
      tmp=$5
      $5=$4
      $4=tmp
    }
    $4-=1
  }
  {print}
' "${REF}" \
> "${TMPDIR}/ref_zero_based.gff3"


bedtools intersect \
  -loj \
  -s \
  -a "${TMPDIR}/ref_zero_based.gff3" \
  -b ${HINTS} \
> "${TMPDIR}/ref_with_inter.tsv"


awk -F '\t' -v extra_col="${EXTRA_COL}" '
  BEGIN {OFS="\t"}
  $3 == "CDS" {
    gstart=$4
    gend=$5

    if (extra_col=="0") {
      hseqid=$10
      hstart=$11
      hend=$12
      hname=$13
      hscore=$14
      hstrand=$15
      hsource=$16
    } else {
      hseqid=$11
      hstart=$12
      hend=$13
      hname=$14
      hscore=$15
      hstrand=$16
      hsource=$17
    }

    if (hstart >= gend || hend <= gstart || hstart == hend || gstart == gend) {
      ninter = 0
    } else {
      l=(gstart > hstart ? gstart : hstart)
      r=(gend < hend ? gend : hend)
      ninter = r - l
    }

    if (hsource == ".") {
      hint="."
    } else {
      hint=hsource" "hname" "hstart + 1" "hend" "hscore" "ninter
    }
    print $1, $2, $3, $4 + 1, $5, $6, $7, $8, $9, hint
  }
  $3 != "CDS" { print $1, $2, $3, $4 + 1, $5, $6, $7, $8, $9, "." }
' < "${TMPDIR}/ref_with_inter.tsv" \
> "${TMPDIR}/ref_with_inter_hints.tsv"

sort -k1,1 -k2,2 -k3,3 -k4,4 -k 5,5 -k6,6 -k7,7 -k8,8 -k9,9 -k10,10 "${TMPDIR}/ref_with_inter_hints.tsv" \
| uniq > "${TMPDIR}/ref_with_inter_sorted.tsv"

bedtools groupby -i "${TMPDIR}/ref_with_inter_sorted.tsv" -g 1-9 -c 10 -o distinct > "${TMPDIR}/ref_grouped.tsv"

awk -F'\t' '
  BEGIN {OFS="\t"}
  $10 == "." { print $1, $2, $3, $4, $5, $6, $7, $8, $9 }
  $10 != "." { print $1, $2, $3, $4, $5, $6, $7, $8, $9";hint="$10 }
' "${TMPDIR}/ref_grouped.tsv" \
| sort -k1,1 -k4,4n -k5,5n -k7,7 -k3,3 -k2,2 -k9,9 \
| uniq \
> "${OUTFILE}"

rm -rf -- "${TMPDIR}"
