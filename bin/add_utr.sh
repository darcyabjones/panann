#!/usr/bin/env bash

set -euo pipefail

# Default temp directory
TMPDIR="tmp$$"
OUTFILE="/dev/stdout"

# Set some defaults
TIDY=false
GFF3_OUT=false

usage() {
  echo 'Runs

USAGE: add_utr.sh -g GENOME -p PROTEINS -b BED -n NCPUS -t TMPDIR > results.gff

Arguments:
-g -- The gene annotations in GFF3 format [REQUIRED].
-e -- The est matches in GFF3 gene format [REQUIRED].
-i -- Tidy the input and add introns. Required genometools installed on your PATH.
-f -- Write output as GFF3 instead of GTF. Requires genometools installed on your PATH.
-t -- Where to store temporary files [DEFAULT: "tmp$$"].
-o -- Where to write the output gtf/gff to [DEFAULT: stdout].


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
  < "${GFF}" \
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


while getopts ":hg:e:t:o:if" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    g ) GENE="${OPTARG}" ;;
    e ) EST="${OPTARG}" ;;
    i ) TIDY=true ;;
    f ) GFF3_OUT=true ;;
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


# Find transcripts with introns not in the gff.
EST_NOVEL_INTRONS="${TMPDIR}/est_novel_introns.gff3"
awk '$10 == "."' \
< "${EST_SHARED_INTRONS}" \
| cut -f-9 \
> "${EST_NOVEL_INTRONS}"

EST_NOVEL_INTRON_IDS="${TMPDIR}/est_novel_intron_ids.txt"
awk '{b=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $9); print b;}' \
< "${EST_SHARED_INTRONS}" \
| sort -u \
> "${EST_NOVEL_INTRON_IDS}"


# Find genes with introns not in the ests.
GENE_NOVEL_INTRONS="${TMPDIR}/gene_novel_introns.gff3"
awk '$10 == "."' \
< "${GENE_SHARED_INTRONS}" \
| cut -f-9 \
> "${GENE_NOVEL_INTRONS}"

GENE_NOVEL_INTRON_IDS="${TMPDIR}/gene_novel_intron_ids.txt"
awk '{b=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $9); print b;}' \
< "${GENE_NOVEL_INTRONS}" \
| sort -u \
> "${GENE_NOVEL_INTRON_IDS}"


# Find the number of introns per gene in ref gff.
GENE_NUM_INTRONS="${TMPDIR}/gene_num_introns.tsv"
cut -f9 "${GENE_INTRON}" \
| sort \
| uniq -c \
| awk '
  BEGIN {OFS="\t"}
  {
    $2=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $2);
    print $2, $1;
  }
' \
| sort \
> "${GENE_NUM_INTRONS}"

# Find the number of shared introns per gene in ref gff.
SHARED_NUM_INTRONS="${TMPDIR}/shared_num_introns.tsv"
awk '
  BEGIN {OFS="\t"}
  {
    a=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $9);
    b=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $18);
    print a, b;
  }
' < "${EST_SHARED_INTRONS}" \
| sort \
| uniq -c \
| awk 'BEGIN {OFS="\t"} { print $3, $2, $1 }' \
| sort \
> "${SHARED_NUM_INTRONS}"


# Find cases where the number of shared introns is the same as
# the number of introns in the ref gene.
# Output the id of the est.
ESTS_WITH_SAME_INTRONS_IDS="${TMPDIR}/transcript_same_num_introns_ids.tsv"
join \
  -j 1 \
  "${GENE_NUM_INTRONS}" \
  "${SHARED_NUM_INTRONS}" \
| awk '$2 == $4 {print $3}' \
> "${ESTS_WITH_SAME_INTRONS_IDS}"


# Select the ests features with the same introns.
# Excluding those with novel introns.
EST_MRNA_OVERLAP_NO_NOVEL="${TMPDIR}/est_mrna_overlap_nonovel.gff3"
grep \
  -f "${EST_NOVEL_INTRON_IDS}" \
  -F \
  -v \
  "${EST_MRNA_OVERLAP}" \
| grep \
  -f "${ESTS_WITH_SAME_INTRONS_IDS}" \
  -F \
> "${EST_MRNA_OVERLAP_NO_NOVEL}"


# Select ref genes with the same introns.
GENE_MRNA_NO_NOVEL="${TMPDIR}/gene_mrna_no_novel.gff3"
grep \
  -f "${GENE_NOVEL_INTRON_IDS}" \
  -F \
  -v \
  "${GENE_MRNA}" \
> "${GENE_MRNA_NO_NOVEL}"


# Get the overlapping ests and genes with same introns side by side
EST_MRNA_OVERLAP_NO_NOVEL_JOINED="${TMPDIR}/est_mrna_no_novel_joined.gffish"
bedtools intersect \
  -a "${EST_MRNA_OVERLAP_NO_NOVEL}" \
  -b "${GENE_MRNA_NO_NOVEL}" \
  -F 1.0 \
  -s \
  -wa \
  -wb \
> "${EST_MRNA_OVERLAP_NO_NOVEL_JOINED}"


# Select the ids of the joined ests+genes and the score column.
MATCH_IDS="${TMPDIR}/match_ids.tsv"
awk '
  BEGIN {OFS="\t"}
  {
    a=gensub(/.*ID=([^;]+).*/, "\\1", "g", $9);
    b=gensub(/.*ID=([^;]+).*/, "\\1", "g", $18);
    print a, b, $6;
  }
' \
< "${EST_MRNA_OVERLAP_NO_NOVEL_JOINED}" \
| sort -u \
> "${MATCH_IDS}"


# For each gene select the highest scoring matching est.
BEST_MATCH_IDS="${TMPDIR}/best_match_ids.tsv"
sort \
  -u \
  -k2,2 \
  -k3,3gr \
  "${MATCH_IDS}" \
| sort -k2,2 -u \
> "${BEST_MATCH_IDS}"


# Select the features with the best matching ids.
EST_SELECTED="${TMPDIR}/est_selected.gff3"
grep \
  -f <(cut -f1 "${BEST_MATCH_IDS}") \
  -F \
  "${EST}" \
> "${EST_SELECTED}"


# Select just the exons and mRNAs remove everything from attributes but ids.
# Join with the matching ids table to add a new column
# Transfer the gene id to the est attribute column and write out as gtf file.
EST_SELECTED_EXONS="${TMPDIR}/est_shared_exons.gtf"
awk -F'\t' '
  BEGIN {OFS="\t"}
  $3 == "mRNA" {
    $9=gensub(/.*ID=([^;]+).*/, "\\1", "g", $9);
    print
  }
  $3 == "exon" {
    $9=gensub(/.*Parent=([^;]+).*/, "\\1", "g", $9);
    print
  }
' < "${EST_SELECTED}" \
| sort -k9,9 \
| join -1 9 -2 2 - <(sort -k2,2 "${BEST_MATCH_IDS}") \
| awk '
    BEGIN {OFS="\t"}
    $4 == "exon" {
      print $2, $3, $4, $5, $6, $7, $8, $9, "transcript_id \""$10"\"; gene_id \""$10"\";";
    }
  ' \
> "${EST_SELECTED_EXONS}"


# Select just the CDSs from the ref gff as a gtf file.
GENE_CDS="${TMPDIR}/gene_cds.gtf"
awk -F'\t' '
  BEGIN {OFS="\t"}
  $3 == "CDS" {
    $9=gensub(/.*Parent=([^;]+).*/, "transcript_id \"\\1\"; gene_id \"\\1\";", "g", $9);
    print
  }
' "${GENE}" \
> "${GENE_CDS}"


# Combine the Exons and CDSs into a gtf
GENES_WITH_UTRS_GTF="genes_with_utrs.gtf"
cat "${EST_SELECTED_EXONS}" "${GENE_CDS}" \
| sork -k1,1 -k4,4n -k9,9 \
> "${GENES_WITH_UTRS_GTF}"


if ${GFF3_OUT}
then
    # Convert the gtf to gff, this is the final output.
    gt gtf_to_gff3 -tidy "${GENES_WITH_UTRS_GTF}" > "${OUTFILE}"
else
    cat "${GENES_WITH_UTRS_GTF}" > "${OUTFILE}"
fi
