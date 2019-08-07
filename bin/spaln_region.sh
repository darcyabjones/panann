#!/usr/bin/env bash

set -euo pipefail

# Default output directory
OUT_PREFIX="tmp$$"

usage() {
  echo 'Runs spaln on a region of a genome with a subset of proteins.

USAGE: genome_region.sh -g GENOME -p PROTEINS -c SEQID -s START -e END -i protein1[,protein2[,...]]

Arguments:
-g -- the genome as a tab-separated file
-p -- the proteins to search as a tab-separated file e.g. "id\tsequence"
-c -- the seqid/chromosome to search in the genome
-s -- the start position within the seqid of the genome to subset 1 indexed
-e -- the end position within the seqid of the genome to subset 1 indexed
-i -- A comma separated list of protein ids to search from the -p file.
      These ids must the first column of -p EXACTLY.
-o -- The prefix to save temp files and the output file as [Optional].
      By default will use "tmp$$" where "$$" is the process id (a uniqueish number).
      Temp files will be saved in a directory with this name, and the final gff
      will be saved with this name and the ".gff" extension.

Example:
```
spaln_region.sh -g genome.tsv -p proteins.tsv -c chr1 -s 1000 -e 50000 -i uniref1,uniref2  -o my_results

# Delete the temp files if you want to.
rm -rf -- my_results

# View the results
head my_results.gff
```

NB because -i is provided as a comma separated list, any real commas in
ids (for whatever reason) should be escaped as "%2C".
This will be converted back to a comma internally.
'
}


select_seqs() {
  TSV=$1
  ID_FILE=$2
  grep -F -f "${ID_FILE}" "${TSV}" | tsv_to_fasta.sh
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi

while getopts ":hg:p:c:s:e:i:o:" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    g ) GENOME="${OPTARG}" ;;
    p ) PROTEINS="${OPTARG}" ;;
    c ) SEQID="${OPTARG}" ;;
    s ) START="${OPTARG}" ;;
    e ) END="${OPTARG}" ;;
    i ) PROTEIN_IDS="${OPTARG}" ;;
    o ) OUT_PREFIX="${OPTARG}" ;;
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
  || [ -z "${SEQID-}" ] \
  || [ -z "${START-}" ] \
  || [ -z "${END-}" ] \
  || [ -z "${PROTEIN_IDS-}" ]
then
  usage
  echo "ERROR: One or more required arguments have not been provided"
  exit 1
fi

TMP_DIR="${OUT_PREFIX}_dir"
mkdir -p "${TMP_DIR}"

# Get protein ids as a file and unescape the commas.
echo "${PROTEIN_IDS}" \
| tr ',' '\n' \
| sed 's/%2C/,/g' \
> "${TMP_DIR}/protein_ids.txt"

# Fetch the protein sequences.
# The grep method (see function) turns out to be much faster than
# using samtools when there are lots of sequences.
select_seqs \
  "${PROTEINS}" \
  "${TMP_DIR}/protein_ids.txt" \
> "${TMP_DIR}/protein_queries.fasta"

# Fetch the chromosome/contig of interest.
select_seqs \
  "${GENOME}" \
  <(echo "${SEQID}") \
> "${TMP_DIR}/genome_target.fasta"

TRANS_TABLE=1

# Allows std and at/an
ALLOW_ALT_SPLICE=1
MIN_INTRON_LEN=20

# -M specifies single best match
# -LS specifies local alignment. Doesn't try to fit the whole protein.
# -yS use "salvage" procedure.

# -H min alignment score for inclusion, was 35.
MIN_SCORE=50


  #-Q0 \

spaln \
  -C${TRANS_TABLE} \
  -ya${ALLOW_ALT_SPLICE} \
  -yL${MIN_INTRON_LEN} \
  -M \
  -LS \
  -yS \
  -H${MIN_SCORE} \
  -KP \
  -O12 \
  -o${OUT_PREFIX} \
  "${TMP_DIR}/genome_target.fasta ${START} ${END}" \
  "${TMP_DIR}/protein_queries.fasta"

  #--softmasktarget

# Filter out the exonerate junk and fix the seqid, start, and end.
#grep -v -E '^#|^-|^Command line|^Hostname' "${OUT_PREFIX}/exonerate.gff" \
#| awk -F '\t' -v seqid="${SEQID}" -v start="${START}" '
#  BEGIN { OFS="\t" }
#  {
#    $1=seqid;
#    $4=($4 + start);
#    $5=($5 + start);
#    print
#  }
#' > "${OUT_PREFIX}.gff"
