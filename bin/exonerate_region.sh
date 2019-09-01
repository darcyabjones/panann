#!/usr/bin/env bash

set -Euo pipefail

# Default output directory
OUT_PREFIX="tmp$$"
MIN_INTRON=20
MAX_INTRON=200000
GENCODE=1

usage() {
  echo 'Runs exonerate on a region of a genome with a subset of proteins.
Automatically fixes start and stop positions.

USAGE: exonerate_region.sh -g GENOME -p PROTEINS -c SEQID -s START -e END -i protein1[,protein2[,...]]

Arguments:
-g -- the genome as a fasta file
-p -- the proteins to search as a tab-separated file e.g. "id\tsequence"
-c -- the seqid/chromosome to search in the genome
-s -- the start position within the seqid of the genome to subset 1 indexed
-e -- the end position within the seqid of the genome to subset 1 indexed
-i -- A comma separated list of protein ids to search from the -p file.
      These ids must the first column of -p EXACTLY.
-m -- The minimum intron length in bp (default 20).
-x -- The maximum intron length in bp (default 200000).
-r -- The ncbi genetic code number to use (default 1).
-o -- The prefix to save temp files and the output file as [Optional].
      By default will use "tmp$$" where "$$" is the process id (a uniqueish number).
      Temp files will be saved in a directory with this name, and the final gff
      will be saved with this name and the ".gff" extension.

Example:
```
exonerate_region.sh -g genome.fasta -p proteins.tsv -c chr1 -s 1000 -e 50000 -i uniref1,uniref2  -o my_results

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


get_proteins() {
  TSV=$1
  ID_FILE=$2
  grep -F -f "${ID_FILE}" "${TSV}" \
  | awk -F '\t' '{ printf(">%s\n%s\n", $1, $2) }'
}


run_exonerate() {
  PROTEINS="$2"
  GENOME="$1"
  MIN_INTRON="$3"
  MAX_INTRON="$4"
  GENCODE="$5"

  exonerate \
    --query "${PROTEINS}" \
    --target "${GENOME}" \
    --querytype protein \
    --targettype dna \
    --model protein2genome \
    --refine region \
    --percent 70 \
    --score 100 \
    --geneseed 250 \
    --bestn 2 \
    --minintron "${MIN_INTRON}" \
    --maxintron "${MAX_INTRON}" \
    --geneticcode "${GENCODE}" \
    --showtargetgff yes \
    --showalignment no \
    --showvulgar no
}

run_exonerate_loop() {
  PROTEINS="$2"
  GENOME="$1"
  MIN_INTRON="$3"
  MAX_INTRON="$4"
  GENCODE="$5"

  NSEQS=$(wc -l < "${PROTEINS}")

  set +e
  for i in $(seq 1 2 ${NSEQS})
  do
    tail -n+${i} "${PROTEINS}" | head -n2 > "${PROTEINS}_${i}.fasta"
    run_exonerate \
      "${PROTEINS}_${i}.fasta" \
      "${GENOME}" \
      "${MIN_INTRON}" \
      "${MAX_INTRON}" \
      "${GENCODE}" || true
  done

  set -e  
}

filter_exonerate() {
  GFF="$1"
  SEQID="$2"
  START="$3"

  # Filter out the exonerate junk and fix the seqid, start, and end.
  gawk -F '\t' -v seqid="${SEQID}" -v start="${START}" '
    BEGIN { OFS="\t" }
    !/^#|^-|^Command line|^Hostname/ {
      $1=seqid;
      $4=($4 + start);
      $5=($5 + start);
      $9=gensub(/sequence\s+([^\s;]+)\s+;/, "sequence \\1_" seqid "_" start " ;", "g", $9 );
      $9=gensub(/gene_id\s+([^\s;]+)\s+;/, "gene_id \\1_" seqid "_" start " ;", "g", $9 );
      print;
      next;
    }
    { print }
  ' < "${GFF}"
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi

while getopts ":hg:p:c:s:e:i:m:x:r:o:" opt
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
    m ) MIN_INTRON="${OPTARG}" ;;
    x ) MAX_INTRON="${OPTARG}" ;;
    r ) GENCODE="${OPTARG}" ;;
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

mkdir -p "${OUT_PREFIX}"

# Get protein ids as a file and unescape the commas.
echo "${PROTEIN_IDS}" \
| tr ',' '\n' \
| sed 's/%2C/,/g' \
> "${OUT_PREFIX}/protein_ids.txt"

# Fetch the genomic region.
samtools faidx \
  "${GENOME}" \
  "${SEQID}:${START}-${END}" \
> "${OUT_PREFIX}/genome_target.fasta"

# Fetch the protein sequences.
# The grep method (see function) turns out to be much faster than
# using samtools when there are lots of sequences.
get_proteins \
  "${PROTEINS}" \
  "${OUT_PREFIX}/protein_ids.txt" \
> "${OUT_PREFIX}/protein_queries.fasta"


firsttry() {
  OUT_PREFIX="$1"
  MIN_INTRON="$2"
  MAX_INTRON="$3"
  GENCODE="$4"

  run_exonerate \
    "${OUT_PREFIX}/genome_target.fasta" \
    "${OUT_PREFIX}/protein_queries.fasta" \
    "${MIN_INTRON}" \
    "${MAX_INTRON}" \
    "${GENCODE}" \
  > "${OUT_PREFIX}/exonerate.gff"

  filter_exonerate "${OUT_PREFIX}/exonerate.gff" "${SEQID}" "${START}" > "${OUT_PREFIX}.gff"
}

secondtry() {
  OUT_PREFIX="$1"
  MIN_INTRON="$2"
  MAX_INTRON="$3"
  GENCODE="$4"

  run_exonerate_loop \
    "${OUT_PREFIX}/genome_target.fasta" \
    "${OUT_PREFIX}/protein_queries.fasta" \
    "${MIN_INTRON}" \
    "${MAX_INTRON}" \
    "${GENCODE}" \
  > "${OUT_PREFIX}/exonerate.gff"

  filter_exonerate "${OUT_PREFIX}/exonerate.gff" "${SEQID}" "${START}" > "${OUT_PREFIX}.gff"
}



# Run exonerate
firsttry "${OUT_PREFIX}" "${MIN_INTRON}" "${MAX_INTRON}" "${GENCODE}"

ECODE=$?
if [ ${ECODE} -eq 139 ]
then
  secondtry "${OUT_PREFIX}" "${MIN_INTRON}" "${MAX_INTRON}" "${GENCODE}"
else
  exit "${ECODE}"
fi
