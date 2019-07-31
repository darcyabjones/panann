#!/usr/bin/env bash

set -euo pipefail

# Default output directory
OUT_PREFIX="tmp$$"

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
  | awk -F '\t' '{ printf(">%s\n%s", $1, $2) }'
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

# Run exonerate
exonerate \
  --query "${OUT_PREFIX}/protein_queries.fasta" \
  --target "${OUT_PREFIX}/genome_target.fasta" \
  --querytype protein \
  --targettype dna \
  --model protein2genome \
  --refine region \
  --percent 70 \
  --score 100 \
  --geneseed 250 \
  --bestn 2 \
  --minintron 20 \
  --maxintron 50000 \
  --geneticcode 1 \
  --showtargetgff yes \
  --showalignment no \
  --showvulgar no \
> "${OUT_PREFIX}/exonerate.gff"

  #--softmasktarget

# Filter out the exonerate junk and fix the seqid, start, and end.
grep -v -E '^#|^-|^Command line|^Hostname' "${OUT_PREFIX}/exonerate.gff" \
| awk -F '\t' -v seqid="${SEQID}" -v start="${START}" '
  BEGIN { OFS="\t" }
  {
    $1=seqid;
    $4=($4 + start);
    $5=($5 + start);
    print
  }
' > "${OUT_PREFIX}.gff"
