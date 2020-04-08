#!/usr/bin/env bash

set -euo pipefail

INFILE="/dev/stdin"
OUTFILE="/dev/stdout"
FIELD="ID"
TYPE="CDS"
ISGFF2=false
SOURCE=""

usage() {
  echo 'USAGE: gff2bed.sh -o out.bed -f Parent -t CDS in.gff3 > results.gff

Arguments:
-o -- Where to write the output to [DEFAULT: stdout].
-f -- The field to take from the attributes column [DEFAULT: "ID"]
-t -- The type of features to select. Should match column 3 in gff.
-2 -- Parse the attribute column as if it is from gff2 instead of gff3.
-s -- Replace the source with this string [DEFAULT: use column 2 in gff].

Requires standard linux utilities: GNU awk and bash

The output will have 7 columns in a bed-like format.
seqid start end {field} score strand {source}

Where field and source can be altered by the -f and -s parameters.
The only non-standard bit is the source column at the end.
This is added because the script is intended to be used alongside the
get_hint_coverage.sh script which needs that field to annotate genes.
'
}


if [ $# -eq 0 ]
then
  usage
  echo "No arguments provided"
  exit 0
fi


while getopts ":ho:f:t:s:2" opt
do
  case "${opt}" in
    h )
      usage
      exit 0
      ;;
    o ) OUTFILE="${OPTARG}" ;;
    f ) FIELD="${OPTARG}" ;;
    t ) TYPE="${OPTARG}" ;;
    s ) SOURCE="${OPTARG}" ;;
    2 ) ISGFF2=true ;;
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

if [ $# -eq 0 ]
then
  INFILES="/dev/stdin"
else
  INFILES=${@}
fi


get_gff3() {
  TYPE=$1
  FIELD=$2
  SOURCE=$3
  FILE=$4

  awk -F '\t' -v field="${FIELD}" -v type="${TYPE}" -v source="${SOURCE}" '
    BEGIN {
      OFS = "\t"
      field_regex =  ".*"field"=([^;[:space:]]+).*"
    }
    $3 == type {
      if (match($9, field_regex)) {
        id=gensub(field_regex, "\\1", "g", $9);
      } else {
        print "ERROR: Encountered a line without the chosen field." >> "/dev/stderr"
        print "The line was: "$0 >> "/dev/stderr"
        exit 2
      }

      if (source == "") {
        line_source = $2
      } else {
        line_source = source
      }

      if ($4 < 1) {
        start = 0
      } else {
        start = $4 - 1
      }
      print $1, start, $5, id, $6, $7, line_source
    }
  ' < "${FILE}"
}

get_gff2() {
  TYPE=$1
  FIELD=$2
  SOURCE=$3
  FILE=$4

  awk -F '\t' -v field="${FIELD}" -v type="${TYPE}" -v source="${SOURCE}" '
    BEGIN {
      OFS = "\t"
      field_regex =  ".*"field"[[:space:]]+([^;[:space:]]+).*"
    }
    $3 == type {
      if (match($9, field_regex)) {
      	id=gensub(field_regex, "\\1", "g", $9);
      } else {
        print "ERROR: Encountered a line without the chosen field." >> "/dev/stderr"
        print "The line was: "$0 >> "/dev/stderr"
        exit 2
      }

      if (source == "") {
        line_source = $2
      } else {
        line_source = source
      }

      if ($4 < 1) {
        start = 0
      } else {
        start = $4 - 1
      }

      print $1, start, $5, id, $6, $7, line_source
    }
  ' < "${FILE}"
}

for f in ${INFILES}
do
  if ${ISGFF2}
  then
    get_gff2 "${TYPE}" "${FIELD}" "${SOURCE}" "${f}" >> "${OUTFILE}"
  else
    get_gff3 "${TYPE}" "${FIELD}" "${SOURCE}" "${f}" >> "${OUTFILE}"
  fi
done
