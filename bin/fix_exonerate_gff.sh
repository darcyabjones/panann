#!/usr/bin/env bash

awk -F '\t' '
  BEGIN {
    OFS="\t"
    gnum = 0
    gname = "null"
    query = "null"
  }
  /# --- START OF GFF DUMP / {
    gnum++
    gname = "null"
    query = "null"

  }
  !/^#/ && $3 == "gene" {
    gname="gene" gnum
    query=gensub(/.*sequence ([^[:space:];]+).*/, "\\1", "g", $9)
    ori=gensub(/.*gene_orientation ([^[:space:];]+).*/, "\\1", "g", $9)
    identity=gensub(/.*identity ([^[:space:];]+).*/, "\\1", "g", $9)

    print
  }
  !/^#/ {
    $9="gene_id " gid " ; " $9
    print
  }
' < /dev/stdin
