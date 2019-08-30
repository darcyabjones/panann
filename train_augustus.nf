#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # panann

    A pipeline to predict genes in a multiple fungal genomes.

    ## Usage

    ```bash
    ```

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.


    ## Output
    aligned_reads - Fastq reads aligned to each genome as bams.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

params.genome = false
params.annotations = false
params.add_utr = false
params.nocrf = false


// Tidy gff
// Tidy utrs

process addUtrs {

    label "bedtools"
    label "small_task"

    when:
    params.add_utr

    input:
    file genome genes.gff3
    file add_utr ests.gff3

    output:
    file genes, best_match_ids

    script:
    """
    add_utr.sh \
      -g "genes.gff3" \
      -e "ests.gff3" \
      -t ./tmp \
      -o "genes_with_utrs.gtf"

    mv tmp/best_match_ids.tsv ./
    """
}


// Extract proteins

process filterProteins {
    label "posix"
    label "small_task"

    input:
    proteins

    output:
    filtered_proteins

    script:
    """
    # This removes proteins with internal stop codons.
    # Convert to tab-separated file
    awk '
      /^>/ {
        b=gensub(/^>\\s*(\\S+).*$/, "\\\\1", "g", \$0);
        printf("%s%s\\t", (N>0?"\\n":""), b);
        N++;
        next;
      }
      {
        printf("%s", \$0)
      }
      END {
        printf("\\n");
      }
    ' < "${PROTEIN_FILE}" \
    | sed 's/\\*\$//g' \
    | awk -F '\\t' '!($2 ~ /\\*/) {printf(">%s\\n", \$0)}' \
    | tr '\\t' '\\n' \
    > "${OUTDIR}/complete.faa"
    """
}


process findSelfMatches {

    label "mmseqs"
    label "big_task"

    input:
    filtered_proteins

    output:

    script:
    """
    mkdir -p "proteins"
    mmseqs createdb "${FILTERED_PROTEINS}" "proteins/db"

    mkdir -p "clustered"
    mkdir -p "tmp"
    mmseqs cluster \
      "proteins/db" \
      "clustered/db" \
      "tmp" \
      --min-seq-id 0.6 \
      -c 0.5 \
      --cov-mode 0

    mmseqs createtsv \
      "proteins/db" \
      "proteins/db" \
      "clustered/db" \
      "clustered.tsv"

    rm -rf -- "tmp"

    awk '$1 != $2' "clustered.tsv" \
    | grep \
      -F \
      -f <(cut -f2 01-get_utrs/best_match_ids.tsv) \
    | sort -u -k1,1 \
    > "clustered_utr.tsv"

    cut -f 2 "clustered_utr.tsv" > "clustered_utr_ids.txt" 

    grep \
      -F \
      -f <(cut -f 1 clustered_utr.tsv) \
      -v \
    < "clustered.tsv" \
    | cut -f 1 \
    | sort -u \
    > "clustered_non_utr_ids.txt"

    cat "clustered_utr_ids.txt" "clustered_non_utr_ids.txt" \
    | grep -v "^#" \
    | sort -u \
    > "non_redundant_ids.txt"
    """
}


process getTrainingSet {

}
