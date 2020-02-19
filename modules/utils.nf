#!/usr/bin/env nextflow


def symmetric_difference(a, b) {
    return (a + b) - a.intersect(b)
}

process get_univec {

    label "download"
    label "small_task"

    time '3h'

    output:
    path "univec.fasta"

    script:
    """
    wget -O univec.fasta ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
    """
}

process get_augustus_config {

    label "augustus"
    label "small_task"

    time '1h'

    output:
    path "config"

    script:
    """
    cp -r \${AUGUSTUS_CONFIG_PATH} ./config
    """
}


process clean_transcripts {

    label "pasa"
    label "small_task"
    time '2h'

    input:
    path "transcripts.fasta"
    path "univec.fasta"

    output:
    tuple path("transcripts.fasta"),
          path("transcripts.fasta.cln"),
          path("transcripts.fasta.clean")

    script:
    """
    # this user thing is needed for seqclean. Unknown reasons
    export USER="root"
    seqclean "transcripts.fasta" -v "univec.fasta"
    """
}


process tidy_genome {

    label "posix"
    label "small_task"

    time '1h'

    tag "${name}"

    input:
    val min_contig_length
    tuple val(name), path("in.fa")

    output:
    tuple val(name), path("${name}.fasta")

    script:
    """
    # braker panics if the genome has descriptions
    sed -r 's/^(>[^[:space:]]*).*\$/\\1/' in.fa \
    | fasta_to_tsv.sh \
    | awk 'length(\$2) >= ${min_contig_length}' \
    | tsv_to_fasta.sh \
    > "${name}.fasta"
    """
}

process get_faidx {

    label "samtools"
    label "small_task"

    time '1h'

    tag "${name}"

    input:
    tuple val(name), path("orig.fa")

    output:
    tuple val(name), path("${name}.fasta.fai")

    script:
    """
    samtools faidx "orig.fa"

    mv orig.fa.fai "${name}.fasta.fai"
    """
}


process fasta_to_tsv {

    label "posix"
    label "small_task"
    time "1h"

    input:
    tuple val(name), path("seqs.fasta")

    output:
    tuple val(name), path("${name}.tsv")

    script:
    """
    awk '
      /^>/ {
        b=gensub(/^>\\s*(\\S+).*\$/, "\\\\1", "g", \$0);
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
    ' < seqs.fasta \
    > "${name}.tsv"
    """
}


/*
 * Mostly this is just to add intron features.
 */
process tidy_gff3 {

    label "aegean"
    label "small_task"
    time '1h'

    tag "${name}"

    input:
    val analysis
    val source
    tuple val(name), path("in.gff3")

    output:
    tuple val(name), path("${name}_${analysis}_tidied.gff3")

    script:
    """
    grep -v "^#" in.gff3 \
    | gt gff3 \
      -tidy \
      -sort \
      -retainids \
      -addintrons \
      -setsource "${source}" \
    | canon-gff3 -i - \
    > "${name}_${analysis}_tidied.gff3"
    """
}


/*
 * This is a naive merge, it assumes that the contigs were intact.
 */
process combine_and_tidy_gff3 {

    label "aegean"
    label "small_task"
    time '1h'

    tag "${name}"

    input:
    val analysis
    val source
    tuple val(name), path("*chunks.gff")

    output:
    tuple val(name), path("${name}_${analysis}_tidied.gff3")

    script:
    """
    for f in *chunks.gff
    do
      if [ -s "\${f}" ]
      then
        gt gff3 -tidy -sort -addintrons -setsource "${source}" -o "\${f}_tidied.gff3" "\${f}"
      fi
    done

      gt merge -tidy *_tidied.gff3 \
    | canon-gff3 -i - \
    > "${name}_${analysis}_tidied.gff3"
    """
}


process merge_gffs {

    label "genometools"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    val analysis
    tuple val(name), path("to_merge/*gff3")

    output:
    tuple val(name), path("${name}_${analysis}.gff3")

    script:
    """
    gt merge -tidy to_merge/*gff3 > "${name}_${analysis}.gff3"
    """
}


/*
 * Deduplicate identical user provided proteins and concat into single file.
 */
process combine_fastas {

    label "seqrenamer"
    label "small_task"
    time '2h'

    input:
    path "*fasta"

    output:
    path "combined.fasta"
    path "combined.tsv"

    script:
    """
    sr encode \
      --format fasta \
      --column id \
      --deduplicate \
      --upper \
      --drop-desc \
      --strip "*-" \
      --map "combined.tsv" \
      --outfile "combined.fasta" \
      *fasta
    """
}


/*
 * Augustus is pretty slow so we split it into ~16 roughly
 * equally sized chunks to run in parallel.
 * NB this doesnt split within chromosomes/scaffolds/contigs.
 * They are kept intact
 */
process chunkify_genomes {

    label "python3"
    label "small_task"
    time '1h'

    tag "${name}"

    input:
    val nchunks
    tuple val(name),
        path("input.fasta")

    output:
    tuple val(name), path("${name}_chunkfied_*.fasta")

    script:
    """
    chunk_genomes.py -n "${nchunks}" --prefix "${name}_chunkfied_" input.fasta
    """
}


/*
 * Extracts protein and nucleotide sequences from predictions.
 */
process extract_seqs {

    label "genometools"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    input:
    val trans_table
    tuple val(name),
        val(analysis),
        path(gff3),
        path(fasta)

    output:
    tuple val(name),
        val(analysis),
        path("${name}_${analysis}.faa")

    tuple val(name),
        val(analysis),
        path("${name}_${analysis}.fna")

    script:
    """
    gt gff3 \
      -sort \
      -retainids \
      "${gff3}" \
    > sorted.gff3

    gt extractfeat \
      -type CDS \
      -join \
      -translate \
      -retainids \
      -gcode "${trans_table}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "sorted.gff3" \
    > "${name}_${analysis}.faa"

    gt extractfeat \
      -type CDS \
      -join \
      -retainids \
      -gcode "${trans_table}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "sorted.gff3" \
    > "${name}_${analysis}.fna"
    """
}
