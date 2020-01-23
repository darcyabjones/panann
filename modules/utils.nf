#!/usr/bin/env nextflow

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
}


process get_faidx {

    label "samtools"
    label "small_task"

    time '1h'

    tag "${name}"

    input:
    val min_contig_length
    tuple val(name), path("orig.fa")

    output:
    tuple val(name), path("${name}.fasta"), path("${name}.fasta.fai")

    script:
    """
    # braker panics if the genome has descriptions
    sed -r 's/^(>[^[:space:]]*).*\$/\\1/' orig.fa \
    | fasta_to_tsv.sh \
    | awk 'length(\$2) >= ${min_contig_length}' \
    | tsv_to_fasta.sh \
    > "${name}.fasta"


    samtools faidx "${name}.fasta"
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

    tag "${name} - ${paramset}"

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
      if [ -s "\${f}.tmp" ]
      then
        gt gff3 -tidy -sort -addintrons -setsource "${source}" -o "\${f}_tidied.gff3" "\${f}.tmp"
      fi
    done

      gt merge -tidy *_tidied.gff3 \
    | canon-gff3 -i - \
    > "${name}_${analysis}_tidied.gff3"
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
    tuple path("combined.fasta"), path("combined.tsv")

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
    set val(name),
        val(analysis),
        path("${name}_${analysis}.faa")

    tuple val(name),
        val(analysis),
        path("${name}_${analysis}.fna")

    script:
    """
    gt extractfeat \
      -type CDS \
      -join \
      -translate \
      -retainids \
      -gcode "${trans_table}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "${gff3}" \
    > "${name}_${analysis}.faa"

    gt extractfeat \
      -type CDS \
      -join \
      -retainids \
      -gcode "${trans_table}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      "${gff3}" \
    > "${name}_${analysis}.fna"
    """
}
