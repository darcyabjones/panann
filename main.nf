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
    cat <<EOF > known_sites.tsv
    genome	gff	confidence
    genome1.fasta	genome1.gff3	8
    EOF

    nextflow run main.nf \
      --genomes "genomes/*.fasta" \
      --known_sites known_sites.tsv \
      --transcripts transcripts.fasta \
      --proteins proteins.fasta \
      --remote_proteins fungal_proteins.fasta
    ```

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false

params.known_sites = false

params.transcripts = false

params.proteins = false
params.remote_proteins = false

params.genome_alignment = false

params.busco_lineage = false

/*
 * Sanitise the input
 */

if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.baseName, g] }
        .set { genomes }
} else {
    log.error "Please provide some genomes to predict genes for with `--genomes`."
    exit 1
}


if ( params.transcripts ) {
    Channel
        .fromPath(params.transcripts, checkIfExists: true, type: "file")
        .collectFile(name: "transcripts.fasta", newLine: true, sort: "deep")
        .first()
        .set { transcripts }
} else {
    transcripts = false
}

if ( params.proteins ) {
    Channel
        .fromPath(params.proteins, checkIfExists: true, type: "file")
        .collectFile(name: "proteins.fasta", newLine: true, sort: "deep")
        .first()
        .set { proteins }
} else {
    proteins = false
}

if ( params.remote_proteins ) {
    Channel
        .fromPath(params.remote_proteins, checkIfExists: true, type: "file")
        .collectFile(name: "remote_proteins.fasta", newLine: true, sort: "deep")
        .first()
        .set { remoteProteins }
}

if ( params.busco_lineage ) {
    Channel
        .fromPath(params.busco_lineage, checkIfExists: true, type: "file")
        .first()
        .set { buscoLineage }
}


genomes.into {
    genomes4Busco;
    genomes4SpalnIndex;
}


if ( params.busco_lineage ) {

    process runBusco {
        label "busco"
        label "medium_task"
        tag { name }

        publishDir "${params.outdir}/busco"

        input:
        set val(name), file(fasta) from genomes4Busco
        file "lineage" from buscoLineage

        output:
        file "${name}" into buscoResults

        script:
        """
        run_BUSCO.py \
          --in "${fasta}" \
          --out "${name}" \
          --cpu ${task.cpus} \
          --mode "genome" \
          --lineage_path "lineage"

        mv "run_${name}" "${name}"
        """
    }
}


// If genome alignment not provided, run sibelliaz
// Should there be a step to run transcript assembly with trinity?

process getSpalnIndex {

    label "spaln"
    label "small_task"

    tag { name }

    input:
    set val(name), file(genome) from genomes4SpalnIndex

    output:
    set val(name),
        file("${genome.baseName}.bkn"),
        file("${genome.baseName}.ent"),
        file("${genome.baseName}.idx"),
        file("${genome.baseName}.bkp"),
        file("${genome.baseName}.grp"),
        file("${genome.baseName}.seq") into spalnIndices

    """
    makeidx.pl -inp ${genome}
    """
}


// 1 align transcripts to all genomes

process AlignSpalnTranscripts {
    label "spaln"
    label "medium_task"

    cpus 4

    publishDir "aligned"

    tag { name }

    when:
    params.transcripts

    input:
    file transcripts

    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices

    output:
    set val(name), file("${name}.gff3")

    """
    spaln \
      -L \
      -M3 \
      -O0 \
      -Q7 \
      -S3 \
      -TDothideo \
      -yX \
      -yS \
      -ya3 \
      -t ${task.cpus} \
      -d "${name}" \
      "${transcripts}" \
    > "${name}.gff3"
    """
}

// Align with minimap2


// 2 align proteins to all genomes

// 3 align remote proteins to all genomes

// 4 run busco on all genomes

// 5 Run genemark, pasa, braker2, codingquarry on all genomes using previous steps.

// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats
