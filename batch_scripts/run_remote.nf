#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include get_faidx from '../modules/utils'
include fasta_to_tsv from '../modules/utils'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include align_remote_proteins from '../modules/workflows'

def helpMessage() {
    log.info"""
    # panann


    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.min_intron_hard = 5
params.max_intron_hard = 15000
params.max_gene_hard = 20000
params.trans_table = 1

params.table = false
params.genomes = false


workflow {

    main:

    if ( params.remote_proteins ) {
        remote_proteins = Channel
            .fromPath(params.remote_proteins, checkIfExists: true, type: "file")
            .collectFile(newLine: true)
            .map { ["remote", it] }
    } else {
        log.error "Please provide some remote proteins to use with --remote_proteins"
        exit 1
    }

    if ( params.genomes ) {
        genomes = Channel
            .fromPath(params.genomes, checkIfExists: true, type: "file")
            .map { g -> [g.baseName, g] }

    } else {
        genomes = Channel.empty()
    }

    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(genomes, table)
        in_genomes = input_channels.genomes
    } else if ( !params.genomes ){
        log.error "Please provide some genomes to align to."
        exit 1
    } else {
        in_genomes = genomes
    }

    genomes_faidx = get_faidx(in_genomes)
    mmseqs_genome_indices = get_mmseqs_genome_db(in_genomes)

    remote = align_remote_proteins(
        params.min_intron_hard,
        params.max_intron_hard,
        params.max_gene_hard,
        params.trans_table,
        remote_proteins,
        in_genomes,
        genomes_faidx,
        mmseqs_genome_indices
    )

    publish:
    remote.remote_protein_matches to: "${params.outdir}/alignments"
    remote.clustered_remote_protein_matches to: "${params.outdir}/alignments"
    remote.exonerate_matches to: "${params.outdir}/alignments"
    remote.exonerate_augustus_hints to: "${params.outdir}/hints"
    remote.exonerate_evm_hints to: "${params.outdir}/hints"
}
