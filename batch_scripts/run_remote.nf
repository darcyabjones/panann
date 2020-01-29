#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

include tidy_genome from './modules/utils'
include get_faidx from './modules/utils'
include fasta_to_tsv from './modules/utils'

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


// Target genomes to annotate as softmasked fasta files.
// Note the file basename (filename up to but excluding the
// last extension) is used to match genomes with other hints.
params.min_contig_length = 500
params.min_intron_soft = 20
params.min_intron_hard = 5
params.max_intron_hard = 15000
params.max_gene_hard = 20000
params.genomes = false

// Remote proteins aligned to each genome with exonerate.
// As a tsv file, with name and remote columns
params.remote_alignments = false

params.trans_table = 1


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
        log.error "Please provide some genomes to predict genes for with `--genomes`."
        exit 1
    }

    mmseqs_remote_index = get_mmseqs_protein_db(remote_proteins)
    remote_proteins_tsv = fasta_to_tsv(
        params.min_contig_length,
        remote_proteins
    )

    tidied_genomes = tidy_genome(params.min_contig_length, genomes)
    genomes_faidx = get_faidx(tidied_genomes)

    mmseqs_genome_indices = get_mmseqs_genome_db(tidied_genomes)


    publish:
    remote_protein_matches
    clustered_remote_protein_matches
    exonerate_matches
}
