#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # pante


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
params.reads = false
params.repbase = false


/*
 * Sanitise the input
 */

if ( params.reads ) {
    reads = Channel.fromPath(params.reads, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.genome != null && it.read1_file != null && it.read2_file != null }
        .map {[
            it.genome,
            file(it.read1_file, checkIfExists: true),
            file(it.read2_file, checkIfExists: true)
        ]}
}


if ( params.genomes ) {
    genomes = Channel
        .fromPath(params.assemblies, checkIfExists: true, type: "file")
        .map { file -> [file.name, file] }
} else if ( params.table ) {
    genomes = Channel.fromPath(params.reads, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.genome != null }
        .unique { it.genome }
        .map { file(it.genome, checkIfExists: true) }
        .map { file -> [file.name, file] }
} else {
    log.info "Hey I need some genomes to assess please."
    exit 1
}





/*
 * De-novo repeat assemblers
 */

// REPdenovo

/*
 * De-novo repeat finders
 */


/*
 * RED
process runRed {
    label "red"

    input:
    set val(name), file(genome) from ..

    output:
    ...

    """
    """
}
*/

// Repeat modeller
// CARP ? Might need to pipeline this myself, currently exists as separate tools and description of workflow.
// MGEScan-non-ltr

/*
 * Homology based repeat finders
 */

// Repeat masker

/*
 * LTR finders
 */


// LTR_retriever
// LTR_harvest + digest // See ltr_retriever paper for "good" parameters
// LTR_detector

/*
 * SINE finders
 */

// SINE_Scan


/*
 * MITE finders
 */

// MITEfinderII
// MITEHunter

/*
 * Helitron finders
 */

// HELsearch


/*
 * Combine and remove redundancy.
 * Meshclust looks like good tool for clustering, (alternative to cdhit-cdna)
 * Possibly need to look at overlapping genomic annotations?
 * Filtering out false positives? negative Database matches?
 */


/*
 * Multiple sequence alignment and annotation
 * Mafft is always solid aligner choice.
 * Paste probably good option for classification.
 * Classification might go hand in hand with filtering false positives?
 * Maybe fp is first step is to exclude members of clusters,
 *  whereas here it is to exclude entire clusters?
 */
