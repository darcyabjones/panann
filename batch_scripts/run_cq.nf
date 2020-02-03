#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include run_codingquarry from '../modules/workflows'
include get_file from '../modules/workflows'
include handle_table from '../modules/workflows'

params.genomes = false
params.table = false
params.signalp = false

def is_null = { f -> (f == null || f == '') }


workflow {

    main:

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

    } else {
        log.error "Running codingquarry requires stringtie input."
        exit 1
    }

    cq = run_codingquarry(
        params.signalp,
        input_channels.genomes,
        input_channels.stringtie
    )

    publish:
    cq.cq_gff3          to: "${params.outdir}/annotations"
    cq.cq_faa           to: "${params.outdir}/annotations"
    cq.cq_fna           to: "${params.outdir}/annotations"
    cq.cq_dubious       to: "${params.outdir}/annotations"
    cq.cq_fusions       to: "${params.outdir}/annotations"
    cq.cq_overlap       to: "${params.outdir}/annotations"
    cq.cqpm_gff3        to: "${params.outdir}/annotations"
    cq.cqpm_fusions     to: "${params.outdir}/annotations"
    cq.cqpm_overlap     to: "${params.outdir}/annotations"
    cq.cq_gff3_tidied   to: "${params.outdir}/annotations"
    cq.cqpm_gff3_tidied to: "${params.outdir}/annotations"
}
