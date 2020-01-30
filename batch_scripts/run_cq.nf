#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include run_codingquarry from './modules/workflows'

params.genomes = false
params.stringtie = false
params.signalp = false

def is_null = { f -> (f == null || f == '') }


workflow {

    main:

    if ( params.genomes ) {
        genomes = Channel
            .fromPath(params.genomes, checkIfExists: true, type: "file")
            .map { g -> [g.baseName, g] }

    } else {
        log.error "Please provide some genomes to predict genes for with `--genomes`."
        exit 1

    }

    if ( params.stringtie ) {
        stringtie = Channel
            .fromPath(params.stringtie, checkIfExists: true, type: 'file')
            .splitCsv(by: 1, sep: '\t', header: true)
            .filter { (!is_null(it.name) && !is_null(it.stringtie)) }
            .map {[it.name, file(it.stringtie, checkIfExists: true)]}
            .unique()

    } else {
        log.error "Running codingquarry requires stringtie input."
        exit 1
    }

    cq = run_codingquarry(params.signalp, genomes, stringtie)

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
