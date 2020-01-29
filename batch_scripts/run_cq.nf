#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include codingquarry from './modules/predictors'
include codingquarrypm from './modules/predictors'
include signalp from './modules/predictors'
include deepsig from './modules/predictors'

include tidy_gff3 as tidy_codingquarry_gff3 from './modules/utils'
include tidy_gff3 as tidy_codingquarrypm_gff3 from './modules/utils'

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

    (cq_fixed_gff3, cq_gff3, cq_faa, cq_fna, cq_dubious, cq_fusions, cq_overlap) = codingquarry(
        stringtie.join(genomes)
    )

    if (params.signalp) {
        cq_secreted = signalp(cq_faa)
    } else {
        cq_secreted = deepsig(cq_faa)
    }

    (cqpm_gff3, cqpm_fusions, cqpm_overlap) = codingquarrypm(
        stringtie.join(genomes).join(cq_gff3).join(cq_secreted)
    )


    cq_gff3_tidied = tidy_codingquarry_gff3(
        "codingquarry",
        "codingquarry",
        cq_fixed_gff3
    )

    cqpm_gff3_tidied = tidy_codingquarrypm_gff3(
        "codingquarrypm",
        "codingquarrypm",
        cqpm_gff3
    )

    publish:
    cq_gff3          to: "${params.outdir}/annotations"
    cq_faa           to: "${params.outdir}/annotations"
    cq_fna           to: "${params.outdir}/annotations"
    cq_dubious       to: "${params.outdir}/annotations"
    cq_fusions       to: "${params.outdir}/annotations"
    cq_overlap       to: "${params.outdir}/annotations"
    cqpm_gff3        to: "${params.outdir}/annotations"
    cqpm_fusions     to: "${params.outdir}/annotations"
    cqpm_overlap     to: "${params.outdir}/annotations"
    cq_gff3_tidied   to: "${params.outdir}/annotations"
    cqpm_gff3_tidied to: "${params.outdir}/annotations"
}
