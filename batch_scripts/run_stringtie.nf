nextflow.preview.dsl=2

include run_stringtie from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'

params.genomes = false
params.table = false

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
        log.error "Running stringtie requires cram input."
        exit 1
    }

    (combined, indiv) = run_stringtie(
        input_channels.genome,
        input_channels.known,
        input_channels.cram
    )

    publish:
    combined to: "${params.outdir}/assemblies"
    indiv to: "${params.outdir}/assemblies"
}
