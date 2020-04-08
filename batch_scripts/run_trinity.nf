nextflow.preview.dsl=2

include trinity_assemble_denovo from '../modules/assemblers'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'

params.not_fungus = false
params.table = false


workflow {

    main:
    genomes = Channel.empty()

    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(genomes, table)

    } else {
        log.error "Running trinity requires fastq input."
        exit 1
    }

    assemblies = trinity_assemble_denovo(params.not_fungus, input_channels.fastq)

    publish:
    assemblies to: "${params.outdir}/assemblies"
}
