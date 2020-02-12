nextflow.preview.dsl=2

include get_augustus_hints from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'


params.table = false

def is_null = { f -> (f == null || f == '') }


workflow {
    main:
    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(genomes, table)

    } else {
        log.error "Getting augustus hints required a table file."
        exit 1
    }

    hints = get_augustus_hints(
        input_channels.known,
        input_channels.spaln_transcripts,
        input_channels.spaln_proteins,
        input_channels.exonerate,
        input_channels.genemark,
        input_channels.pasa,
        input_channels.codingquarry,
        input_channels.codingquarrypm,
        input_channels.gemoma,
        input_channels.augustus,
        input_channels.gemoma_comparative
    )


    publish:
    hints to: "${params.outdir}/augustus_hints"
}
