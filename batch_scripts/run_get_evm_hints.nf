nextflow.preview.dsl=2

include get_evm_hints from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'

// params.min_intron_hard = false
// params.max_intron_hard = false
// params.table = false

def is_null = { f -> (f == null || f == '') }


workflow {
    main:
    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(Channel.empty(), table)

    } else {
        log.error "Getting evm hints requires a table file."
        exit 1
    }

    hints = get_evm_hints(
        params.min_intron_hard,
        params.max_intron_hard,
        input_channels.spaln_transcripts,
        input_channels.gmap,
        input_channels.spaln_proteins,
        input_channels.exonerate.map { n, f -> [n, "remote", f] }
    )

    publish:
    hints to: "${params.outdir}/evm_hints"
}
