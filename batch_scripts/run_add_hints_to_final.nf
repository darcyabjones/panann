nextflow.preview.dsl=2

include filter_preds from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include param_unexpected_error from '../modules/cli'


//params.augustus_species = false
//params.augustus_utr = false
//params.not_fungus = false
//params.min_intron_hard = false
//params.valid_splicesites = false
//params.genomes = false
//params.table = false
//params.augustus_config = false
//params.augustus_gapfiller_weights = "data/extrinsic_gapfiller.cfg"
//params.evm_config = "data/evm.cfg"

def is_null = { f -> (f == null || f == '') }


workflow {
    main:

    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(Channel.empty(), table)

    } else {
        log.error "Need a table to filter predictions."
        exit 1
    }


    filtered = filter_preds(
        input_channels.complete,
        input_channels.spaln_transcripts,
        input_channels.spaln_proteins,
        input_channels.gmap,
        input_channels.exonerate,
        input_channels.genemark,
        input_channels.pasa,
        input_channels.codingquarry,
        input_channels.codingquarrypm,
        input_channels.gemoma,
        input_channels.augustus,
        input_channels.gemoma_comparative,
    )

    publish:
    filtered to: "${params.outdir}/final"
}
