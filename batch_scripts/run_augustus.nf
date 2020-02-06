include run_augustus from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include get_augustus_config from '../modules/utils'
include param_unexpected_error from '../modules/cli'


params.augustus_species = false
params.augustus_utr = false
params.not_fungus = false
params.min_intron_hard = false
params.valid_splicesites = false
params.genomes = false
params.genomes = false
params.table = false
params.augustus_config = false
params.augustus_hint_weights = "data/extrinsic_hints.cfg"

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
        log.error "Running augustus requires hints."
        exit 1
    }

    if ( params.augustus_config ) {
        config_dir = get_file(params.augustus_config)
    } else {
        config_dir = get_augustus_config()
    }

    if ( params.augustus_hint_weights ) {
        hint_weights = get_file(params.augustus_hint_weights)
    } else {
        param_unexpected_error()
    }


    au = run_augustus(
        params.augustus_species,
        params.augustus_utr,
        params.not_fungus,
        params.min_intron_hard,
        params.valid_splicesites,
        input_channels.genome,
        input_channels.known,
        config_dir,
        hint_weights,
        input_channels.augustus_hints
    )

    publish:
    au.augustus_gff3_tidied    to: "${params.outdir}/annotations"
    au.augustus_augustus_hints to: "${params.outdir}/hints"
}
