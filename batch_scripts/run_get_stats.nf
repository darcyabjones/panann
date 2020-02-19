nextflow.preview.dsl=2

include run_stats from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include get_augustus_config from '../modules/utils'
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
        log.error "Running evidence modeller and the augustus gapfiller requires hints."
        exit 1
    }

    if ( params.augustus_config ) {
        augustus_config_dir = get_file(params.augustus_config)
    } else {
        augustus_config_dir = get_augustus_config()
    }

    if ( params.busco_lineage ) {
        busco_lineage = get_file(params.busco_lineage)
    } else {
        busco_lineage = Channel.empty()
    }

    gffs = input_channels.pasa.map { n, f -> [n, "pasa", f] } .mix(
        input_channels.genemark.map { n, f -> [n, "genemark", f] },
        input_channels.codingquarry.map { n, f -> [n, "codingquarry", f] },
        input_channels.codingquarrypm.map { n, f -> [n, "codingquarrypm", f] },
        input_channels.augustus.map { n, f -> [n, "augustus", f] },
        input_channels.gemoma.map { n, f -> [n, "gemoma", f] },
        input_channels.gemoma_comparative.map { n, f -> [n, "gemoma_comparative", f] },
        input_channels.evm.map { n, f -> [n, "evm", f] },
        input_channels.augustus_gapfiller.map { n, f -> [n, "augustus_gapfiller", f] },
        input_channels.complete.map { n, f -> [n, "final", f] },
    )

    stats = run_stats(
        params.trans_table,
        input_channels.genome,
        gffs,
        input_channels.known,
        busco_lineage,
        augustus_config_dir,
    )

    publish:
    stats to: "${params.outdir}/stats"
}
