nextflow.preview.dsl=2

include run_evm from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include get_faidx from '../modules/utils'
include get_augustus_config from '../modules/utils'
include param_unexpected_error from '../modules/cli'


params.augustus_species = false
params.augustus_utr = false
params.not_fungus = false
params.min_intron_hard = false
params.valid_splicesites = false
params.genomes = false
params.table = false
params.augustus_config = false
params.augustus_gapfiller_weights = "data/extrinsic_gapfiller.cfg"
params.evm_config = "data/evm.cfg"

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

    if ( params.augustus_gapfiller_weights ) {
        augustus_hint_weights = get_file(params.augustus_gapfiller_weights)
    } else {
        param_unexpected_error()
    }

    if ( params.evm_config ) {
        evm_weights = get_file(params.evm_config)
    } else {
        param_unexpected_error()
    }

    evm_gene_hints = input_channels.pasa
        .mix(
            input_channels.genemark,
            input_channels.codingquarry,
            input_channels.codingquarrypm,
            input_channels.augustus,
            input_channels.gemoma,
            input_channels.gemoma_comparative,
        )

    genome_faidxs = get_faidx(input_channels.genome)

    evm = run_evm(
        params.augustus_species,
        params.not_fungus,
        params.augustus_utr,
        params.valid_splicesites,
        params.min_intron_hard,
        augustus_config_dir,
        augustus_hint_weights,
        evm_weights,
        input_channels.genome,
        genome_faidxs,
        input_channels.evm_transcript_hints,
        input_channels.evm_protein_hints,
        evm_gene_hints,
        input_channels.known,
        input_channels.augustus_intron_hints.mix(input_channels.augustus_hints)
    )

    publish:
    evm.evm_gff3_tidied to: "${params.outdir}/annotations"
    evm.augustus_gff3_tidied to: "${params.outdir}/annotations"
    evm.final_gff3 to: "${params.outdir}/annotations"
}
