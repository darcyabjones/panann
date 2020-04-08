nextflow.preview.dsl=2

include align_rnaseq_reads from '../modules/workflows'
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
        log.error "Running star requires fastq input."
        exit 1
    }

    (crams, aug_hints, gemoma_hints) = align_rnaseq_reads(
        params.min_intron_hard,
        params.max_intron_hard,
        params.star_novel_params,
        params.star_align_params,
        params.valid_splicesites,
        input_channels.genomes,
        input_channels.known,
        input_channels.fastq
    )

    publish:
    crams to: "${params.outdir}/alignments"
    aug_hints to: "${params.outdir}/hints"
    gemoma_hints to: "${params.outdir}/hints"
}
