nextflow.preview.dsl=2

include run_gemoma from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include get_mmseqs_genome_db from '../modules/aligners'


params.trans_table = false
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
        log.error "Running gemoma requires some known sites and introns input."
        exit 1
    }

    genome_mmseqs_indices = get_mmseqs_genome_db(input_channels.genome)

    gm = run_gemoma(
        params.trans_table,
        input_channels.genome,
        genome_mmseqs_indices,
        input_channels.gemoma_intron_hints,
        input_channels.known
    )

    publish:
    gm.gemoma_gff3 to: "${params.outdir}/annotations"
    gm.gemoma_gff3_tidied to: "${params.outdir}/annotations"
    gm.gemoma_augustus_hints to: "${params.outdir}/hints"
}
