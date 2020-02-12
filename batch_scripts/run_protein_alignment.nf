nextflow.preview.dsl=2

include align_proteins from '../modules/workflows'
include get_spaln_index from '../modules/aligners'

include get_file from '../modules/cli'
include handle_table from '../modules/cli'

params.trans_table = false
params.min_intron_soft = false
params.max_gene_hard = false
params.proteins
params.genomes = false


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

    if ( params.proteins ) {
        proteins = Channel
            .fromPath(params.proteins, checkIfExists: true, type: "file")
    } else {
        log.error "Running protein alignments requires the protein fasta(s)."
        exit 1
    }

    spaln_indices = get_spaln_index(genomes)

    aligned_proteins = align_proteins(
        params.trans_table,
        params.min_intron_soft,
        params.max_gene_hard,
        proteins,
        spaln_indices,
    )

    publish:
    aligned_proteins.combined_fasta to: "${params.outdir}/proteins"
    aligned_proteins.combined_tsv to: "${params.outdir}/proteins"
    aligned_proteins.spaln_aligned_tidied to: "${params.outdir}/alignments"
    aligned_proteins.spaln_augustus_hints to: "${params.outdir}/hints"
    aligned_proteins.spaln_evm_hints to: "${params.outdir}/hints"
}
