nextflow.preview.dsl=2

include align_transcripts from '../modules/workflows'
include get_spaln_index from '../modules/aligners'
include get_gmap_index from '../modules/aligners'

include get_file from '../modules/cli'
include handle_table from '../modules/cli'
include get_univec from '../modules/utils'

params.species = false
params.min_intron_soft = false
params.min_intron_hard = false
params.max_intron_hard = false
params.max_gene_hard = false
params.univec = false
params.transcripts = false
params.transcripts_cln = false
params.transcripts_clean = false
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

    if ( params.transcripts ) {
        transcripts = Channel.value(get_file(params.transcripts))
    } else {
        log.error "Running transcript alignments requires the transcript fasta."
        exit 1
    }

    if ( params.univec ) {
        univec = Channel.value( get_file(params.univec) )
    } else {
        univec = get_univec()
    }

    spaln_indices = get_spaln_index(genomes)
    gmap_indices = get_gmap_index(genomes)

    aligned_transcripts = align_transcripts(
        params.species,
        params.min_intron_soft,
        params.min_intron_hard,
        params.max_intron_hard,
        params.max_gene_hard,
        univec,
        transcripts,
        spaln_indices,
        gmap_indices
    )

    publish:
    aligned_transcripts.combined_fasta to: "${params.outdir}/assemblies"
    aligned_transcripts.combined_tsv to: "${params.outdir}/assemblies"
    aligned_transcripts.cleaned_transcripts to: "${params.outdir}/assemblies"
    aligned_transcripts.spaln_aligned to: "${params.outdir}/alignments"
    aligned_transcripts.spaln_augustus_hints to: "${params.outdir}/hints"
    aligned_transcripts.spaln_evm_hints to: "${params.outdir}/hints"
    aligned_transcripts.gmap_aligned to: "${params.outdir}/alignments"
    aligned_transcripts.gmap_evm_hints to: "${params.outdir}/hints"
}
