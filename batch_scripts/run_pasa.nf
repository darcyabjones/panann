#!/usr/bin/env nextflow

nextflow.preview.dsl=2

include run_pasa from '../modules/workflows'
include get_file from '../modules/cli'
include handle_table from '../modules/cli'

params.not_fungus = false
params.max_intron_hard = false
params.genomes = false
params.table = false
params.transcripts = false
params.transcripts_cln = false
params.transcripts_clean = false

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
        transcripts = get_file(params.transcripts)
    } else {
        log.error "Running pasa requires the transcript fasta."
        log.error "This can be found in the output of the " +
                  "run_transcript_alignment pipeline as transcripts.fasta."
        exit 1
    }

    if ( params.transcripts_cln ) {
        transcripts_cln = get_file(params.transcripts_cln)
    } else {
        log.error "Running pasa requires the transcript cln fasta."
        log.error "This can be found in the output of the " +
                  "run_transcript_alignment pipeline as transcripts.fasta.cln."
        exit 1
    }

    if ( params.transcripts_clean ) {
        transcripts_clean = get_file(params.transcripts_clean)
    } else {
        log.error "Running pasa requires the transcript clean fasta."
        log.error "This can be found in the output of the " +
                  "run_transcript_alignment pipeline as transcripts.fasta.clean."
        exit 1
    }

    if ( params.table ) {
        table = get_file(params.table)
        input_channels = handle_table(genomes, table)

    } else {
        log.error "Running pasa requires transcript alignments input."
        exit 1
    }

    ps = run_pasa(
        params.not_fungus,
        params.max_intron_hard,
        Channel.value( [transcripts, transcripts_cln, transcripts_clean] ),
        input_channels.genome,
        input_channels.known,
        input_channels.stringtie,
        input_channels.gmap
    )

    publish:
    ps.pasa_gff3           to: "${params.outdir}/annotations"
    ps.pasa_gff3_tidied    to: "${params.outdir}/annotations"
    ps.pasa_augustus_hints to: "${params.outdir}/hints"
}
