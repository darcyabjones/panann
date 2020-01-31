#!/usr/bin/env nextflow

include get_mmseqs_protein_db as get_mmseqs_remote_protein_db from './aligners'
include mmseqs_search_genome_against_proteins as mmseqs_search_genome_against_remote_proteins from './aligners'
include cluster_genome_vs_protein_matches as cluster_genome_vs_remote_protein_matches from './aligners'
include exonerate_regions as exonerate_remote_proteins from './aligners'

include codingquarry from './predictors'
include codingquarrypm from './predictors'
include signalp from './predictors'
include deepsig from './predictors'

include extract_gemoma_comparative_cds_parts from './predictors'
include cluster_gemoma_cds_parts from './predictors'
include gemoma as gemoma_comparative from './predictors'
include gemoma_combine as gemoma_comparative_combine from './predictors'

include fasta_to_tsv as remote_proteins_fasta_to_tsv from './utils'
include tidy_gff3 as tidy_gemoma_gff3 from './utils'
include tidy_gff3 as tidy_codingquarry_gff3 from './utils'
include tidy_gff3 as tidy_codingquarrypm_gff3 from './utils'

include extract_augustus_rnaseq_hints from './hints'
include extract_gemoma_rnaseq_hints from './hints'
include get_star_index from './aligners'
include star_find_splicesites from './aligners'
include star_align_reads from './aligners'

is_null = { f -> (f == null || f == '') }


valid_table_analyses = [
    "genome",
    "known",
    "fastq_forward",
    "fastq_reverse",
    "gmap",
    "spaln_transcripts",
    "spaln_proteins",
    "exonerate",
    "cram",
    "stringtie",
    "augustus_intron_hints",
    "genemark",
    "pasa",
    "codingquarry",
    "codingquarrypm",
    "gemoma_intron_hints",
    "gemoma_forward_coverage",
    "gemoma_reverse_coverage",
    "gemoma",
    "augustus",
    "gemoma_comparative",
    "evm",
    "augustus_gapfiller"
]


def param_unexpected_error() {
    log.error "We encountered an error while validating input arguments that " +
        "should be possible. Please raise an issue on github or contact the " +
        "authors."
    exit 1
}


def validate_params(params) {
    def should_run = [:]

    def error = false

    if ( false ) {
        log.error "message"
        error = true
    }

    if (error) {
        exit 1
    }

    return should_run
}


def get_file(filepath) {
    if ( filepath ) {
        handle = Channel.value( file(filepath, checkIfExists: true) )
    } else {
        // The expectation is that you would check that filepath is not false
        // before you use this.
        param_unexpected_error()
    }

    return handle
}


def table_not_null(inval, field, line) {
    if (is_null(inval)) {
        log.error "A line supplied to `--table` had no ${field} field specified."
        log.error "This value is required."
        log.error "The line was: ${line}"

        exit 1
    }

    return inval
}


def table_get_file(filepath, line) {
    if ( filepath ) {
        handle = file(filepath, checkIfExists: true)
    } else {
        log.error 'A line supplied to `--table` had no file field.'
        log.error 'The line was: ${it}'

        exit 1
    }

    return handle
}


def table_one_of(inval, options, field, line) {
    if ( is_null ) {
        return inval
    }

    if ( inval in options ) {
        return inval
    } else {
        log.error "A line supplied to `--table` used an invalid choice for the ${field} field."
        log.error "The line was: ${line}"

        exit 1
    }
}


def handle_table(genomes, table) {

    branched = table
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { it ->
            def line = [:]

            // Analysis is reqiured.
            line.analysis = table_one_of(
                table_not_null(it.analysis, "analysis", it),
                valid_table_analyses,
                "analysis",
                it
            )

            // Name is required for all analyses except fastq files.
            if ( line.analysis in ["fastq_forward", "fastq_reverse"] ) {
                line.name = it.name
            } else {
                line.name = table_not_null(it.name, "name", it)
            }

            line.file = table_get_file(
                table_not_null(it.file, "file", it),
                it
            )

            if ( line.analysis in ["fastq_forward", "fastq_reverse", "cram"]) {
                line.strand = table_one_of(
                    table_not_null(it.strand, "strand", it),
                    ["fr", "rf"],
                    "strand",
                    it
                )
            }

            if ( line.analysis in ["fastq_forward", "fastq_reverse"]) {
                line.read_group = table_not_null(it.read_group, "read_group", it)
            }

            return line
        }
        .branch {
            genome: it.analysis == "genome"
            known: it.analysis == "known"
            fastq_forward: it.analysis == "fastq_forward"
            fastq_reverse: it.analysis == "fastq_reverse"
            gmap: it.analysis == "gmap"
            spaln_transcripts: it.analysis == "spaln_transcripts"
            spaln_proteins: it.analysis == "spaln_proteins"
            exonerate: it.analysis == "exonerate"
            cram: it.analysis == "cram"
            stringtie: it.analysis == "stringtie"
            augustus_intron_hints: it.analysis == "augustus_intron_hints"
            genemark: it.analysis == "genemark"
            pasa: it.analysis == "pasa"
            codingquarry: it.analysis == "codingquarry"
            codingquarrypm: it.analysis == "codingquarrypm"
            gemoma_intron_hints: it.analysis == "gemoma_intron_hints"
            gemoma_forward_coverage: it.analysis == "gemoma_forward_coverage"
            gemoma_reverse_coverage: it.analysis == "gemoma_reverse_coverage"
            gemoma: it.analysis == "gemoma"
            augustus: it.analysis == "augustus"
            gemoma_comparative: it.analysis == "gemoma_comparative"
            evm: it.analysis == "evm"
            augustus_gapfiller: it.analysis == "augustus_gapfiller"
        }


    // At this point we know that files exist and no fields that should be filled
    // are null.
    table_genomes = genomes.mix(branched.genome.map { [it.name, it.file] })

    known = branched.known.map { [it.name, it.file] }

    fastq_strands = branched.fastq_forward
        .mix(branched.fastq_reverse)
        .map { [it.read_group, it.strand] }
        .unique()
        .groupTuple(by: 0)
        .map { rg, strands ->
            if (strands.size() > 1) {
                log.error "The fastq read group ${rg} in the table had more than one strandedness specified."
                log.error "Please make sure they are all the same or split them into multiple read groups."
                exit 1
            }

            return [rg, strands[0]]
        }

    fastq = branched.fastq_forward
        .map { [it.read_group, it.file] }
        .groupTuple(by: 0)
        .join(
            branched.fastq_reverse
                .map { [it.read_group, it.file] }
                .groupTuple(by: 0),
            by: 0
        )
        .combine(fastq_strands, by: 0)

    cram = branched.cram.map { [it.name, it.file, it.strand] }

    out_channels = ["genomes": table_genomes, "fastq": fastq, "cram": cram]
    valid_table_analyses.each { analysis ->
        if ( !(analysis in ["genomes", "fastq", "cram"]) ) {
            out_channels[analysis] = branched
                .getProperty(analysis)
                .map { [it.name, it.file] }
        }
    }

    return out_channels
}


workflow align_rnaseq_reads {

    get:
    min_intron_len  // integer, minimum allowed intron length
    max_intron_len  // integer, maximum allowed intron length
    novel_extra_params  // String of extra options to provide to STAR
    align_extra_params  // String of extra options to provide to STAR
    valid_splicesites
    genomes  // tuple(val(name), path(fasta))
    gffs  // tuple(val(name), path(gff))
    fastq  // tuple(val(read_group), path(list(r1)), path(list(r2)), val(strand))

    main:
    index = get_star_index(genomes.join(gffs, by: 0))
    splice_sites = star_find_splicesites(
        min_intron_len,
        max_intron_len,
        novel_extra_params,
        index.combine(fastq.map { rg, r1, r2, s -> [rg, r1, r2] })
    )

    aligned = star_align_reads(
        min_intron_len,
        max_intron_len,
        align_extra_params,
        index
            .combine(fastq)
            .map {n, i, rg, r1, r2, s -> [n, rg, i, r1, r2] }
            .combine(splice_sites.groupTuple(by: [0, 1]), by: [0, 1])
    )

    aligned_with_strand = aligned
        .map { n, rg, c -> [rg, n, c]}
        .combine(fastq.map { rg, r1, r2, s -> [rg, s] }, by: 0)
        .map { rg, n, c, s -> [n, rg, c, s] }

    augustus_rnaseq_hints = extract_augustus_rnaseq_hints(
        min_intron_len,
        max_intron_len,
        valid_splicesites,
        genomes.combine(aligned, by: 0).map { n, f, rg, c -> [n, rg, f, c] }
    )

    gemoma_rnaseq_hints = extract_gemoma_rnaseq_hints(
        genomes
            .combine(aligned_with_strand, by: 0)
            .map { n, f, rg, c, s -> [n, rg, f, c, s] }
    )

    emit:
    aligned_with_strand
    augustus_rnaseq_hints
    gemoma_rnaseq_hints
}


workflow run_stringtie {

    get:
    genomes
    aligned

    main:
    individually_assembled = stringtie_assemble(
        genomes
            .combine(aligned, by: 0)
            .map { n, f, rg, c -> [rg, n, f, c]}
            .combine(fastq.map { rg, r1, r2, s -> [rg, s] }, by: 0)
            .map { rg, n, f, c, s -> [n, rg, f, c, s] }
            .combine(gffs, by: 0)
    )

    stringtie = stringtie_merge(
        individually_assembled
            .map { n, rg, g -> [n, g] }
            .groupTuple(by: 0)
            .join(gffs)
    )

    emit:
    stringtie
    individually_assembled
}


workflow align_remote_proteins {

    // align_remote_proteins(5, 20000, 20000, 1, remote_proteins,
    //                       genomes, genome_faidxs,
    //                       genome_mmseqs_indices)

    get:
    min_intron_hard  // An integer value channel for the minimum intron length
    max_intron_hard  // Integer for max intron length
    max_gene_hard  // Absolute Max gene length
    trans_table  // The numeric translation table for proteins.
    remote_proteins  // tuple(val(name), path("proteins.fasta"))
    genomes  // tuple(val(name), path("genome.fasta"))
    genome_faidxs  // tuple(val(name), path("genome.fasta.fai"))
    genome_mmseqs_indices  // tuple(val(name), path("genome_db"))

    main:
    index = get_mmseqs_remote_protein_db(remote_proteins)
    tsv = remote_proteins_fasta_to_tsv(remote_proteins)

    mmseqs_matches = mmseqs_search_genome_against_remote_proteins(
        trans_table,
        genome_mmseqs_indices.combine(index)
    )

    clustered = cluster_genome_vs_remote_protein_matches(
        max_gene_hard,
        1000,
        genome_faidxs.join(mmseqs_matches, by: 0)
    )

    exonerate_matches = exonerate_remote_proteins(
        trans_table,
        min_intron_hard,
        max_intron_hard,
        genomes
            .join(clustered, by: 0)
            .combine(tsv.map { n, r -> r })
    )

    emit:
    mmseqs_matches
    clustered
    exonerate_matches
}


workflow run_codingquarry {

    get:
    signalp // bool
    genomes // tuple(val(name), path(fasta))
    stringtie // tuple(val(name), path(gtf))

    main:
    (cq_fixed_gff3, cq_gff3, cq_faa, cq_fna, cq_dubious, cq_fusions, cq_overlap) = codingquarry(
        stringtie.join(genomes)
    )

    if (signalp) {
        cq_secreted = signalp(cq_faa)
    } else {
        cq_secreted = deepsig(cq_faa)
    }

    (cqpm_gff3, cqpm_fusions, cqpm_overlap) = codingquarrypm(
        stringtie.join(genomes).join(cq_gff3).join(cq_secreted)
    )


    cq_gff3_tidied = tidy_codingquarry_gff3(
        "codingquarry",
        "codingquarry",
        cq_fixed_gff3
    )

    cqpm_gff3_tidied = tidy_codingquarrypm_gff3(
        "codingquarrypm",
        "codingquarrypm",
        cqpm_gff3
    )

    emit:
    cq_gff3
    cq_faa
    cq_fna
    cq_dubious
    cq_fusions
    cq_overlap
    cqpm_gff3
    cqpm_fusions
    cqpm_overlap
    cq_gff3_tidied
    cqpm_gff3_tidied
}


workflow run_gemoma_comparative {

    get:
    trans_table  // val, should be ncbi table integer.
    genomes // tuple(val(name), file(fasta))
    introns  // tuple(val(name), file(introns.gff))
    pasa
    cq
    cqpm
    augustus

    main:
    gffs = pasa
        .map {n, f -> [n, "pasa", f]}
        .mix(
            cq.map {n, f -> [n, "codingquarry", f]},
            cqpm.map {n, f -> [n, "codingquarrypm", f]},
            augustus.map {n, f -> [n, "augustus", f]}
        )


    cds_parts = extract_gemoma_comparative_cds_parts(
        genomes
            .combine(gffs, by: 0)
            .map { n, f, a, g -> [n, a, f, g] }
    )


    clustered_cds_parts = cluster_gemoma_cds_parts(
        cds_parts
            .map {n, an, c, a, p -> [c, a, p]}
            .collect()
    )

    mmseqs_search_gemoma_cds_parts(
        trans_table,
        clustered_cds_parts
            .map { c, a, p -> ["comparative", c, a, p] }
            .combine(genomes, by: 0)
    )

    gemoma_matches = gemoma(
        genomes.combine(
            clustered_cds_parts.map {c, a, p -> ["comparative", c, a, p]},
            by: 0
        )
        .map { t, f, r, c, a, p -> [t, r, f, c, a, p] }
        .join(mmseqs_matches, by: 0)
        .join(introns, by: 0)
    )

    gemoma_predictions = gemoma_combine(
        gemoma_matches
            .groupTuple(by: 0)
            .join(genomes, by: 0)
            .join(introns, by: 0)
    )

    gemoma_predictions_tidied = tidy_gemoma(
        "gemoma",
        "gemoma",
        gemoma_predictions
    )

    emit:
    gemoma_predictions
    gemoma_predictions_tidied
}
