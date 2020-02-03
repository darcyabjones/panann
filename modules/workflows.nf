#!/usr/bin/env nextflow

include get_mmseqs_protein_db as get_mmseqs_remote_protein_db from './aligners'
include mmseqs_search_genome_against_proteins as mmseqs_search_genome_against_remote_proteins from './aligners'
include cluster_genome_vs_protein_matches as cluster_genome_vs_remote_protein_matches from './aligners'
include exonerate_regions as exonerate_remote_proteins from './aligners'
include spaln_align_transcripts from './aligners'
include gmap_align_transcripts from './aligners'

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
include tidy_gff3 as tidy_pasa_gff3 from './utils'
include tidy_gff3 as tidy_codingquarry_gff3 from './utils'
include tidy_gff3 as tidy_codingquarrypm_gff3 from './utils'
include clean_transcripts from './utils'
include combine_fastas from './utils'
include symmetric_difference from './utils'

include extract_augustus_rnaseq_hints from './hints'
include extract_gemoma_rnaseq_hints from './hints'
include extract_augustus_hints as extract_spaln_transcript_augustus_hints from './hints'
include extract_augustus_hints as extract_codingquarry_augustus_hints from './hints'
include extract_augustus_hints as extract_codingquarrypm_augustus_hints from './hints'
include extract_augustus_split_hints as extract_pasa_augustus_hints from './hints'
include extract_spaln_transcript_evm_hints from './hints'
include extract_gmap_evm_hints from './hints'

include get_star_index from './aligners'
include star_find_splicesites from './aligners'
include star_align_reads from './aligners'

def is_null = { f -> (f == null || f == '') }

/**
 * Aligns short-read RNASeq to genomes and get intron hints for
 * augustus and gemoma.
 *
 * @param min_intron_hard An integer value for the minimum intron length
 * @param max_intron_hard An integer value for the maximum intron length
 * @param novel_extra_params String of extra options to provide to STAR
 *                           during first alignment pass.
 * @param align_extra_params String of extra options to provide to STAR
 *                           during second alignment pass.
 * @param valid_splicesites String. comma separated list of valid splice
 *                          boundaries (e.g. "GTAG,CAAG").
 * @param genomes A channel containing fasta genomes to align to.
 *                Should have the structure: tuple(val(name), file("genome.fasta"))
 * @param known A channel containing known gffs for genomes we're aligning to.
 *             Should have the structure: tuple(val(name), file("genome.gff3"))
 * @param fastq A Channel containing fastq reads to align to the genomes.
 *              Should have the structure:
 *              tuple(val(read_group),
 *                    list(file("r1.fastq")),
 *                    list(file("r2.fastq")),
 *                    val(strand))
 * @return A Channel with crams and strands.
 *         tuple(val(name), val(read_group), path(cram), val(strand))
 * @return A Channel with augustus rnaseq hints.
 *         tuple(val(name), val(read_group), path("introns.gff3"))
 * @return A Channel with gemoma rnaseq hints.
 *         tuple(val(name), val(read_group), path("introns.gff3"),
 *               path("forward.bedgraph"), path("reverse.bedgraph"))
 */
workflow align_rnaseq_reads {

    get:
    min_intron_hard
    max_intron_hard
    novel_extra_params
    align_extra_params
    valid_splicesites
    genomes
    known
    fastq

    main:
    // The WAS_NULL business allows us to run where we don't have a known
    // gff dataset.
    index = get_star_index(
        genomes
            .join(known, by: 0, remainder: true)
            .map { n, f, g -> is_null(g) ? [n, f, file('WAS_NULL')]: [n, f, g] }
    )

    // This is the first pass, where we find all splice sites across all
    // read groups.
    splice_sites = star_find_splicesites(
        min_intron_hard,
        max_intron_hard,
        novel_extra_params,
        index.combine(fastq.map { rg, r1, r2, s -> [rg, r1, r2] })
    )

    // Combine the splice sites from the first pass and align the reads
    // in earnest this time.
    aligned = star_align_reads(
        min_intron_hard,
        max_intron_hard,
        align_extra_params,
        index
            .join(genomes, by: 0)
            .combine(fastq)
            .map {n, i, f, rg, r1, r2, s -> [n, rg, f, i, r1, r2] }
            .combine(splice_sites.groupTuple(by: [0, 1]), by: [0, 1])
    )

    aligned_with_strand = aligned
        .map { n, rg, c -> [rg, n, c]}
        .combine(fastq.map { rg, r1, r2, s -> [rg, s] }, by: 0)
        .map { rg, n, c, s -> [n, rg, c, s] }

    // Get the intron hints for augustus.
    augustus_rnaseq_hints = extract_augustus_rnaseq_hints(
        min_intron_hard,
        max_intron_hard,
        valid_splicesites,
        genomes.combine(aligned, by: 0).map { n, f, rg, c -> [n, rg, f, c] }
    )

    // Get the intron and coverage hints for gemoma.
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


/**
 * Assembles transcripts against the genome based on read alignments.
 *
 * @param genomes A channel containing fasta genomes to align to.
 *                Should have the structure: tuple(val(name), file("genome.fasta"))
 * @param known A channel containing known gffs for genomes we're aligning to.
 *             Should have the structure: tuple(val(name), file("genome.gff3"))
 * @param crams A Channel containing cram alignments to the genomes.
 *              Should have the structure:
 *              tuple(val(read_group),
 *                    path("in.cram"),
 *                    val(strand))
 * @return A Channel with the final merged gtf files.
 *         tuple(val(name), val(gtf))
 * @return A Channel with individual gtf assemblies for each read_group.
 *         tuple(val(name), val(read_group), path(gtf))
 */
workflow run_stringtie {

    get:
    genomes
    known
    crams

    main:

    known_with_null = genomes
        .join(known, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('WAS_NULL')]: [n, g] }

    individually_assembled = stringtie_assemble(
        genomes
            .combine(aligned, by: 0)
            .map { n, f, rg, c -> [n, rg, f, c, s]}
            .combine(known_with_null, by: 0)
    )

    stringtie = stringtie_merge(
        individually_assembled
            .map { n, rg, g -> [n, g] }
            .groupTuple(by: 0)
            .join(known_with_null, by: 0)
    )

    emit:
    stringtie
    individually_assembled
}


/**
 * Aligns proteins from potentially distantly related organsisms to the genomes.
 *
 * Example:
 * > align_remote_proteins(5, 20000, 20000, 1, remote_proteins,
 * >                       genomes, genome_faidxs,
 * >                       genome_mmseqs_indices)
 *
 * @param min_intron_hard An integer value for the minimum intron length
 * @param max_intron_hard An integer value for the maximum intron length
 * @param max_gene_hard An integer value for the maximum gene length
 * @param trans_table The numeric id of the NCBI translation table to use.
 * @param remote_proteins The protein set(s) to align to the genomes.
 *                        tuple(val(name), path("proteins.faa"))
 * @param genomes The genomes to align the proteins to.
 *                tuple(val(name), path("genome.fasta"))
 * @param genome_faidxs The genome fasta indices.
 *                tuple(val(name), path("genome.fasta.fai"))
 * @param genome_mmseqs_indices The mmseqs indexed genomes.
 * @return A channel containing mmseqs protein matches as a tsv file.
 * @return A channel containing bedfiles showing matches in a genomic
 *         context.
 * @return The exonerate GFF2 matches.
 */
workflow align_remote_proteins {

    get:
    min_intron_hard
    max_intron_hard
    max_gene_hard
    trans_table
    remote_proteins
    genomes
    genome_faidxs
    genome_mmseqs_indices

    main:
    // Get the mmseqs index for the remote protein dataset.
    index = get_mmseqs_remote_protein_db(remote_proteins)

    // Convert the fasta to a tsv, which we use later to realign using exonerate.
    tsv = remote_proteins_fasta_to_tsv(remote_proteins)

    // Search the proteins against the genome, and get a tsv.
    mmseqs_matches = mmseqs_search_genome_against_remote_proteins(
        trans_table,
        genome_mmseqs_indices.combine(index)
    )

    // Find regions of genome to align to and proteins that should
    // be realigned to those regions.
    clustered = cluster_genome_vs_remote_protein_matches(
        max_gene_hard,
        1000,
        genome_faidxs.join(mmseqs_matches, by: 0)
    )

    // Align each subset of proteins to the genome regions that
    // they matched with the less sensitive methods.
    exonerate_matches = exonerate_remote_proteins(
        trans_table,
        min_intron_hard,
        max_intron_hard,
        genomes
            .join(clustered, by: 0)
            .combine(tsv.map { n, r -> r })
    )

    exonerate_augustus_hints = extract_exonerate_hints(
        min_intron_hard,
        max_intron_hard,
        exonerate_matches
    )

    exonerate_evm_hints = extract_exonerate_evm_hints(
        min_intron_hard,
        max_intron_hard,
        exonerate_matches
    )

    // Do hint extraction here?

    emit:
    mmseqs_matches
    clustered
    exonerate_matches
    exonerate_augustus_hints
    exonerate_evm_hints
}


/**
 * Align transcripts to the genomes using spaln and gmap.
 *
 * @param species The "species" to specify for spaln splice site motifs.
 * @param min_intron_soft The rough minimum length an intron should be,
 *                        some introns may still be shorter than this.
 * @param min_intron_hard The minimum length an intron should be.
 * @param max_intron_hard The maximum length an intron should be.
 * @param max_gene_hard The maximum length a gene should be.
 * @param univec A value channel containing the univec fasta file.
 * @param transcripts A channel containing transcript fasta files.
 * @param spaln_indices Spaln formatted databases.
 * @param gmap_indices gmap formatted databases.
 */
workflow align_transcripts {

    get:
    species
    min_intron_soft
    min_intron_hard
    max_intron_hard
    max_gene_hard
    univec
    transcripts
    spaln_indices
    gmap_indices

    main:
    (combined_fasta, combined_tsv) = combine_fastas(transcripts.collect())
    cleaned_transcripts = clean_transcripts(combined_fasta)

    spaln_aligned = spaln_align_transcripts(
        species,
        max_gene_hard,
        min_intron_soft,
        spaln_indices.combine(cleaned_transcripts.map { f, cln, clean -> clean })
    )

    spaln_augustus_hints = extract_spaln_transcript_augustus_hints(
        "spaln_transcripts",
        "spaln",
        "E",
        3,
        6,
        0,
        0,
        0,
        false,
        spaln_aligned
    )

    spaln_evm_hints = extract_spaln_transcript_evm_hints(spaln_aligned)

    gmap_aligned = gmap_align_transcripts(
        min_intron_hard,
        max_intron_hard,
        gmap_indices.combine(cleaned_transcripts.map { f, cln, clean -> clean })
    )

    gmap_evm_hints = extract_gmap_evm_hints(gmap_aligned)

    emit:
    combined_fasta
    combined_tsv
    cleaned_transcripts
    spaln_aligned
    spaln_augustus_hints
    spaln_evm_hints
    gmap_aligned
    gmap_evm_hints
}


/**
 * Run PASA/transdecoder
 *
 * @param not_fungus bool, don't use fungal optimised parameters.
 * @param max_intron_hard the maximum length that an intron should be.
 * @param transcripts the cleaned transcripts that were aligned with gmap.
 *                    structure tuple("in.fasta", "in.fasta.cln", "in.fasta.clean")
 * @param genomes The fasta genomes to predict genes in.
 * @param known GFF3s for known sites in each genome (if any).
 * @param stringtie Stringtie gtfs to use as evidence.
 * @param gmap Gmap gff3 alignments to the genomes.
 */
workflow run_pasa {

    get:
    not_fungus
    max_intron_hard
    transcripts
    genomes
    known
    stringtie
    gmap

    main:
    known_with_null = genomes
        .join(known, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('KNOWN_WAS_NULL')]: [n, g] }

    stringtie_with_null = genomes
        .join(stringtie, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('STRINGTIE_WAS_NULL')]: [n, g] }

    genome_names = genomes.map { n, f -> n }.toList()
    gmap_names = gmap.map { n, g -> n }.toList()
    sym_diff = symmetric_difference(genome_names, gmap_names)
    if ( sym_diff.length != 0 ) {
        log.error "Some names in the genomes and gmap files don't match up."
        log.error "They are: {sym_diff}"
        exit 1
    }

    pasa_gff3 = pasa(
        not_fungus,
        max_intron_hard,
        genomes
            .join(known_with_null, by: 0)
            .join(stringtie_with_null, by: 0)
            .join(gmap_aligned, by: 0)
            .combine(transcripts)
    )

    pasa_gff3_tidied = tidy_pasa_gff3(
        "pasa",
        "pasa",
        pasa_gff3
    )

    pasa_augustus_hints = extract_pasa_augustus_hints(
        "pasa",
        "pasa",
        "transdecoder",
        "PASA",
        "TD",
        4, // exon_priority
        3, // CDS priority
        9, // exon_trim
        9, // cds_trim
        6, // utr_trim,
        6, // gene_trim
        false,
        pasa_gff3_tidied
    )

    emit:
    pasa_gff3
    pasa_gff3_tidied
    pasa_augustus_hints
}


/**
 * Run codingquarry and codingquarrypm.
 *
 * @param signalp A boolean value, whether to use signalp or not (deepsig).
 * @param genomes A channel containing the genomes to predict genes for.
 *                Should have the structure: tuple(val(name), path("genome.fasta"))
 * @param stringtie A channel containing the stringtie assemblies.
 *                  Should have the structure: tuple(val(name), path("str.gtf"))
 * @return cq_gff3 The raw codingquarry gff3 file (warts and all).
 * @return cq_faa The predicted proteins from cq.
 * @return cq_fna " for CDSs
 * @return cq_dubious A set of dubious gene predictions.
 * @return cq_fusions A list of stringtie assemblies that were split by cq.
 * @return cq_overlap The cq overlap report.
 * @return cqpm_gff3 The raw codingquarrypm gff3 file (wards and all).
 * @return cqpm_fusions
 * @return cqpm_overlap
 * @return cq_gff3_tidied A tidied gff3 of codingquarry predictions.
 * @return cqpm_gff3_tidied
 */
workflow run_codingquarry {

    get:
    signalp
    genomes
    stringtie

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

    cqpm_augustus_hints = extract_codingquarry_augustus_hints(
        "codingquarry",
        "codingquarry",
        "CQ",
        4, // priority
        6, // exon_trim
        6, // cds_trim
        6, // utr trim
        6, // gene trim
        false,
        cq_gff3_tidied
    )

    cqpm_gff3_tidied = tidy_codingquarrypm_gff3(
        "codingquarrypm",
        "codingquarrypm",
        cqpm_gff3
    )

    cqpm_augustus_hints = extract_codingquarrypm_augustus_hints(
        "codingquarrypm",
        "codingquarrypm",
        "CQPM",
        4, // priority
        6, // exon_trim
        6, // cds_trim
        6, // utr trim
        6, // gene trim
        false,
        cqpm_gff3_tidied
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
    cq_augustus_hints
    cqpm_augustus_hints
}


/**
 * Runs the comparative gemoma pipeline.
 *
 * @param trans_table The integer corresponding to the NCBI translation table.
 * @param genomes A channel of genomes to run. Structure: tuple(val(name), path("in.fasta"))
 * @param introns A channel of intron hints to use for gemoma.
 *                Structure: tuple(val(name), path("hints.gff3"))
 * @param pasa PASA/transdecoder predictions in gff3 format.
 * @param cq Codingquarry predictions in gff3 format.
 * @param cqpm CodingQuarryPM predictions in gff3 format.
 * @param augustus Augustus predictions in gff3 format.
 * @return gemoma_predictions The raw genoma predictions.
 *         Structure tuple(val(name), path(gff3)).
 * @return gemoma_predictions_tidied The tidied gemoma predictions.
 *         Structure tuple(val(name), path(gff3)).
 */
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
