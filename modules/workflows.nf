#!/usr/bin/env nextflow

include get_mmseqs_protein_db as get_mmseqs_remote_protein_db from './aligners'
include mmseqs_search_genome_against_proteins as mmseqs_search_genome_against_remote_proteins from './aligners'
include cluster_genome_vs_protein_matches as cluster_genome_vs_remote_protein_matches from './aligners'
include exonerate_regions as exonerate_remote_proteins from './aligners'
include spaln_align_transcripts from './aligners'
include gmap_align_transcripts from './aligners'
include spaln_align_proteins from './aligners'
include fix_spaln_proteins_stop from './aligners'
include get_star_index from './aligners'
include star_find_splicesites from './aligners'
include star_align_reads from './aligners'

include codingquarry from './predictors'
include codingquarrypm from './predictors'
include signalp from './predictors'
include deepsig from './predictors'
include pasa from './predictors'

include cluster_gemoma_cds_parts from './predictors'
include extract_gemoma_cds_parts from './predictors'
include extract_gemoma_comparative_cds_parts from './predictors'
include mmseqs_search_gemoma_cds_parts from './predictors'
include mmseqs_search_gemoma_cds_parts as mmseqs_search_gemoma_comparative_cds_parts from './predictors'
include gemoma from './predictors'
include gemoma as gemoma_comparative from './predictors'
include gemoma_combine from './predictors'
include gemoma_combine as gemoma_comparative_combine from './predictors'
include augustus_hints from './predictors'
include evm from './predictors'
include find_missing_evm_predictions from './predictors'
include augustus_gap_filler from './predictors'

include stringtie_assemble from './assemblers'
include stringtie_merge from './assemblers'

include fasta_to_tsv as remote_proteins_fasta_to_tsv from './utils'
include tidy_gff3 as tidy_gemoma_gff3 from './utils'
include tidy_gff3 as tidy_gemoma_comparative_gff3 from './utils'
include tidy_gff3 as tidy_pasa_gff3 from './utils'
include tidy_gff3 as tidy_codingquarry_gff3 from './utils'
include tidy_gff3 as tidy_codingquarrypm_gff3 from './utils'
include tidy_gff3 as tidy_known_gff3 from './utils'
include tidy_gff3 as tidy_evm_gff3 from './utils'
include tidy_gff3 as tidy_spaln_transcripts_gff3 from './utils'
include tidy_gff3 as tidy_spaln_proteins_gff3 from './utils'
include combine_and_tidy_gff3 as combine_and_tidy_augustus_gff3 from './utils'
include clean_transcripts from './utils'
include combine_fastas from './utils'
include chunkify_genomes from './utils'
include merge_gffs from './utils'
include extract_seqs from './utils'
include exonerate_to_gff3 from './utils'
include gff_to_bed as exonerate_gff_to_bed from './utils'
include gff_to_bed as spaln_transcripts_gff_to_bed from './utils'
include gff_to_bed as spaln_proteins_gff_to_bed from './utils'
include gff_to_bed as gmap_transcripts_gff_to_bed from './utils'
include gff_to_bed as pasa_gff_to_bed from './utils'
include gff_to_bed as genemark_gff_to_bed from './utils'
include gff_to_bed as cq_gff_to_bed from './utils'
include gff_to_bed as cqpm_gff_to_bed from './utils'
include gff_to_bed as gemoma_gff_to_bed from './utils'
include gff_to_bed as augustus_gff_to_bed from './utils'
include gff_to_bed as gemoma_comparative_gff_to_bed from './utils'
include get_hint_coverage from './utils'

include extract_augustus_rnaseq_hints from './hints'
include extract_gemoma_rnaseq_hints from './hints'
include combine_gemoma_rnaseq_hints from './hints'
include extract_augustus_hints as extract_spaln_transcript_augustus_hints from './hints'
include extract_augustus_hints as extract_spaln_protein_augustus_hints from './hints'
include extract_augustus_hints as extract_genemark_augustus_hints from './hints'
include extract_augustus_hints as extract_codingquarry_augustus_hints from './hints'
include extract_augustus_hints as extract_codingquarrypm_augustus_hints from './hints'
include extract_augustus_hints as extract_known_augustus_hints from './hints'
include extract_augustus_hints as extract_augustus_augustus_hints from './hints'
include extract_augustus_split_hints as extract_pasa_augustus_hints from './hints'
include extract_augustus_split_hints as extract_gemoma_augustus_hints from './hints'
include extract_augustus_split_hints as extract_gemoma_comparative_augustus_hints from './hints'
include extract_exonerate_hints from './hints'
include extract_exonerate_evm_hints from './hints'

include extract_spaln_transcript_evm_hints from './hints'
include extract_spaln_protein_evm_hints from './hints'
include extract_gmap_evm_hints from './hints'

include assert_same_names from './cli'
include assert_two_covers_one from './cli'

include busco_proteins from './evaluation'
include get_stats from './evaluation'
include get_splice_site_info from './evaluation'
include get_known_stats from './evaluation'

def is_null = { f -> (f == null || f == '') }


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
    // Require that we have the same indices for genomes.
    // This is kind of hacky since we don't take `genomes` as a parameter.
    assert_same_names(
        spaln_indices.map { it[0] },
        gmap_indices.map { it[0] },
        "spaln index",
        "gmap index",
        "We really want to make sure you get all of the results you'd expect :)"
    )

    (combined_fasta, combined_tsv) = combine_fastas(transcripts.collect())
    cleaned_transcripts = clean_transcripts(combined_fasta, univec)

    spaln_aligned = spaln_align_transcripts(
        species,
        max_gene_hard,
        min_intron_soft,
        spaln_indices.combine(cleaned_transcripts.map { f, cln, clean -> clean })
    )

    spaln_aligned_tidied = tidy_spaln_transcripts_gff3(
        "spaln_transcripts",
        "spaln",
        spaln_aligned
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
        spaln_aligned_tidied
    )

    spaln_evm_hints = extract_spaln_transcript_evm_hints(spaln_aligned_tidied)

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
    spaln_aligned_tidied
    spaln_augustus_hints
    spaln_evm_hints
    gmap_aligned
    gmap_evm_hints
}


/**
 * Align proteins to the genomes using spaln.
 */
workflow align_proteins {

    get:
    trans_table
    min_intron_soft
    max_gene_hard
    proteins
    spaln_indices

    main:
    (combined_fasta, combined_tsv) = combine_fastas(proteins.collect())

    spaln_aligned = spaln_align_proteins(
        trans_table,
        min_intron_soft,
        max_gene_hard,
        spaln_indices.combine(combined_fasta)
    )

    spaln_aligned_fixed = fix_spaln_proteins_stop(spaln_aligned)
    spaln_aligned_tidied = tidy_spaln_proteins_gff3(
        "spaln_proteins",
        "spaln",
        spaln_aligned_fixed
    )

    spaln_augustus_hints = extract_spaln_protein_augustus_hints(
        "spaln_proteins",
        "spaln",
        "P",
        3,
        9,
        9,
        9,
        9,
        false,
        spaln_aligned_tidied
    )

    spaln_evm_hints = extract_spaln_protein_evm_hints(spaln_aligned_tidied)

    emit:
    combined_fasta
    combined_tsv
    spaln_aligned_tidied
    spaln_augustus_hints
    spaln_evm_hints
}


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
    gemoma_rnaseq_hints_indiv = extract_gemoma_rnaseq_hints(
        genomes
            .combine(aligned_with_strand, by: 0)
            .map { n, f, rg, c, s -> [n, rg, f, c, s] }
    )

    gemoma_rnaseq_hints = combine_gemoma_rnaseq_hints(
        gemoma_rnaseq_hints_indiv.groupTuple(by: 0)
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
    // Require that all genomes have a cram.
    // This isn't perfect, because we don't check that each genome
    // has all read_groups
    assert_two_covers_one(
        genomes.map { n, f -> n },
        crams.map { n, rg, c, s -> n },
        "genomes",
        "crams",
        "We really want to be able to assemble stringtie transcripts for all genomes!",
    )

    known_with_null = genomes
        .join(known, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('WAS_NULL')]: [n, g] }

    individually_assembled = stringtie_assemble(
        genomes
            .combine(crams, by: 0)
            .map { n, f, rg, c, s -> [n, rg, f, c, s]}
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
    // Require that we have faidx for each genome
    assert_two_covers_one(
        genomes.map { n, f -> n },
        genome_faidxs.map { n, f -> n },
        "genome",
        "fasta faidx",
        "We really need faidx files for each genome."
    )

    // Require that we have mmseqs_index for each genome
    assert_two_covers_one(
        genomes.map { n, f -> n },
        genome_mmseqs_indices.map { n, f -> n },
        "genome",
        "mmseqs index",
        "We really need mmseqs indices for each genome."
    )

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
    // Require that we have gmap alignments for all genomes.
    assert_two_covers_one(
        genomes.map { n, f -> n },
        gmap.map { n, g -> n },
        "genomes",
        "gmap alignments",
        "We really want to make sure you get all of the results you'd expect :)"
    )

    known_with_null = genomes
        .join(known, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('KNOWN_WAS_NULL')]: [n, g] }

    stringtie_with_null = genomes
        .join(stringtie, by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, file('STRINGTIE_WAS_NULL')]: [n, g] }

    pasa_gff3 = pasa(
        not_fungus,
        max_intron_hard,
        genomes
            .join(known_with_null, by: 0)
            .join(stringtie_with_null, by: 0)
            .join(gmap, by: 0)
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
    assert_two_covers_one(
        genomes.map { n, f -> n },
        stringtie.map { n, g -> n },
        "genomes",
        "stringtie assemblies",
        "We require assemblies for each genome."
    )

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
        "codingquarrypm",
        "codingquarry",
        cq_fixed_gff3
    )

    cq_augustus_hints = extract_codingquarry_augustus_hints(
        "codingquarrypm",
        "codingquarrypm",
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
 * Runs the regular gemoma pipeline
 *
 * @param trans_table The integer corresponding to the NCBI translation table.
 * @param genomes A channel of genomes to run. Structure: tuple(val(name), path("in.fasta"))
 * @param genome_mmseqs_indices A channel of MMSeqs formatted genome databases.
 * @param introns A channel of intron hints to use for gemoma.
 *                Structure: tuple(val(name), path("hints.gff3"))
 * @param known known loci in gff3 format.
 */
workflow run_gemoma {

    get:
    trans_table
    genomes
    genome_mmseqs_indices
    introns
    known

    main:
    cds_parts = extract_gemoma_cds_parts(
        genomes.join(known, by: 0, remainder: false)
    )

    cds_matches = mmseqs_search_gemoma_cds_parts(
        trans_table,
        cds_parts
            .combine(genome_mmseqs_indices)
            .filter { rn, c, a, p, tn, g -> rn != tn }
    )

    gemoma_indiv_preds = gemoma(
        genomes
            .combine(cds_matches, by: 0)
            .map { tn, f, rn, m -> [rn, tn, f, m] }
            .combine(cds_parts, by: 0)
            .map { rn, tn, f, m, c, a, p -> [tn, rn, f, c, a, p, m] }
            .combine(introns, by: 0)
    )

    gemoma_gff3 = gemoma_combine(
        "gemoma",
        gemoma_indiv_preds
            .groupTuple(by: 0)
            .join(genomes, by: 0, remainder: false)
            .combine(introns, by: 0)
    )

    gemoma_gff3_tidied = tidy_gemoma_gff3(
        "gemoma",
        "gemoma",
        gemoma_gff3
    )

    gemoma_augustus_hints = extract_gemoma_augustus_hints(
        "gemoma",
        "gemoma_exon",
        "gemoma_cds",
        "GEMOMA",
        "GEMOMA",
        4, // exon_priority
        3, // CDS priority
        9, // exon_trim
        6, // cds_trim
        9, // utr_trim,
        6, // gene_trim
        false,
        gemoma_gff3_tidied
    )

    emit:
    gemoma_gff3
    gemoma_gff3_tidied
    gemoma_augustus_hints
}


/**
 * @param d
 */
workflow get_augustus_hints {

    get:
    min_intron_hard
    max_intron_hard
    known
    spaln_transcripts
    spaln_proteins
    exonerate  // tuple(val(name), val(protein_name), path("evm.gff3"))
    genemark
    pasa
    codingquarry
    codingquarrypm
    gemoma
    augustus
    gemoma_comparative

    main:
    spaln_transcripts_augustus_hints = extract_spaln_transcript_augustus_hints(
        "spaln_transcripts",
        "spaln",
        "E",
        3,
        6,
        0,
        0,
        0,
        false,
        spaln_transcripts
    )

    spaln_proteins_augustus_hints = extract_spaln_protein_augustus_hints(
        "spaln_proteins",
        "spaln",
        "P",
        3,
        9,
        9,
        9,
        9,
        false,
        spaln_proteins
    )

    exonerate_augustus_hints = extract_exonerate_hints(
        min_intron_hard,
        max_intron_hard,
        exonerate
    )

    genemark_augustus_hints = extract_genemark_augustus_hints(
        "genemark",
        "genemark",
        "GM",
        3, // priority
        9, // exon_trim
        9, // cds_trim
        9, // utr trim
        9, // gene trim
        false,
        genemark
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
        pasa
    )

    cq_augustus_hints = extract_codingquarry_augustus_hints(
        "codingquarry",
        "codingquarry",
        "CQ",
        4, // priority
        6, // exon_trim
        6, // cds_trim
        6, // utr trim
        6, // gene trim
        false,
        codingquarry
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
        codingquarrypm
    )

    gemoma_augustus_hints = extract_gemoma_augustus_hints(
        "gemoma",
        "gemoma_exon",
        "gemoma_cds",
        "GEMOMA",
        "GEMOMA",
        4, // exon_priority
        3, // CDS priority
        9, // exon_trim
        6, // cds_trim
        9, // utr_trim,
        6, // gene_trim
        false,
        gemoma
    )

    augustus_augustus_hints = extract_augustus_augustus_hints(
        "augustus",
        "augustus",
        "AUG",
        4, // priority
        6, // exon_trim
        6, // cds_trim
        6, // utr_trim,
        6, // gene_trim
        false,
        augustus
    )

    gemoma_comparative_augustus_hints = extract_gemoma_comparative_augustus_hints(
        "gemoma_comparative",
        "gemoma_comparative_exon",
        "gemoma_comparative_cds",
        "COMPGEMOMA",
        "COMPGEMOMA",
        4, // exon_priority
        3, // CDS priority
        9, // exon_trim
        6, // cds_trim
        9, // utr_trim,
        6, // gene_trim
        false,
        gemoma_comparative
    )

    emit:
    spaln_transcripts_augustus_hints
    spaln_proteins_augustus_hints
    exonerate_augustus_hints
    genemark_augustus_hints
    pasa_augustus_hints
    cq_augustus_hints
    cqpm_augustus_hints
    gemoma_augustus_hints
    augustus_augustus_hints
    gemoma_comparative_augustus_hints
}


/**
 * @param species The augustus species to use. Should be available in the config_dir
 * @param predict_utrs Boolean, whether augustus should predict utrs.
 * @param not_fungus Boolean, should we not use fungal specific parameters?
 * @param min_intron_hard the minimum length an intron should be.
 * @param valid_splicesites a comma separated list of valid splice sites e.g. "atac,gtag"
 * @param genomes A channel containing the genomes to predict, tuple(val(name), path("in.fasta"))
 * @param known A channel containing gffs with known sites. tuple(val(name), path("known.gff3"))
 * @param config_dir A value channel containing the path to the augustus config directory.
 * @param extrinsic the extrinsic config file, containing weights for hints.
 * @param hints A channel of hints files for augustus. tuple(val(name), path("hint.gff3"))
 *
 * @return augustus_gff3_tidied The tidied augustue predictions.
 * @return augustus_augustus_hints Hints to use for augustus during gapfiller stage.
 */
workflow run_augustus {

    get:
    species
    predict_utrs
    not_fungus
    min_intron_hard
    valid_splicesites
    genomes
    known
    config_dir
    extrinsic
    hints

    main:
    genome_chunks = chunkify_genomes(16, genomes)

    if ( predict_utrs && !not_fungus ) {
        genome_chunks_with_strand = genome_chunks
            .flatMap { n, fs -> fs.collect { f -> [n, f] } }
            .flatMap { n, f -> [[n, "forward", f], [n, "reverse", f]] }
    } else {
        genome_chunks_with_strand = genome_chunks
            .flatMap { n, fs -> fs.collect { f -> [n, "both", f] } }
    }

    known_gff3_tidied = tidy_known_gff3(
        "known",
        "known",
        known
    )

    known_augustus_hints = extract_known_augustus_hints(
        "known",
        "known",
        "M",
        50, // priority
        0, // exon_trim
        0, // cds_trim
        0, // utr_trim,
        0, // gene_trim
        false,
        known_gff3_tidied
    )

    augustus_gff3_chunks = augustus_hints(
        species,
        predict_utrs,
        not_fungus,
        min_intron_hard,
        valid_splicesites,
        genome_chunks_with_strand.combine(
            hints.mix(known_augustus_hints.map { n, a, g -> [n, g] }).groupTuple(by: 0),
            by: 0
        ),
        config_dir,
        extrinsic
    )

    augustus_gff3_tidied = combine_and_tidy_augustus_gff3(
        "augustus",
        "augustus",
        augustus_gff3_chunks.map { n, s, g -> [n, g] }.groupTuple(by: 0)
    )

    augustus_augustus_hints = extract_augustus_augustus_hints(
        "augustus",
        "augustus",
        "AUG",
        4, // priority
        6, // exon_trim
        6, // cds_trim
        6, // utr_trim,
        6, // gene_trim
        false,
        augustus_gff3_tidied
    )

    emit:
    augustus_gff3_tidied
    augustus_augustus_hints
}


/**
 * Runs the comparative gemoma pipeline.
 *
 * @param trans_table The integer corresponding to the NCBI translation table.
 * @param genomes A channel of genomes to run. Structure: tuple(val(name), path("in.fasta"))
 * @param genome_mmseqs_indices A channel of MMSeqs formatted genome databases.
 * @param introns A channel of intron hints to use for gemoma.
 *                Structure: tuple(val(name), path("hints.gff3"))
 * @param predictions A channel of predictions to transfer.
 *                Structure tuple(val(name), val(analysis), path("in.gff3")
 * @return gemoma_predictions The raw genoma predictions.
 *         Structure tuple(val(name), path(gff3)).
 * @return gemoma_predictions_tidied The tidied gemoma predictions.
 *         Structure tuple(val(name), path(gff3)).
 */
workflow run_gemoma_comparative {

    get:
    trans_table  // val, should be ncbi table integer.
    genomes // tuple(val(name), file(fasta))
    genome_mmseqs_indices
    introns  // tuple(val(name), file(introns.gff))
    predictions // tuple(val(name), val(analysis), file(gff3))

    main:

    cds_parts = extract_gemoma_comparative_cds_parts(
        genomes
            .combine(predictions, by: 0)
            .map { n, f, a, g -> [n, a, f, g] }
    )

    clustered_cds_parts = cluster_gemoma_cds_parts(
        cds_parts
            .map {n, an, c, a, p -> [c, a, p]}
            .toList()
            .map { l -> l.transpose() }
            .collect()
    )


    mmseqs_matches = mmseqs_search_gemoma_comparative_cds_parts(
        trans_table,
        clustered_cds_parts
            .map { c, a, p -> ["comparative", c, a, p] }
            .combine(genome_mmseqs_indices)
    )

    gemoma_matches = gemoma_comparative(
        genomes.combine(
            clustered_cds_parts.map {c, a, p -> ["comparative", c, a, p]},
        )
        .map { t, f, r, c, a, p -> [t, r, f, c, a, p] }
        .join(mmseqs_matches.view(), by: [0, 1])
        .join(introns, by: 0)
    )

    gemoma_gff3 = gemoma_comparative_combine(
        "gemoma_comparative",
        gemoma_matches
            .groupTuple(by: 0)
            .join(genomes, by: 0)
            .join(introns, by: 0)
    )

    gemoma_gff3_tidied = tidy_gemoma_comparative_gff3(
        "gemoma_comparative",
        "gemoma_comparative",
        gemoma_gff3
    )

    gemoma_augustus_hints = extract_gemoma_comparative_augustus_hints(
        "gemoma_comparative",
        "gemoma_comparative_exon",
        "gemoma_comparative_cds",
        "COMPGEMOMA",
        "COMPGEMOMA",
        4, // exon_priority
        3, // CDS priority
        9, // exon_trim
        6, // cds_trim
        9, // utr_trim,
        6, // gene_trim
        false,
        gemoma_gff3_tidied
    )

    emit:
    gemoma_gff3
    gemoma_gff3_tidied
    gemoma_augustus_hints
}


workflow get_evm_hints {

    get:
    min_intron_hard
    max_intron_hard
    spaln_transcripts
    gmap
    spaln_proteins
    exonerate // tuple(val(name), val(protein_name), path(hints.gff3))

    main:
    spaln_transcript_evm_hints = extract_spaln_transcript_evm_hints(spaln_transcripts)
    gmap_evm_hints = extract_gmap_evm_hints(gmap)
    spaln_protein_evm_hints = extract_spaln_protein_evm_hints(spaln_proteins)

    exonerate_evm_hints = extract_exonerate_evm_hints(
        min_intron_hard,
        max_intron_hard,
        exonerate
    )

    emit:
    spaln_transcript_evm_hints
    gmap_evm_hints
    spaln_protein_evm_hints
    exonerate_evm_hints
}


workflow run_evm {

    get:
    species
    not_fungus
    augustus_utr
    valid_splicesites
    min_intron_hard
    augustus_config
    augustus_hint_weights
    evm_weights
    genomes
    genome_faidxs
    transcripts
    proteins
    genes
    known
    hints
    evm_done
    augustus_gapfiller_done

    main:
    // Config expects source to be "manual"
    known_gff3_tidied = tidy_known_gff3(
        "manual",
        "manual",
        known
    )

    transcripts_with_null = genomes
        .join(transcripts.groupTuple(by: 0), by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, []]: [n, g] }

    proteins_with_null = genomes
        .join(proteins.groupTuple(by: 0), by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, []]: [n, g] }

    genes_with_null = genomes
        .join(genes.mix(known_gff3_tidied).groupTuple(by: 0), by: 0, remainder: true)
        .map { n, f, g -> is_null(g) ? [n, []]: [n, g] }


    evm_gff3 = evm(
        min_intron_hard,
        evm_weights,
        genomes
            .join(evm_done, by: 0, remainder: true)
            .filter { n, f, d -> is_null(d) }
            .map { n, f, d -> [n, f] }
            .join(transcripts_with_null, by: 0)
            .join(proteins_with_null, by: 0)
            .join(genes_with_null, by: 0)
    )

    evm_gff3_tidied = tidy_evm_gff3(
        "evm",
        "evm",
        evm_gff3
    )

    missing_preds = find_missing_evm_predictions(
        evm_gff3_tidied.mix(evm_done)
            .join(genes.groupTuple(by: 0), by: 0)
            .join(genome_faidxs, by: 0)
    )

    known_augustus_hints = extract_known_augustus_hints(
        "known",
        "known",
        "M",
        50, // priority
        0, // exon_trim
        0, // cds_trim
        0, // utr_trim,
        0, // gene_trim
        false,
        known_gff3_tidied
    )

    augustus_gff3s = augustus_gap_filler(
        species,
        not_fungus,
        augustus_utr,
        valid_splicesites,
        min_intron_hard,
        genomes
            .join(augustus_gapfiller_done, by: 0, remainder: true)
            .filter { n, f, d -> is_null(d) }
            .map { n, f, d -> [n, f] }
            .join(missing_preds, by: 0)
            .join(
                hints
                    .mix(known_augustus_hints.map { n, a, g -> [n, g] })
                    .groupTuple(by: 0),
                by: 0
            ),
        augustus_config,
        augustus_hint_weights
    )

    augustus_gff3_tidied = combine_and_tidy_augustus_gff3(
        "augustus_gapfiller",
        "augustus",
        augustus_gff3s
    )

    final_gff3 = merge_gffs(
        "final",
        evm_gff3_tidied
            .mix(evm_done, augustus_gff3_tidied, augustus_gapfiller_done)
            .groupTuple(by: 0)
    )

    emit:
    evm_gff3_tidied
    augustus_gff3_tidied
    final_gff3
}


workflow filter_preds {

    get:
    predictions
    spaln_transcripts
    spaln_proteins
    gmap_transcripts
    exonerate
    genemark
    pasa
    codingquarry
    codingquarrypm
    gemoma
    augustus
    gemoma_comparative

    main:
    spaln_transcripts_beds = spaln_transcripts_gff_to_bed(
        "Target",
        "exon",
        "spaln_transcript",
        false,
        spaln_transcripts
    )

    spaln_proteins_beds = spaln_proteins_gff_to_bed(
        "Target",
        "CDS",
        "spaln_protein",
        false,
        spaln_proteins
    )

    gmap_transcripts_beds = gmap_transcripts_gff_to_bed(
        "Name",
        "cDNA_match",
        "gmap_transcript",
        false,
        gmap_transcripts
    )

    exonerate_gff3s = exonerate_to_gff3(exonerate)
    exonerate_beds = exonerate_gff_to_bed(
        "query",
        "CDS",
        "exonerate",
        false,
        exonerate_gff3s
    )

    genemark_beds = genemark_gff_to_bed(
        "Parent",
        "CDS",
        "genemark",
        false,
        genemark
    )

    pasa_beds = pasa_gff_to_bed(
        "Parent",
        "CDS",
        "pasa",
        false,
        pasa
    )

    cq_beds = cq_gff_to_bed(
        "Parent",
        "CDS",
        "codingquarry",
        false,
        codingquarry
    )

    cqpm_beds = cqpm_gff_to_bed(
        "Parent",
        "CDS",
        "codingquarrypm",
        false,
        codingquarrypm
    )

    gemoma_beds = gemoma_gff_to_bed(
        "Parent",
        "CDS",
        "gemoma",
        false,
        gemoma
    )

    augustus_beds = augustus_gff_to_bed(
        "Parent",
        "CDS",
        "augustus",
        false,
        augustus
    )

    gemoma_comparative_beds = gemoma_comparative_gff_to_bed(
        "Parent",
        "CDS",
        "gemoma_comparative",
        false,
        gemoma_comparative
    )

    cov = get_hint_coverage(
        "CDS",
        predictions.join(
            spaln_transcripts_beds.mix(
                spaln_proteins_beds,
                gmap_transcripts_beds,
                exonerate_beds,
                genemark_beds,
                pasa_beds,
                cq_beds,
                cqpm_beds,
                gemoma_beds,
                augustus_beds,
                gemoma_comparative_beds,
            ).groupTuple(by: 0),
            by: 0
        )
    )

    emit:
    cov
}


workflow run_stats {

    get:
    trans_table
    genomes
    gffs // tuple(val(name), val(analysis), path(gff))
    known
    busco_lineage
    augustus_config

    main:
    (proteins, cdss) = extract_seqs(trans_table, gffs.combine(genomes, by: 0))

    if (busco_lineage) {
        busco_preds = busco_proteins(proteins, busco_lineage, augustus_config)
    } else {
        busco_preds = Channel.empty()
    }

    stats = get_stats(gffs)

    splice_site_stats = get_splice_site_info(gffs.combine(genomes, by: 0))
    known_stats = get_known_stats(gffs.join(known, by: 0))

    emit:
    proteins
    cdss
    busco_preds
    stats
    splice_site_stats
    known_stats
}
