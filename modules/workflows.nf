include get_mmseqs_protein_db get_mmseqs_protein_db as get_mmseqs_remote_protein_db from './modules/aligners'
include mmseqs_search_genome_against_proteins as mmseqs_search_genome_against_remote_proteins from './modules/aligners'
include cluster_genome_vs_protein_matches as cluster_genome_vs_remote_protein_matches from './modules/aligners'
include exonerate_regions as exonerate_remote_proteins from './modules/aligners'

include fasta_to_tsv as remote_proteins_fasta_to_tsv from './modules/utils'


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
