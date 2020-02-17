#!/usr/bin/env nextflow


/*
 * Evaluate genome completeness with BUSCO on the genomes.
 * Later we evaluate each gene prediction set too.
 * Could compare this number with that one.
 */
process busco {

    label "busco"
    label "medium_task"
    time '12h'

    tag "${name}"

    input:
    tuple val(name), file(fasta), file(faidx) from genomes4Busco
    path "lineage" from buscoLineage
    path "augustus_config" from augustusConfig

    output:
    path "${name}" into buscoResults

    script:
    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    run_BUSCO.py \
      --in "${fasta}" \
      --out "${name}" \
      --cpu ${task.cpus} \
      --mode "genome" \
      --lineage_path "lineage"

    mv "run_${name}" "${name}"
    """
}


/*
 * Evaluate gene predictions using protein comparisons with BUSCO sets.
 */
process busco_proteins {

    label "busco"
    label "medium_task"
    time '6h'

    tag "${name} - ${analysis}"

    input:
    tuple val(name),
          val(analysis),
          path(fasta)
    path "lineage"
    path "augustus_config"

    output:
    path "${analysis}_busco"

    script:
    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    run_BUSCO.py \
      --in "${fasta}" \
      --out "${name}" \
      --cpu ${task.cpus} \
      --mode "proteins" \
      --lineage_path "lineage"

    mv "run_${name}" "${analysis}_busco"
    """
}


/*
 * Get number of genes, exons, distributions of lengths etc.
 */
process get_stats {

    label "genometools"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    input:
    tuple val(name), val(analysis), path(preds)

    output:
    path "${name}_${analysis}_stats.txt"

    script:
    """
    gt stat \
      -addintrons \
      -genelengthdistri \
      -genescoredistri \
      -exonlengthdistri \
      -exonnumberdistri \
      -intronlengthdistri \
      -cdslengthdistri \
      "${preds}" \
    > "${name}_${analysis}_stats.txt"
    """
}


/*
 * Get distributions of splice site pairs.
 */
process get_splice_site_info {

    label "genometools"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    input:
    path val(name),
        val(analysis),
        path(preds),
        path(fasta)

    output:
    path "${name}_${analysis}_splice_sites.txt"

    script:
    """
    gt splicesiteinfo \
      -seqfile "${fasta}" \
      -matchdescstart \
      -addintrons \
      "${preds}" \
    > "${name}_${analysis}_splice_sites.txt"
    """
}


/*
 * Get stats when we have known sites
 */
process get_known_stats {

    label "aegean"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    input:
    tuple val(name),
          val(analysis),
          path(preds),
          path(known)

    output:
    path "${name}_${analysis}_parseval.txt"

    script:
    """
    parseval \
      --nogff3 \
      --outformat "text" \
      --summary \
      --outfile "${name}_${analysis}_parseval.txt" \
      "${known}" \
      "${preds}"
    """
}
