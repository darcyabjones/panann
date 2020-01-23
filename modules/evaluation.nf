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

    publishDir "${params.outdir}/qc/${name}"

    when:
    params.busco_lineage && params.busco_genomes

    input:
    set val(name), file(fasta), file(faidx) from genomes4Busco
    file "lineage" from buscoLineage
    file "augustus_config" from augustusConfig

    output:
    file "${name}" into buscoResults

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
process runBuscoProteins {

    label "busco"
    label "medium_task"
    time '6h'

    tag "${name} - ${analysis}"

    publishDir "${params.outdir}/qc/${name}"

    when:
    params.busco_lineage

    input:
    set val(name), val(analysis), file(fasta) from extractedProteins
    file "lineage" from buscoLineage
    file "augustus_config" from augustusConfig

    output:
    file "${analysis}_busco" into buscoProteinsResults

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
process getStats {

    label "genometools"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    publishDir "${params.outdir}/qc/${name}"

    input:
    set val(name), val(analysis), file(preds) from predictions4Stats

    output:
    file "${name}_${analysis}_stats.txt"

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
process getSpliceSiteInfo {

    label "genometools"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    publishDir "${params.outdir}/qc/${name}"

    input:
    set val(name),
        val(analysis),
        file(preds),
        file(fasta),
        file(faidx) from predictions4GetSpliceSiteInfo
            .combine(genomes4GetSpliceSiteInfo, by: 0)

    output:
    file "${name}_${analysis}_splice_sites.txt"

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
process getKnownStats {

    label "aegean"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    publishDir "${params.outdir}/qc/${name}"

    input:
    set val(name), val(analysis), file(preds), file(known) from predictions4KnownStats
        .combine(genomes4GetKnownStats.map { n, f, i, g -> [n, g] }, by: 0)
        .filter { n, a, p, k -> k.name != "WAS_NULL" }

    output:
    file "${name}_${analysis}_parseval.txt"

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
