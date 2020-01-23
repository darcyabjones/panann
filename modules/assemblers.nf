#!/usr/bin/env nextflow


/*
 * Assemble transcripts from RNAseq with Stringtie
 * Note that using the known sites option just completely ignores
 * rnaseq info at the loci so it may appear that No UTRs are included
 * if they aren't in the known sites file.
 */
process stringtie_assemble {

    label "stringtie"
    label "medium_task"
    time '6h'

    tag "${name} - ${read_group}"

    input:
    tuple val(name),
        val(read_group),
        path(fasta),
        path(faidx),
        path(cram),
        val(strand),
        path(gff)

    output:
    tuple val(name),
        val(read_group),
        path("${name}_${read_group}_stringtie.gtf")

    script:
    def strand_flag = strand == "fr" ? "--fr" : "--rf"
    def known = gff.name != 'WAS_NULL' ? "-G ${gff}" : ''

    """
    # Convert cram to bam.
    samtools view \
        -b \
        -T "${fasta}" \
        -@ "${task.cpus}" \
        -o "tmp.bam" \
        "${cram}"

    stringtie \
      -p "${task.cpus}" \
      ${strand_flag} \
      ${known} \
      -o "${name}_${read_group}_stringtie.gtf" \
      -m 150 \
      "tmp.bam"

    rm -f tmp.bam
    """
}


/*
 * Combine stringtie annotations from multiple bams.
 */
process stringtie_merge {

    label "stringtie"
    label "medium_task"
    time '6h'

    tag "${name}"

    input:
    tuple val(name), path("*gtf"), path(gff)

    output:
    tuple val(name), path("${name}_stringtie.gtf")

    script:
    def known = gff.name != 'WAS_NULL' ? "-G ${gff}" : ''

    """
    stringtie \
      -p "${task.cpus}" \
      ${known} \
      --merge \
      -o "${name}_stringtie.gtf" \
      *gtf
    """
}


/*
 * Assemble reads into transcripts with trinity.
 */
process trinity_assemble_denovo {

    label "trinity"
    label "big_task"
    time '1d'

    tag "${read_group}"

    input:
    val not_fungus
    tuple val(read_group),
        path(r1s),
        path(r2s),
        val(strand)

    output:
    tuple val(read_group), path("${read_group}_trinity_denovo.fasta")

    script:
    def r1_joined = r1s.join(',')
    def r2_joined = r2s.join(',')
    def use_jaccard = not_fungus ? '' : "--jaccard_clip "
    def strand_flag = strand == "fr" ? "--SS_lib_type FR " : "--SS_lib_type RF "

    """
    Trinity \
      --seqType fq \
      --max_memory "${task.memory.toGiga()}G" \
      --CPU "${task.cpus}" \
      ${use_jaccard} \
      ${strand_flag} \
      --output trinity_assembly \
      --left "${r1_joined}" \
      --right "${r2_joined}"

    mv trinity_assembly/Trinity.fasta "${read_group}_trinity_denovo.fasta"
    rm -rf -- trinity_assembly
    """
}
