#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # panann

    A pipeline to predict genes in a multiple fungal genomes.

    ## Usage

    ```bash
    cat <<EOF > known_sites.tsv
    genome	gff	confidence
    genome1.fasta	genome1.gff3	8
    EOF

    nextflow run main.nf \
      --genomes "genomes/*.fasta" \
      --known_sites known_sites.tsv \
      --transcripts transcripts.fasta \
      --proteins proteins.fasta \
      --remote_proteins fungal_proteins.fasta
    ```

    ## Parameters
    genomes - fasta genomes
    known_sites - table
    transcripts - already assembled transcripts (skip assembly?)
    proteins - predicted proteins from closely related species/isolate

    remote_proteins - protein sequences from remotely related isolates, to be weighted differently.
    genome_alignment - A multiple genome alignment in MAF format, e.g. from cactus or sibelliaz.
      if not provided, align the genomes with sibelliaz

    busco_lineage - The busco lineage directory (un-tarred).
    augustus_config - The augustus config directory. If not provided, assumes
      that AUGUSTUS_CONFIG_PATH is set where the tasks will be executed.

    fastq - table, rnaseq fastq files to assemble. need name, r1, r2. If add fasta field, do reference guided assembly.
    bams - table, already aligned bam-files. Need name/fasta and bam.
    cufflinks_gtf - Cufflinks assemblies, need name/fasta and gtf.


    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


params.genomes = false

params.known_sites = false

params.transcripts = false

params.proteins = false
params.remote_proteins = false

params.genome_alignment = false

params.busco_lineage = false

params.augustus_config = false

// RNAseq params
params.fastq = false
params.bams = false
params.cufflinks_gtf = false

/*
 * Sanitise the input
 */

if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.baseName, g] }
        .set { genomes }
} else {
    log.error "Please provide some genomes to predict genes for with `--genomes`."
    exit 1
}


if ( params.transcripts ) {
    Channel
        .fromPath(params.transcripts, checkIfExists: true, type: "file")
        .collectFile(name: "transcripts.fasta", newLine: true, sort: "deep")
        .first()
        .set { transcripts }
} else {
    transcripts = false
}

if ( params.proteins ) {
    Channel
        .fromPath(params.proteins, checkIfExists: true, type: "file")
        .collectFile(name: "proteins.fasta", newLine: true, sort: "deep")
        .first()
        .set { proteins }
} else {
    proteins = false
}

if ( params.remote_proteins ) {
    Channel
        .fromPath(params.remote_proteins, checkIfExists: true, type: "file")
        .collectFile(name: "remote_proteins.fasta", newLine: true, sort: "deep")
        .first()
        .set { remoteProteins }
} else {
    remoteProteins = false
}

if ( params.busco_lineage ) {
    Channel
        .fromPath(params.busco_lineage, checkIfExists: true, type: "file")
        .first()
        .set { buscoLineage }
} else {
    buscoLineage = false
}

if ( params.genome_alignment ) {
    Channel
        .fromPath(params.genome_alignment, checkIfExists: true, type: "file")
        .first()
        .set { genomeAlignment }
} else {
    // We can align it here or can do it later.
    genomeAlignment = false
}

// Because augustus requires the config folder to be editable, it's easier
// to keep track of it in the pipeline.
// Assumes the correct environment variable is set.
if ( params.augustus_config ) {
    Channel
        .fromPath(params.augustus_config, checkIfExists: true, type: "dir")
        .first()
        .set { augustusConfig }
} else {
    process getAugustusConfig {

        label "augustus"
        label "small_task"

        output:
        file "config" into augustusConfig

        """
        cp -r \${AUGUSTUS_CONFIG_PATH} ./config
        """
    }
}


if ( params.fastq ) {
    Channel
        .fromPath(params.fastq, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { it.name != null && it.read1_file != null && it.read2_file != null }
        .map {[
            it.name,
            (it.genome != null) ? file(it.genome).baseName : null,
            file(it.read1_file, checkIfExists: true),
            file(it.read2_file, checkIfExists: true)
        ]}
        .set { fastq }
} else {
    fastq = Channel.empty()
}


// Split up the inputs.

genomes.into {
    genomes4Busco;
    genomes4SpalnIndex;
    genomes4GmapIndex;
}

fastq.into {
    fastq4TrinityAssemble;
    fastq4StarAlign;
}


// Indexing and preprocessing

process getSpalnIndex {

    label "spaln"
    label "small_task"

    tag { name }

    input:
    set val(name), file(genome) from genomes4SpalnIndex

    output:
    set val(name),
        file("${genome.baseName}.bkn"),
        file("${genome.baseName}.ent"),
        file("${genome.baseName}.idx"),
        file("${genome.baseName}.bkp"),
        file("${genome.baseName}.grp"),
        file("${genome.baseName}.seq") into spalnIndices

    """
    makeidx.pl -inp ${genome}
    """
}

spalnIndices.into {
    spalnIndices4AlignSpalnTranscripts;
    spalnIndices4AlignSpalnProteins;
    spalnIndices4AlignSpalnRemoteProteins;
}

process getGmapIndex {

    label "gmap"
    label "medium_task"

    tag { name }

    input:
    set val(name), file(genome) from genomes4GmapIndex

    output:
    set val(name), file("db") into gmapIndices

    """
    gmap_build \
      -k 13 \
      -D db \
      -d "${name}" \
      ${genome}
    """
}


// RNAseq alignments and assembly

process trinityAssembleTranscripts {

}

process starAlignTranscripts {

}

// 1 align transcripts to all genomes
process alignSpalnTranscripts {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/spaln"

    tag { name }

    when:
    params.transcripts

    input:
    file transcripts

    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices4AlignSpalnTranscripts

    output:
    set val(name), file("${name}_transcripts.gff3") into spalnAlignedTranscripts

    """
    spaln \
      -L \
      -M3 \
      -O3 \
      -Q7 \
      -S3 \
      -TDothideo \
      -yX \
      -yS \
      -ya2 \
      -t ${task.cpus} \
      -d "${name}" \
      "${transcripts}" \
    > "${name}_transcripts.gff3"
    """
}

process alignGmapTranscripts {

    label "gmap"
    label "medium_task"

    publishDir "${params.outdir}/aligned/gmap"

    tag { name }

    when:
    params.transcripts

    input:
    file transcripts
    set val(name), file("db") from gmapIndices

    output:
    set val(name), file("${name}_transcripts.gff3") into gmapAlignedTranscripts

    script:
    // Options to parametrise
    //-n, --npaths=INT               Maximum number of paths to show (default 5).  If set to 1, GMAP
    //                                 will not report chimeric alignments, since those imply
    //                                 two paths.  If you want a single alignment plus chimeric
    //                                 alignments, then set this to be 0.
    //--min-intronlength=INT         Min length for one internal intron (default 9).  Below this size,
    //                                 a genomic gap will be considered a deletion rather than an intron.
    //--max-intronlength-middle=INT  Max length for one internal intron (default 500000).  Note: for backward
    //                                 compatibility, the -K or --intronlength flag will set both
    //                                 --max-intronlength-middle and --max-intronlength-ends.
    //                                 Also see --split-large-introns below.
    //--max-intronlength-ends=INT    Max length for first or last intron (default 10000).  Note: for backward
    //                                 compatibility, the -K or --intronlength flag will set both
    //                                 --max-intronlength-middle and --max-intronlength-ends.
    //--trim-end-exons=INT           Trim end exons with fewer than given number of matches
    //                                 (in nt, default 12)
    //--microexon-spliceprob=FLOAT   Allow microexons only if one of the splice site probabilities is
    //                                 greater than this value (default 0.95)
    //--canonical-mode=INT           Reward for canonical and semi-canonical introns
    //                                 0=low reward, 1=high reward (default), 2=low reward for
    //                                 high-identity sequences and high reward otherwise
    //--cross-species                Use a more sensitive search for canonical splicing, which helps especially
    //                                 for cross-species alignments and other difficult cases
    """
    gmap \
      --chimera-margin=50 \
      --npaths=1 \
      --format=gff3_match_cdna \
      --nthreads "${task.cpus}" \
      -D db \
      -d "${name}" \
      ${transcripts} \
    > ${name}_transcripts.gff3
    """
}

// Align with minimap2


// 2 align proteins to all genomes

process alignSpalnProteins {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/spaln"

    tag { name }

    when:
    params.proteins

    input:
    file proteins

    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices4AlignSpalnProteins

    output:
    set val(name), file("${name}_proteins.gff3")

    """
    spaln \
      -L \
      -M3 \
      -O3 \
      -Q7 \
      -TDothideo \
      -ya2 \
      -t ${task.cpus} \
      -a \
      -d "${name}" \
      "${proteins}" \
    > "${name}_proteins.gff3"
    """
}

// 3 align remote proteins to all genomes

process alignSpalnRemoteProteins {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/spaln"

    tag { name }

    when:
    params.remote_proteins

    input:
    file remoteProteins

    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices4AlignSpalnRemoteProteins

    output:
    set val(name), file("${name}_remote_proteins.gff3")

    script:
    """
    spaln \
      -L \
      -M3 \
      -O3 \
      -Q7 \
      -TDothideo \
      -ya1 \
      -t ${task.cpus} \
      -a \
      -d "${name}" \
      "${proteins}" \
    > "${name}_proteins.gff3"
    """
}

// 4 run busco on all genomes

process runBusco {
    label "busco"
    label "medium_task"
    tag { name }

    publishDir "${params.outdir}/busco"

    when:
    params.busco_lineage

    input:
    set val(name), file(fasta) from genomes4Busco
    file "lineage" from buscoLineage
    file "augustus_config" from augustusConfig

    output:
    file "${name}" into buscoResults

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

// 5 Run genemark, pasa, braker2, codingquarry on all genomes using previous steps.

// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats
