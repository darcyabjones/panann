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
    known_sites - table. name, gff3. Note gff3 must have the exons set for star to work. Use genometools to add this.
    transcripts - already assembled transcripts (skip assembly?)
    proteins - predicted proteins from closely related species/isolate

    remote_proteins - protein sequences from remotely related isolates, to be weighted differently.
    genome_alignment - A multiple genome alignment in MAF format, e.g. from cactus or sibelliaz.
      if not provided, align the genomes with sibelliaz

    busco_lineage - The busco lineage directory (un-tarred).
    augustus_config - The augustus config directory. If not provided, assumes
      that AUGUSTUS_CONFIG_PATH is set where the tasks will be executed.

    fastq - table, rnaseq fastq files to assemble. need read_group, r1, r2. If add fasta field, do reference guided assembly.
    bams - table, already aligned bam-files. Need name/fasta and bam.
    cufflinks_gtf - Cufflinks assemblies, need name/fasta and gtf.
    fr -- Use if the stranded RNA seq is not RF.

    notfungus - dont use fungal specific features
    genemark - flag to use genemark (it's installed)
    signalp - flag to use signalp to run CodingQuarry - Pathogen mode.
    softmasked


    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.


    ## Output
    aligned_reads - Fastq reads aligned to each genome as bams.
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
params.fr = false
params.bams = false
params.cufflinks_gtf = false

params.notfungus = false
params.genemark = false
params.signalp = false

/*
 * Sanitise the input
 */

def is_null = { f -> (f == null || f == '') }

if ( params.genomes ) {
    Channel
        .fromPath(params.genomes, checkIfExists: true, type: "file")
        .map { g -> [g.baseName, g] }
        .set { genomes }
} else {
    log.error "Please provide some genomes to predict genes for with `--genomes`."
    exit 1
}

if ( params.known_sites ) {
    // Todo: merge with genomes here and error if get multiple gffs per genome
    // or names that don't match any genome
    Channel
        .fromPath(params.known_sites, checkIfExists: true, type: "file")
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { (!is_null(it.name) && !is_null(it.gff3)) }
        .map {[
            it.name,
            file(it.gff3, checkIfExists: true),
        ]}
        .unique()
        .set { knownSites }

} else {
    knownSites = Channel.empty()
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
        .filter { !is_null(it.read_group) && !is_null(it.r1) && !is_null(it.r2) }
        .map {[
            it.read_group,
            !is_null(it.name) ? it.name : 'WAS_NULL',
            file(it.r1, checkIfExists: true),
            file(it.r2, checkIfExists: true)
        ]}
        .unique()
        .set { fastq }

        // Split these into genome guided and non-genome guided.
} else {
    fastq = Channel.empty()
}


// Split up the inputs.

genomes.into {
    genomes4Busco;
    genomes4SpalnIndex;
    genomes4GmapIndex;
    genomes4StarIndex;
    genomes4GetFaidx;
    genomes4RunCodingQuarry;
    genomes4RunCodingQuarryPM;
}

if ( params.known_sites ) {
    genomes4StarIndex
        .join( knownSites4StarIndex, by: 0, remainder: true)
        .view()
        .filter { n, f, g -> (!is_null(f) && !is_null(n)) }
        .map { n, f, g -> [n, f, !is_null(g) ? g : file('WAS_NULL')] }
        .into {
            genomes4StarIndexWithGff;
            genomes4AssembleStringtie;
            genomes4MergeStringtie;
        }
} else {
    genomes4StarIndex
        .map { n, f -> [n, f, file('WAS_NULL')] }
        .into {
            genomes4StarIndexWithGff;
            genomes4AssembleStringtie;
            genomes4MergeStringtie;
        }
}

fastq
    .tap { fastqNoGenome; fastq4Alignment }
    .filter { rg, n, r1, r2 -> n != 'WAS_NULL' }
    .set { fastq4TrinityAssembleGuided }

fastqNoGenome
    .filter {rg, n, r1, r2 -> n == 'WAS_NULL' }
    .map { rg, n, r1, r2 -> [rg, r1, r2] }
    .unique()
    .set { fastq4TrinityAssemble }

fastq4Alignment
    .map { rg, n, r1, r2 -> [rg, r1, r2] }
    .unique()
    .into {
        fastq4StarFindNovelSpliceSites;
        fastq4StarAlignReads;
    }

knownSites.set { knownSites4StarIndex }

fastq4TrinityAssemble
    .tap { fastq4TrinityAssembleDenovo }


// Indexing and preprocessing

process getFaidx {

    label "samtools"
    label "small_task"

    tag { name }

    input:
    set val(name), file(genome) from genomes4GetFaidx

    output:
    set val(name), file(genome), file("${genome}.fai") into genomesWithFaidx

    """
    samtools faidx "${genome}"
    """
}

genomesWithFaidx.set { genomesWithFaidx4TidyBams }


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


/*
 * Index genome for STAR
 */

process getStarIndex {

    label "star"
    label "medium_task"

    tag { name }

    when:
    params.fastq

    input:
    set val(name), file(fasta), file(gff) from genomes4StarIndexWithGff

    output:
    set val(name), file("index") into starIndices

    script:
    // If no gff was in known_sites, use it, otherwise dont
    def sjdb = gff.name != 'WAS_NULL' ? "--sjdbGTFfile ${gff} --sjdbOverhang 149 " : ''

    // --sjdbGTFfeatureExon should there be an option to make this CDS?
    """
    mkdir -p index
    STAR \
      --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir index \
      --genomeFastaFiles ${fasta} \
      --genomeSAindexNbases 11 \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFfeatureExon exon \
      ${sjdb}
    """
}

starIndices.into {
    starIndices4StarFindNovelSpliceSites;
    starIndices4StarAlignReads;
}


/*
 * Perform first pass for STAR.
 *
 * This helps refine intron splice sites in the alignments.
 */
process starFindNovelSpliceSites {

    label "star"
    label "medium_task"

    tag "${name} - ${read_group}"

    when:
    params.fastq

    input:
    set val(name), file("index"), val(read_group),
        file(r1s), file(r2s) from starIndices4StarFindNovelSpliceSites
            .combine(fastq4StarFindNovelSpliceSites)
            .groupTuple(by: [0, 1, 2])

    output:
    set val(name), val(read_group),
        file("${name}_${read_group}.SJ.out.tab") into starNovelSpliceSites

    script:
    def r1_joined = r1s.join(',')
    def r2_joined = r2s.join(',')

    """
    STAR \
      --runThreadN "${task.cpus}" \
      --readFilesCommand zcat \
      --genomeDir "index" \
      --outSAMtype None \
      --outSAMmode None \
      --outSJfilterReads All \
      --outSJfilterOverhangMin 5 1 1 1 \
      --outSJfilterCountUniqueMin 10 5 5 5 \
      --outSJfilterDistToOtherSJmin 5 0 0 3 \
      --outSJfilterIntronMaxVsReadN 0 1 500 5000 10000 20000 \
      --alignIntronMin 5 \
      --alignIntronMax 10000 \
      --alignSJoverhangMin 1 \
      --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds No \
      --outFileNamePrefix "${name}_${read_group}." \
      --readFilesIn "${r1_joined}" "${r2_joined}"
    """
}


/*
 * Perform second pass STAR alignment
 */
process starAlignReads {

    label "star"
    label "medium_task"

    tag "${name} - ${read_group}"

    when:
    params.fastq

    input:
    set val(name), file("index"), val(read_group),
        file(r1s), file(r2s) from starIndices4StarAlignReads
            .combine(fastq4StarAlignReads)
            .groupTuple(by: [0, 1, 2])
    file "*SJ.out.tab" from starNovelSpliceSites
        .map { n, rg, f -> f }
        .collect()

    output:
    set val(name), val(read_group), file("${name}_${read_group}.bam") into alignedReads

    script:
    def r1_joined = r1s.join(',')
    def r2_joined = r2s.join(',')
    """
    STAR \
      --runThreadN ${task.cpus} \
      --readFilesCommand zcat \
      --genomeDir "index" \
      --sjdbFileChrStartEnd *SJ.out.tab \
      --outSAMtype BAM Unsorted \
      --outBAMcompression 1 \
      --outSJfilterReads All \
      --outSJfilterCountUniqueMin 10 5 5 5 \
      --outSJfilterIntronMaxVsReadN 0 1 500 5000 10000 20000 \
      --alignIntronMin 5 \
      --alignIntronMax 10000 \
      --alignSJoverhangMin 10 \
      --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds No \
      --outFilterType BySJout \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 10 \
      --outFilterMismatchNoverLmax 0.2 \
      --outMultimapperOrder Random \
      --outSAMattributes All \
      --outSAMstrandField intronMotif\
      --outSAMattrIHstart 0 \
      --outSAMmapqUnique 50 \
      --outFileNamePrefix "${name}_${read_group}." \
      --readFilesIn "${r1_joined}" "${r2_joined}"

    mv "${name}_${read_group}.Aligned.out.bam" "${name}_${read_group}.bam"
    """
}


process tidyBams {

    label "samtools"
    label "small_task"
    publishDir "${params.outdir}/aligned_reads"

    tag "${name} - ${read_group}"

    when:
    params.fastq

    input:
    set val(name), val(read_group), file("my.bam"),
        file(genome), file(faidx) from alignedReads
            .combine(genomesWithFaidx4TidyBams, by: 0)

    output:
    set val(name), val(read_group), file("${name}_${read_group}.bam"),
        file("${name}_${read_group}.bam.bai") into tidiedBams

    script:
    """
    samtools view \
        -uT "${faidx}" \
        "my.bam" \
    | samtools sort \
        -O BAM \
        -@ "${task.cpus}" \
        -l 9 \
        -o "${name}_${read_group}.bam"

    samtools index "${name}_${read_group}.bam"
    """
}

tidiedBams.set { bams4AssembleStringtie }


/*
 * Assembly transcripts with Stringtie
 */
process assembleStringtie {

    label "stringtie"
    label "medium_task"

    tag "${name} - ${read_group}"

    when:
    params.fastq

    input:
    set val(name), val(read_group), file(bam), file(bai),
        file(fasta), file(gff) from bams4AssembleStringtie
            .combine(genomes4AssembleStringtie, by: 0)

    output:
    set val(name), val(read_group),
        file("${name}_${read_group}.gtf") into stringtieAssembledTranscripts

    script:
    def stranded = params.fr ? '' : "--rf"
    def known = gff.name != 'WAS_NULL' ? "-G ${gff}" : ''
    """
    stringtie \
      -p "${task.cpus}" \
      ${stranded} \
      ${known} \
      -o "${name}_${read_group}.gtf" \
      -m 200 \
      "${bam}"
    """
}


/*
 * Combine stringtie annotations from multiple bams.
 */
process mergeStrintie {

    label "stringtie"
    label "medium_task"
    publishDir "${params.outdir}/stringtie"

    tag "${name}"

    when:
    params.fastq

    input:
    set val(name), file("*gtf"), file(gff) from stringtieAssembledTranscripts
        .map { n, rg, f -> [n, f] }
        .groupTuple(by: 0)
        .combine( genomes4MergeStringtie.map {n, f, g -> [n, g]}, by: 0 )

    output:
    set val(name), file("${name}.gtf") into stringtieMergedTranscripts

    script:
    def known = gff.name != 'WAS_NULL' ? "-G ${gff}" : ''
    """
    stringtie \
      -p "${task.cpus}" \
      ${known} \
      --merge \
      -o "${name}.gtf" \
      *gtf
    """
}


stringtieMergedTranscripts.into {
    stringtieMergedTranscripts4CodingQuarry;
    stringtieMergedTranscripts4CodingQuarryPM;
}


process runCodingQuarry {

    label "codingquarry"
    label "medium_task"
    publishDir "${params.outdir}/codingquarry"

    tag "${name}"

    when:
    params.fastq && !params.notfungus

    input:
    set val(name), file("transcripts.gtf"),
        file("genome.fasta") from stringtieMergedTranscripts4CodingQuarry
            .combine(genomes4RunCodingQuarry, by: 0)

    output:
    set val(name), file("${name}") into codingQuarryPredictions

    script:
    """
    grep -v "^#" transcripts.gtf > transcripts.tmp.gtf
    CufflinksGTF_to_CodingQuarryGFF3.py transcripts.tmp.gtf > transcripts.gff3

    CodingQuarry -f genome.fasta -t transcripts.gff3 -p "${task.cpus}"

      \${QUARRY_PATH}/scripts/fastaTranslate.py out/Predicted_CDS.fa \
    | sed 's/*\$//g' > Predicted_Proteins.faa

      \${QUARRY_PATH}/scripts/gene_errors_Xs.py Predicted_Proteins.faa out/Predicted_Proteins.faa

    mv out "${name}"
    """
}

codingQuarryPredictions.into {
    codingQuarryPredictions4SignalP;
    codingQuarryPredictions4PM;
}


process getCodingQuarrySignalP {

    label "signalp"
    label "medium_task"
    publishDir "${params.outdir}/codingquarry"

    tag "${name}"

    when:
    params.fastq && !params.notfungus && params.signalp

    input:
    set val(name), file("cq") from codingQuarryPredictions4SignalP

    output:
    set val(name), file("${name}.signalp5") into codingQuarryPredictionsSecreted

    script:
    """
    mkdir tmp
    signalp \
      -fasta "cq/Predicted_Proteins.faa" \
      -prefix "${name}" \
      -org euk \
      -tmp tmp

    mv "${name}_summary.signalp5" "${name}.signalp5"

    rm -rf -- tmp
    """
}


/*
*/
process runCodingQuarryPM {

    label "codingquarry"
    label "medium_task"
    publishDir "${params.outdir}/codingquarry"

    tag "${name}"

    when:
    params.fastq && !params.notfungus && params.signalp

    input:
    set val(name), file("transcripts.gtf"), file("genome.fasta"),
        file("first"), file("secretome.signalp5") from stringtieMergedTranscripts4CodingQuarryPM
            .combine(genomes4RunCodingQuarryPM, by: 0)
            .combine(codingQuarryPredictions4PM, by: 0)
            .combine(codingQuarryPredictionsSecreted, by: 0)

    output:
    set val(name), file("${name}") into codingQuarryPMPredictions

    script:
    """
    grep -v "^#" transcripts.gtf > transcripts.tmp.gtf
    CufflinksGTF_to_CodingQuarryGFF3.py transcripts.tmp.gtf > transcripts.gff3

    gawk '
      BEGIN {
        OFS=" "
      }
      \$2 ~ /^SP/ {
        match(\$0, /CS pos: ([0-9]*)-/, x)
        print \$1, x[1]
      }
    ' < "secretome.signalp5" > secretome.txt

    CodingQuarry \
      -f genome.fasta \
      -t transcripts.gff3 \
      -2 "first/PredictedPass.gff3" \
      -p "${task.cpus}" \
      -g secretome.txt \
      -h

    mv out "${name}"
    """
}

// RNAseq alignments and assembly
/*
process trinityAssembleTranscripts {

    label "trinity"
    label "big_task"

    input:

    output:

    script:
    // NB genome guided uses bams to assemble the reads, needs a separate process.
    // jaccard_clip should be useful for fungi. make it an option to disable?.
    """
    Trinity \
      --seqtype fq \
      --max_memory "${task.memory}" \
      --CPU "${task.cpus}" \
      --jaccard_clip \
      --output trinity_assembly \
    """
}

*/

// 1 align transcripts to all genomes
/*
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
*/

/*
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
*/

// Align with minimap2


// 2 align proteins to all genomes
/*
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
*/

// 3 align remote proteins to all genomes
/*
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
*/

// 4 run busco on all genomes
/*
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
*/

// 5 Run genemark, pasa, braker2, codingquarry on all genomes using previous steps.

// Combine transcript input for genemark`
/*
process runGenemark {

    label "genemark"
    label "small_task"

    tag { name }

    when:
    params.genemark

    input:
    genome
    transcript_intron_coords (e.g. from star)

    output:
    genemark.gtf

    script:
    def use_fungus = params.notfungus ? '' : '--fungus '

    """
    gmes_petap.pl \
      --soft_mask 100 \
      --ET ${transcript_intron_coords} \
      --evidence \
      ${use_fungus}
    """
}
*/

/*
process runBraker {

    label "braker"
    label "small_task"

    tag { name }

    input:
    genome
    genemark_gtf
    spaln_protein_alignments
    transcripts
    bams

    output:
    braker.gff
    species file

    script:
    def species = "Pnodorum"
    def use_softmasking = params.softmasked ? "--softmasking --UTR=on"

    // Input GFFs have specific format requirements.
    // NB braker uses spaln output 0
    // Might need to separate bams into different strands? using --stranded
    // --crf might be good to try for the "reference" isolate?.
    // --AUGUSTUS_CONFIG_PATH can set this instead of env-var.
    """
    braker.pl \
      --cores="${task.cpus}" \
      --species="${species}" \
      --alternatives-from-evidence=true \
      --genome="${genome}" \
      --bam="one.bam,two.bam" \
      --hints="one.gff,two.gff" \
      --prot_aln=proteins.gff \
      --prg=spaln \
      --geneMarkGtf=file.gtf \
      --skipGetAnnoFromFasta \
      --verbosity=3 \
      --splice_sites=GTAG,ATAC \
    """
}
*/

// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats
