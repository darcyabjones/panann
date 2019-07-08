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

params.genome_alignment = false

params.busco_lineage = false
params.augustus_config = false

// RNAseq params
params.fastq = false
params.crams = false
params.fr = false
params.cufflinks_gtf = false

params.notfungus = false
params.genemark = false
params.signalp = false
params.notrinity = false
params.nostar = false


/*
 * Sanitise the input
 */

def stranded = params.fr ? "fr" : "rf"

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
    Channel
        .fromPath(params.known_sites, checkIfExists: true, type: "file")
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { (!is_null(it.name) && !is_null(it.gff3)) }
        .map {[
            it.name,
            file(it.gff3, checkIfExists: true),
        ]}
        .unique()
        .into { knownSites; knownSites4CheckNoDups }

    /*
    knownSites4CheckNoDups
        .groupTuple(by: 0)
        .filter { name, gffs -> gffs.length > 1 }
        .map { n, g -> n }
        .collect()
        .set { ks_duplicates }

    if ( ks_duplicates.length != 0 ) {
        log.error "There is more than one "known" GFF annotation " +
            "provided for name ${name}. Please merge them into one " +
            " file. E.G. using `genometools gff3`"
        exit 1
    }
    */

} else {
    knownSites = Channel.empty()
}


if ( params.transcripts ) {
    Channel
        .fromPath(params.transcripts, checkIfExists: true, type: "file")
        .set { transcripts }

} else {
    transcripts = Channel.empty()
}


if ( params.proteins ) {
    Channel
        .fromPath(params.proteins, checkIfExists: true, type: "file")
        .set { proteins }

} else {
    proteins = Channel.empty()
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

        script:
        """
        cp -r \${AUGUSTUS_CONFIG_PATH} ./config
        """
    }
}


if ( params.fastq ) {
    Channel
        .fromPath(params.fastq, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { !is_null(it.r1) || !is_null(it.r2) }
        .map { it ->
            if ( is_null(it.read_group) ) {
                log.error "A fastq pair has no read_group set."
                log.error "The offending line is ${it}"
                exit 1
            };
            if ( is_null(it.r1) || is_null(it.r2) ) {
                log.error "A fastq pair has one member of the pair unset."
                log.error "The pipeline only currently supports stranded " +
                    "paired-RNAseq."
                exit 1
            };
            [
                it.read_group,
                !is_null(it.name) ? it.name : 'WAS_NULL',
                file(it.r1, checkIfExists: true),
                file(it.r2, checkIfExists: true),
                !is_null(it.stranded) ? it.stranded : stranded
            ]
        }
        .unique()
        .set { fastq }

} else {
    fastq = Channel.empty()
}


if ( params.crams ) {
    Channel
        .fromPath(params.cram, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { !is_null(it.cram) }
        .map { it ->
            if ( is_null(it.read_group) ) {
                log.error "The cram file ${it.cram.name} doesn't have a " \
                    "read_group set."
                exit 1
            };
            if ( is_null(is.name) ) {
                log.error "The cram file ${is.cram.name} with read_group " \
                    "${it.read_group}, doesn't have a genome name set."
                exit 1
            };
            [
                it.name,
                it.read_group,
                file(it.cram, checkIfExists: true),
                !is_null(it.stranded) ? it.stranded : stranded
            ]
        }
        .unique()
        .set { crams }

} else {
    crams = Channel.empty()
}


process getFaidx {

    label "samtools"
    label "small_task"

    tag { name }

    input:
    set val(name), file(genome) from genomes

    output:
    set val(name), file(genome), file("${genome}.fai") into genomesWithFaidx

    script:
    """
    samtools faidx "${genome}"
    """
}


genomesWithFaidx.into {
    genomes4KnownSites;
    genomes4SpalnIndex;
    genomes4GmapIndex;
    genomes4UserCrams;
    genomes4TidyBams;
    genomes4TidyFilteredBams;
    genomes4ExtractSpliceSites;
    genomes4Busco;
    genomes4RunGenemark;
    genomes4RunCodingQuarry;
    genomes4RunCodingQuarryPM;
}


if ( params.known_sites ) {
    genomes4KnownSites
        .join( knownSites, by: 0, remainder: true)
        .map { n, f, i, g ->
            if ( is_null(f) ) {
                log.error "The known site file ${g.name} with name ${n} " +
                    "doesn't match any of the provided genomes."
                log.error "Please check that the name in the table matches " +
                    "the basename of one of the genomes."
                exit 1
            };
            [n, f, i, !is_null(g) ? g : file('WAS_NULL')]
        }
        .set { genomesWithKnownSites }

} else {
    genomes4KnownSites
        .map { n, f, i -> [n, f, i, file('WAS_NULL')] }
        .set { genomesWithKnownSites }
}

genomesWithKnownSites.into {
    genomes4StarIndex;
    genomes4AssembleStringtie;
    genomes4MergeStringtie;
    genomes4RunPASA;
}



crams
    .combine(genomes4UserCrams, by: 0)
    .map { n, rg, c, st, fa, fi ->
        if (is_null(fa)) {
            log.error "The provided cram ${c.name} could not be assigned to a" +
                " genome based on the name column."
            log.error "Please check that the name in the table matches "+
                "the basename of one of the genomes."
            exit 1
        };
        [n, rg, fa, fi, c, st]
    }
    .set { userCrams }


fastq
    .tap { fastqNoGenome; fastq4Alignment }
    .filter { rg, n, r1, r2, st -> n != 'WAS_NULL' }
    .map { rg, n, r1, r2, st -> [n, rg] }
    .set { fastq4TrinityAssembleGuided }

fastqNoGenome
    .map { rg, n, r1, r2, st -> [rg, r1, r2, st] }
    .unique()
    .set { fastq4TrinityAssemble }

fastq4Alignment
    .map { rg, n, r1, r2, st -> [rg, r1, r2, st] }
    .unique()
    .into {
        fastq4StarFindNovelSpliceSites;
        fastq4StarAlignReads;
    }


//
// Indexing and preprocessing
//


process getUnivec {

    label "download"
    label "small_task"

    output:
    file "univec.fasta" into univec

    script:
    """
    wget -O univec.fasta ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
    """
}


process getSpalnIndex {

    label "spaln"
    label "small_task"

    tag { name }

    input:
    set val(name), file(genome), file(faidx) from genomes4SpalnIndex

    output:
    set val(name),
        file("${genome.baseName}.bkn"),
        file("${genome.baseName}.ent"),
        file("${genome.baseName}.idx"),
        file("${genome.baseName}.bkp"),
        file("${genome.baseName}.grp"),
        file("${genome.baseName}.seq") into spalnIndices

    script:
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
    set val(name), file(genome), file(faidx) from genomes4GmapIndex

    output:
    set val(name), file("db") into gmapIndices

    script:
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
    params.fastq && !params.nostar

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file(gff) from genomes4StarIndex

    output:
    set val(name), file("index"), file(fasta), file(faidx) into starIndices

    script:
    // If no gff was in known_sites, use it, otherwise dont
    def sjdb = gff.name != 'WAS_NULL' ? "--sjdbGTFfile ${gff} --sjdbOverhang 149 " : ''
    // Possibly need option to set this to CDS, incase user input doesn't have exons?
    def exon_feature = "exon"

    """
    mkdir -p index
    STAR \
      --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir index \
      --genomeFastaFiles "${fasta}" \
      --genomeSAindexNbases 11 \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFfeatureExon "${exon_feature}" \
      ${sjdb}
    """
}

starIndices.into {
    starIndices4StarFindNovelSpliceSites;
    starIndices4StarAlignReads;
}


//
// RNAseq alignments and assembly
//

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
    !params.nostar

    input:
    set val(name),
        file("index"),
        file(fasta),
        file(faidx),
        val(read_group),
        file(r1s),
        file(r2s),
        val(strand) from starIndices4StarFindNovelSpliceSites
            .combine(fastq4StarFindNovelSpliceSites)
            .groupTuple(by: [0, 1, 2, 3, 4])

    output:
    set val(name),
        val(read_group),
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
    publishDir "${params.outdir}/aligned_reads"

    tag "${name} - ${read_group}"

    when:
    !params.nostar

    input:
    set val(name),
        val(read_group),
        file("index"),
        file(fasta),
        file(faidx),
        file(r1s),
        file(r2s),
        val(strand),
        file("*SJ.out.tab") from starIndices4StarAlignReads
            .combine(fastq4StarAlignReads)
            .map { n, i, fa, fi, rg, r1, r2, st -> [n, rg, i, fa, fi, r1, r2, st]}
            .combine(starNovelSpliceSites, by: [0, 1])
            .groupTuple(by: [0, 1, 2, 3, 4])

    output:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file("${name}_${read_group}.cram"),
        val(strand) into alignedReads

    script:
    // todo assert all strand is same?
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
      --outSAMstrandField intronMotif \
      --outSAMattrIHstart 0 \
      --outSAMmapqUnique 50 \
      --outFileNamePrefix "${name}_${read_group}." \
      --readFilesIn "${r1_joined}" "${r2_joined}"

    mkdir tmp
    samtools view \
        -u \
        -C \
        -T "${fasta}" \
        "${name}_${read_group}.Aligned.out.bam" \
    | samtools sort \
        -O cram \
        -@ "${task.cpus}" \
        -T tmp \
        -l 9 \
        -o "${name}_${read_group}.cram"

    rm -rf -- tmp
    rm -f *.bam
    """
}

alignedReads
    .mix(userCrams)
    .into {
        crams4AssembleStringtie;
        crams4TrinityAssembleGuided;
        crams4FilterCrams;
    }


process filterCrams {

    label "augustus"
    label "medium_task"

    tag "${name} - ${read_group}"

    input:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file(cram),
        val(strand) from crams4FilterCrams

    output:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file("${cram.baseName}_filtered.cram"),
        val(strand) into filteredCrams

    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file("${name}_${read_group}_forward.cram"),
        file("${name}_${read_group}_reverse.cram"),
        file("${name}_${read_group}_unpaired.cram") into splitReadCrams

    script:
    if (strand == "fr") {
        fst = "forward"
        snd = "reverse"
    } else {
        fst = "reverse"
        snd = "forward"
    }

    """
    mkdir tmp

    # Convert cram to bam and sort by read name.
    # The -n is the important bit here
    samtools view \
        -u \
        -b \
        -T "${fasta}" \
        "${cram}" \
    | samtools sort \
        -O BAM \
        -@ "${task.cpus}" \
        -T tmp \
        -n \
        -l 9 \
        -o "my.bam"

    filterBam \
      --uniq \
      --paired \
      --pairwiseAlignment \
      --in "my.bam" \
      --out "my_filtered.bam"

    # Convert bam to cram and sort by position.
    samtools view \
        -u \
        -C \
        -T "${fasta}" \
        "my_filtered.bam" \
    | samtools sort \
        -O cram \
        -@ "${task.cpus}" \
        -T tmp \
        -l 9 \
        -o "${name}_${read_group}_filtered.cram"

    samtools view -f 65 "${name}_${read_group}_filtered.cram" > "${name}_${read_group}_${fst}.cram"
    samtools view -f 128 "${name}_${read_group}_filtered.cram" > "${name}_${read_group}_${snd}.cram"
    samtools view -F 193 "${name}_${read_group}_filtered.cram" > "${name}_${read_group}_unpaired.cram"

    rm -rf -- tmp
    rm -f my.bam my_filtered.bam
    """
}

filteredCrams.set { filteredCrams4ExtractSpliceSites }


process extractSpliceSites {

    label "braker"
    label "small_task"

    tag "${name} - ${read_group}"

    input:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file(cram),
        val(strand) from filteredCrams4ExtractSpliceSites

    output:
    set val(name),
        val(read_group),
        file("${name}_${read_group}_introns.gff3") into spliceSites

    script:
    """
    # Convert cram to bam.
    samtools view \
        -b \
        -T "${fasta}" \
        -@ "${task.cpus}" \
        -o "tmp.bam"
        "${cram}"

    bam2hints \
      --intronsonly \
      --in="tmp.bam" \
      --out="tmp.gff3"

    # braker panics if the genome has descriptions
    sed -r 's/^(>[^[:space:]]*).*\$/\\1/' "${fasta}" > tmp.fasta 

    filterIntronsFindStrand.pl \
      tmp.fasta \
      tmp.gff3 \
      --score \
      > "${name}_${read_group}_introns.gff3"

    rm -f tmp.fasta tmp.gff3 tmp.bam
    """
}

spliceSites.set { spliceSites4RunGenemark }


/*
 * Assemble transcripts with Stringtie
 */
process assembleStringtie {

    label "stringtie"
    label "medium_task"

    tag "${name} - ${read_group}"

    input:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file(cram),
        val(strand),
        file(gff) from crams4AssembleStringtie
            .join(
                genomes4AssembleStringtie.map { n, f, i, g -> [n, g] },
                by: 0,
                remainder: false
            )

    output:
    set val(name),
        val(read_group),
        file("${name}_${read_group}.gtf") into stringtieAssembledTranscripts

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
      -o "${name}_${read_group}.gtf" \
      -m 200 \
      "tmp.bam"

    rm -f tmp.bam
    """
}


/*
 * Combine stringtie annotations from multiple bams.
 */
process mergeStringtie {

    label "stringtie"
    label "medium_task"
    publishDir "${params.outdir}/stringtie"

    tag "${name}"

    input:
    set val(name), file("*gtf"), file(gff) from stringtieAssembledTranscripts
        .map { n, rg, f -> [n, f] }
        .groupTuple(by: 0)
        .join( genomes4MergeStringtie.map {n, f, i, g -> [n, g]}, by: 0 )

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
    stringtieMergedTranscripts4PASA;
}


/*
 * Assemble reads into transcripts with trinity.
 */
process trinityAssemble {

    label "trinity"
    label "big_task"
    publishDir "${params.outdir}/trinity"

    tag "${read_group}"

    when:
    !params.notrinity

    input:
    set val(read_group),
        file(r1s),
        file(r2s),
        val(strand) from fastq4TrinityAssemble
            .groupTuple(by: 0)

    output:
    set val(read_group), file("${read_group}.fasta") into assembledTranscripts

    script:
    def r1_joined = r1s.join(',')
    def r2_joined = r2s.join(',')
    def use_jaccard = params.notfungus ? '' : "--jaccard_clip "
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

    mv trinity_assembly/Trinity.fasta "${read_group}.fasta"
    rm -rf -- trinity_assembly
    """
}


/*
 * Perform genome guided assembly with Trinity.
 */
process trinityAssembleGuided {

    label "trinity"
    label "big_task"
    publishDir "${params.outdir}/trinity"

    tag "${name}"

    when:
    !params.notrinity

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file("*cram"),
        val(strand) from fastq4TrinityAssembleGuided
            .join(crams4TrinityAssembleGuided, by: [0, 1])
            .map { n, rg, fa, fi, c, st -> [n, fa, fi, c, st] }
            .groupTuple(by: [0, 1, 2])

    output:
    set val(name), val(read_group), file("${name}_${read_group}.fasta") into guidedAssembledTranscripts

    script:
    def use_jaccard = params.notfungus ? '' : "--jaccard_clip "
    def strand_flag = strand == "fr" ? "--SS_lib_type FR " : "--SS_lib_type RF "

    """
    # Convert each cram to a bam
    for c in *cram; do
        samtools view \
            -b \
            -T "${fasta}" \
            -@ "${task.cpus}" \
            -o "\${c%cram}bam"
            "\${c}"
    done

    samtools merge -r "merged.bam" *bam
    samtools index "merged.bam"

    Trinity \
      --max_memory "${task.memory.toGiga()}G" \
      --CPU "${task.cpus}" \
      ${use_jaccard} \
      ${strand_flag} \
      --output trinity_assembly \
      --genome_guided_max_intron 15000 \
      --genome_guided_bam "merged.bam"

    mv trinity_assembly/Trinity-GG.fasta "${name}.fasta"
    rm -rf -- *bam merged.bam.bai trinity_assembly
    """
}


//
// 1 align transcripts to all genomes
//


// Might be good to keep track of which assemblies come from trinity denovo?
// Pasa can use this somehow.
process combineTranscripts {

    label "python3"
    label "small_task"

    input:
    file "*fasta" from transcripts
        .mix(assembledTranscripts.map { rg, f -> f } )
        .mix(guidedAssembledTranscripts.map { n, rg, f -> f } )
        .collect()

    output:
    file "transcripts.fasta" into combinedTranscripts
    file "transcripts.tsv" into combinedTranscriptsMap

    script:
    """
    unique_rename_fasta.py \
      --infiles *fasta \
      --outfile "transcripts.fasta" \
      --map "transcripts.tsv"
    """
}


/*
 */
process cleanTranscripts {
    label "pasa"
    label "small_task"

    input:
    file "transcripts.fasta" from combinedTranscripts
    file "univec.fasta" from univec

    output:
    set file("transcripts.fasta"),
        file("transcripts.fasta.cln"),
        file("transcripts.fasta.clean") into cleanedTranscripts

    script:
    """
    # this user thing is needed for seqclean. Unknown reasons
    export USER="root"
    seqclean "transcripts.fasta" -v "univec.fasta"
    """
}

cleanedTranscripts.into {
    transcripts4AlignSpalnTranscripts;
    transcripts4AlignGmapTranscripts;
    transcripts4RunPasa;
}


/*
 */
process alignSpalnTranscripts {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/spaln"

    tag { name }

    input:
    set file(fasta),
        file(fasta_cln),
        file(fasta_clean) from transcripts4AlignSpalnTranscripts

    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices4AlignSpalnTranscripts

    output:
    set val(name), file("${name}_transcripts.gff3") into spalnAlignedTranscripts

    script:
    def species = "Dothideo"

    """
    spaln \
      -L \
      -M3 \
      -O2 \
      -Q7 \
      -S3 \
      -T${species} \
      -yX \
      -yS \
      -ya2 \
      -t ${task.cpus} \
      -d "${name}" \
      "${fasta_clean}" \
    > "${name}_transcripts.gff3"
    """
}


/*
 */
process alignGmapTranscripts {

    label "gmap"
    label "medium_task"

    publishDir "${params.outdir}/aligned/gmap"

    tag { name }

    input:
    set file(fasta),
        file(fasta_cln),
        file(fasta_clean) from transcripts4AlignGmapTranscripts

    set val(name), file("db") from gmapIndices

    output:
    set val(name), file("${name}_transcripts.gff3") into gmapAlignedTranscripts

    script:
    def min_intronlength = 6
    def max_intronlength_middle = 20000
    def max_intronlength_ends = 10000
    def trim_end_exons = 12
    def microexon_spliceprob = 0.95
    def canonical_mode = 1
    def cross_species = "" // "--cross-species "

    """
    gmap \
      --npaths=1 \
      --chimera-margin=50 \
      --min-intronlength="${min_intronlength}" \
      --max-intronlength-middle="${max_intronlength_middle}" \
      --max-intronlength-ends="${max_intronlength_ends}" \
      --trim-end-exons="${trim_end_exons}" \
      --microexon-spliceprob="${microexon_spliceprob}" \
      --canonical-mode="${canonical_mode}" \
      ${cross_species} \
      --format=gff3_match_cdna \
      --nthreads "${task.cpus}" \
      -D db \
      -d "${name}" \
      ${fasta_clean} \
    > ${name}_transcripts.gff3
    """
}

gmapAlignedTranscripts.set { gmapAlignedTranscripts4RunPASA }


//
// 2 align proteins to all genomes
//


process combineProteins {

    label "python3"
    label "small_task"

    input:
    file "*fasta" from proteins.collect()

    output:
    file "proteins.fasta" into combinedProteins
    file "proteins.tsv" into combinedProteinsMap

    script:
    """
    unique_rename_fasta.py \
      --infiles *fasta \
      --outfile "proteins.fasta" \
      --map "proteins.tsv"
    """
}

/*
 */
process alignSpalnProteins {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/spaln"

    tag { name }

    input:
    file "proteins.fasta" from combinedProteins
    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from spalnIndices4AlignSpalnProteins

    output:
    set val(name), file("${name}_proteins.gff3")

    script:
    """
    spaln \
      -L \
      -M3 \
      -O2 \
      -Q7 \
      -TDothideo \
      -ya2 \
      -t ${task.cpus} \
      -a \
      -d "${name}" \
      "proteins.fasta" \
    > "${name}_proteins.gff3"
    """
}


//
// 3 run busco on all genomes
//

/*
 */
process runBusco {
    label "busco"
    label "medium_task"
    tag { name }

    publishDir "${params.outdir}/busco"

    when:
    params.busco_lineage

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

//
// 5 Run genemark, pasa, braker2, codingquarry on all genomes using previous steps.
//

/*
 * Add the stringtie assembled and gmap aligned transcripts to the input channel.
 */
if ( params.fastq || params.crams ) {
    genomes4RunPASA
        .join(stringtieMergedTranscripts4PASA, by: 0, remainder: true)
        .map { n, f, i, ks, sg ->
            [ n, f, i, ks, !is_null(sg) ? sg : file("WAS_NULL") ]
        }
        .join(gmapAlignedTranscripts4RunPASA, by: 0, remainder: true)
        .set { processed4RunPASA }

} else {
    genomes4RunPASA
        .map { n, f, i, ks -> [ n, gf, gi, ks, file("WAS_NULL") ]}
        .join(gmapAlignedTranscripts4RunPASA, by: 0, remainder: true)
        .set { processed4RunPASA }
}


/*
 * Run pasa
 */
process runPASA {

    label "pasa"
    label "medium_task"

    tag { name }

    input:
    set val(name),
        file(genome_fasta),
        file(genome_faidx),
        file(known_sites),
        file(stringtie_gtf),
        file(gmap_aligned) from processed4RunPASA

    set file(transcripts_fasta),
        file(transcripts_fasta_cln),
        file(transcripts_fasta_clean) from transcripts4RunPasa

    output:
    set val(name),
        file("${name}.gff3"),
        file("${name}_cds.fna"),
        file("${name}_protein.faa") into pasaPredictions

    script:
    def use_stringent = params.notfungus ? '' : "--stringent_alignment_overlap 30.0 "
    def use_stringtie = stringtie_gtf.name == "WAS_NULL" ? '' : "--trans_gtf ${stringtie_gtf} "
    def use_known = known_sites.name == "WAS_NULL" ? '' : "-L -annots ${known_sites} "

    """
    echo "DATABASE=\${PWD}/pasa.sqlite" > align_assembly.config
    echo "validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80" >> align_assembly.config
    echo "validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=90" >> align_assembly.config
    echo "subcluster_builder.dbi:-m=50" >> align_assembly.config

    Launch_PASA_pipeline.pl \
      --config align_assembly.config \
      --create \
      --run \
      --genome "${genome_fasta}" \
      --transcripts "${transcripts_fasta_clean}" \
      --IMPORT_CUSTOM_ALIGNMENTS_GFF3 "${gmap_aligned}" \
      -T -u "${transcripts_fasta}" \
      --MAX_INTRON_LENGTH 50000 \
      --ALIGNERS blat \
      --CPU "${task.cpus}" \
      --transcribed_is_aligned_orient \
      --TRANSDECODER \
      ${use_stringent} \
      ${use_stringtie} \
      ${use_known}

    pasa_asmbls_to_training_set.dbi \
      --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta \
      --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3

    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 \${PWD}/${name}.gff3
    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.cds \${PWD}/${name}_cds.fna
    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.pep \${PWD}/${name}_protein.faa
    """
}


// To find "complete" cdss for training augustus...
// awk '\$0 ~ /^>.*type:complete/ {
//   r = gensub(/^>[[:space:]]?([^[:space:]]+).*/, "\\1", "g", $0);
//   print r;
// }' ${name}_cds.fna > complete.orfs
//
// grep -F -f complete.orfs "${name}.gff3" \
// | awk '$3 == "exon" || $3 == "CDS"' \
// | 's/cds\.//; s/\.exon[[:digit:]]*//' \
// | sort -s -n -k 1,1 -k 4 -k 9 \
// > training_set_complete.gff3


/*
 */
process runGenemark {

    label "genemarkes"
    label "small_task"

    tag { name }

    when:
    params.genemark

    input:
    set val(name),
        file(genome),
        file(faidx),
        file("*introns.gff3") from genomes4RunGenemark
            .join(
                spliceSites4RunGenemark
                    .map {n, rg, introns -> [n, introns]}
                    .groupTuple(by: 0),
                by: 0
            )

    output:
    file "genemark.gtf"

    script:
    def use_fungus = params.notfungus ? '' : '--fungus '

    """
    sort -k1,1V -k4,4n -k5,5rn -k3,3r *introns.gff3 > hints.gff3

    gmes_petap.pl \
      --cores "${task.cpus}" \
      --soft_mask 100 \
      --ET "hints.gff3" \
      ${use_fungus} \
      --sequence "${genome}"
    """
}

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


process runCodingQuarry {

    label "codingquarry"
    label "medium_task"
    publishDir "${params.outdir}/codingquarry"

    tag "${name}"

    when:
    !params.notfungus

    input:
    set val(name),
        file("transcripts.gtf"),
        file("genome.fasta"),
        file("genomes.fasta.fai") from stringtieMergedTranscripts4CodingQuarry
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


/*
 * CQPM looks for genes that might not be predicted main set because of
 * genome compartmentalisation.
 * NOTE: This fails if there are fewer than 500 secreted genes to train from.
 */
process getCodingQuarrySignalP {

    label "signalp"
    label "medium_task"
    publishDir "${params.outdir}/codingquarry"

    tag "${name}"

    when:
    !params.notfungus && params.signalp

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
    !params.notfungus && params.signalp

    input:
    set val(name),
        file("transcripts.gtf"),
        file("genome.fasta"),
        file("genome.fasta.fai"),
        file("first"),
        file("secretome.signalp5") from stringtieMergedTranscripts4CodingQuarryPM
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

// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats
