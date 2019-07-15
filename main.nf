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
        .fromPath(params.crams, checkIfExists: true)
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { !is_null(it.cram) }
        .map { it ->
            if ( is_null(it.read_group) ) {
                log.error "The cram file ${it.cram.name} doesn't have a " \
                    "read_group set."
                exit 1
            };
            if ( is_null(it.name) ) {
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
    set val(name), file("orig.fa") from genomes

    output:
    set val(name), file("${name}.fasta"), file("${name}.fasta.fai") into genomesWithFaidx

    script:
    """
    # braker panics if the genome has descriptions
    sed -r 's/^(>[^[:space:]]*).*\$/\\1/' orig.fa > "${name}.fasta"

    samtools faidx "${name}.fasta"
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
    genomes4AlignGemomaCDSParts;
    genomes4RunGemoma;
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
    genomes4ExtractGemomaCDSParts;
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
    publishDir "${params.outdir}/aligned/${name}"

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
        crams4ExtractAugustusRnaseqHints;
        crams4ExtractGemomaRnaseqHints;
    }




/*
 * Assemble transcripts with Stringtie
 */
process assembleStringtie {

    label "stringtie"
    label "medium_task"
    publishDir "${params.outdir}/aligned/${name}"

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
        file("${name}_${read_group}_stringtie.gtf") into stringtieAssembledTranscripts

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
process mergeStringtie {

    label "stringtie"
    label "medium_task"
    publishDir "${params.outdir}/aligned/${name}"

    tag "${name}"

    input:
    set val(name), file("*gtf"), file(gff) from stringtieAssembledTranscripts
        .map { n, rg, f -> [n, f] }
        .groupTuple(by: 0)
        .join( genomes4MergeStringtie.map {n, f, i, g -> [n, g]}, by: 0 )

    output:
    set val(name), file("${name}_stringtie.gtf") into stringtieMergedTranscripts

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
    publishDir "${params.outdir}/assembled"

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
    set val(read_group), file("${read_group}_trinity_denovo.fasta") into assembledTranscripts

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

    mv trinity_assembly/Trinity.fasta "${read_group}_trinity_denovo.fasta"
    rm -rf -- trinity_assembly
    """
}


/*
 * Perform genome guided assembly with Trinity.
process trinityAssembleGuided {

    label "trinity"
    label "big_task"
    publishDir "${params.outdir}/assembled"

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
    set val(name), val(read_group), file("${name}_${read_group}_trinity_guided.fasta") into guidedAssembledTranscripts

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

    mv trinity_assembly/Trinity-GG.fasta "${name}_trinity_guided.fasta"
    rm -rf -- *bam merged.bam.bai trinity_assembly
    """
}
 */


//
// 1 align transcripts to all genomes
//


/*
 * Might be good to keep track of which assemblies come from trinity denovo?
 * Pasa can use this somehow.
 */
process combineTranscripts {

    label "python3"
    label "small_task"

    input:
    file "*fasta" from transcripts
        .mix(assembledTranscripts.map { rg, f -> f } )
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
 * Filters transcripts against univec and extracts poly-A sequences.
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
 * Align transcripts using Spaln
 */
process alignSpalnTranscripts {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/${name}"

    tag { name }

    input:
    set file(fasta),
        file(fasta_cln),
        file(fasta_clean),
        val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq) from transcripts4AlignSpalnTranscripts
            .combine(spalnIndices4AlignSpalnTranscripts)

    output:
    set val(name), file("${name}_spaln_transcripts.gff3") into spalnAlignedTranscripts

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
    > "${name}_spaln_transcripts.gff3"
    """
}


/*
 * Align transcripts using gmap
 */
process alignGmapTranscripts {

    label "gmap"
    label "medium_task"

    publishDir "${params.outdir}/aligned/${name}"

    tag { name }

    input:
    set file(fasta),
        file(fasta_cln),
        file(fasta_clean),
        val(name),
        file("db") from transcripts4AlignGmapTranscripts
            .combine(gmapIndices)

    output:
    set val(name), file("${name}_gmap_transcripts.gff3") into gmapAlignedTranscripts

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
    > ${name}_gmap_transcripts.gff3
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
 * Align all proteins to genome with Spaln
 */
process alignSpalnProteins {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/${name}"

    tag { name }

    input:
    set val(name),
        file(bkn),
        file(ent),
        file(idx),
        file(bkp),
        file(grp),
        file(seq),
        file("proteins.fasta") from spalnIndices4AlignSpalnProteins
            .combine(combinedProteins)

    output:
    set val(name), file("${name}_spaln_proteins.gff3")

    script:
    def species = "Dothideo"

    """
    spaln \
      -L \
      -M3 \
      -O2 \
      -Q7 \
      -T${species} \
      -ya2 \
      -t ${task.cpus} \
      -a \
      -d "${name}" \
      "proteins.fasta" \
    > "${name}_spaln_proteins.gff3"
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
    publishDir "${params.outdir}/annotations/${name}"

    tag { name }

    input:
    set val(name),
        file(genome_fasta),
        file(genome_faidx),
        file(known_sites),
        file(stringtie_gtf),
        file(gmap_aligned),
        file(transcripts_fasta),
        file(transcripts_fasta_cln),
        file(transcripts_fasta_clean) from processed4RunPASA
            .combine(transcripts4RunPasa)

    output:
    set val(name),
        file("${name}_pasa.gff3"),
        file("${name}_pasa_cds.fna"),
        file("${name}_pasa_protein.faa") into pasaPredictions

    script:
    def use_stringent = params.notfungus ? '' : "--stringent_alignment_overlap 30.0 "
    // Don't use stringtie if it is fungus
    def use_stringtie = (stringtie_gtf.name == "WAS_NULL" || !params.notfungus) ? '' : "--trans_gtf ${stringtie_gtf} "
    def use_known = known_sites.name == "WAS_NULL" ? '' : "-L --annots ${known_sites} "

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

    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 \${PWD}/${name}_pasa.gff3
    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.cds \${PWD}/${name}_pasa_cds.fna
    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.pep \${PWD}/${name}_pasa_protein.faa
    """
}


/*
 * Extract hints to be used for augustus and genemark.
 */
process extractAugustusRnaseqHints {

    label "braker"
    label "medium_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name} - ${read_group}"

    input:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file(cram),
        val(strand) from crams4ExtractAugustusRnaseqHints

    output:
    set val(name),
        val(read_group),
        file("${name}_${read_group}_intron_hints.gff3"),
        file("${name}_${read_group}_exon_hints.gff3") into augustusRnaseqHints

    script:
    if (strand == "fr") {
        fst = "forward"
        snd = "reverse"
    } else {
        fst = "reverse"
        snd = "forward"
    }

    """
    # Convert cram to bam.
    # `-F 3328`  excludes these flags
    # not primary alignment (0x100)
    # read is PCR or optical duplicate (0x400)
    # supplementary alignment (0x800)
    samtools view \
        -b \
        -T "${fasta}" \
        -F 3328 \
        -q 25 \
        -@ "${task.cpus}" \
        -o "tmp.bam" \
        "${cram}"

    # Extract introns
    bam2hints \
      --intronsonly \
      --in="tmp.bam" \
      --out="tmp.gff3"

    filterIntronsFindStrand.pl \
      "${fasta}" \
      tmp.gff3 \
      --score \
      > "${name}_${read_group}_intron_hints.gff3"

    rm -f tmp.gff3

    # Extract exons
    bam2hints () {
          bam2wig "\${1}.bam" \
        | wig2hints.pl \
            --width=10 \
            --margin=10 \
            --minthresh=2 \
            --minscore=4 \
            --prune=0.1 \
            --src=W \
            --type=ep \
            --UCSC="\${1}.track" \
            --radius=4.5 \
            --pri=4 \
            --strand="\${2}" \
        > "\${1}_hints.gff3"
    }

    if [ ${fst} == forward ]
    then
        samtools view -b -f 65 -@ "${task.cpus}" "tmp.bam" > "forward.bam"
        bam2hints "forward" "+"
        rm -f "forward.bam"

        samtools view -b -f 128 -@ "${task.cpus}" "tmp.bam" > "reverse.bam"
        bam2hints "reverse" "-"
        rm -f "reverse.bam"
    else
        samtools view -b -f 128 -@ "${task.cpus}" "tmp.bam" > "forward.bam"
        bam2hints "forward" "+"
        rm -f "forward.bam"

        samtools view -b -f 65 -@ "${task.cpus}" "tmp.bam" > "reverse.bam"
        bam2hints "reverse" "-"
        rm -f "reverse.bam"
    fi

    samtools view -b -F 193 -@ "${task.cpus}" "tmp.bam" > "unstranded.bam"
    bam2hints "unstranded" "."
    rm -f "unstranded.bam"

    cat forward_hints.gff3 reverse_hints.gff3 unstranded_hints.gff3 \
      > "${name}_${read_group}_exon_hints.gff3"

    rm -f tmp.bam forward_hints.gff3 reverse_hints.gff3 unstranded_hints.gff3
    """
}

augustusRnaseqHints.set { augustusRnaseqHints4RunGenemark }


/*
 * Do denovo gene prediction with intron hints.
 */
process runGenemark {

    label "genemarkes"
    label "small_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag { name }

    when:
    params.genemark

    input:
    set val(name),
        file(genome),
        file(faidx),
        file("*introns.gff3") from genomes4RunGenemark
            .join(
                augustusRnaseqHints4RunGenemark
                    .map {n, rg, introns, exons -> [n, introns]}
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
 * First pass of coding quarry
 */
process runCodingQuarry {

    label "codingquarry"
    label "medium_task"
    publishDir "${params.outdir}/annotations/${name}"

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
    set val(name),
        file("${name}_codingquarry.gff3"),
        file("${name}_codingquarry_cds.fna"),
        file("${name}_codingquarry_proteins.faa"),
        file("${name}_codingquarry_dubiousset.gff3"),
        file("${name}_codingquarry_fusions.txt"),
        file("${name}_codingquarry_overlapreport.txt") into codingQuarryPredictions

    script:
    """
    grep -v "^#" transcripts.gtf > transcripts.tmp.gtf
    CufflinksGTF_to_CodingQuarryGFF3.py transcripts.tmp.gtf > transcripts.gff3

    CodingQuarry -f genome.fasta -t transcripts.gff3 -p "${task.cpus}"

    \${QUARRY_PATH}/scripts/fastaTranslate.py out/Predicted_CDS.fa \
    | sed 's/*\$//g' > Predicted_Proteins.faa

    \${QUARRY_PATH}/scripts/gene_errors_Xs.py Predicted_Proteins.faa out/Predicted_Proteins.faa
    rm Predicted_Proteins.faa

    mv out/DubiousSet.gff3 "${name}_codingquarry_dubiousset.gff3"
    mv out/PredictedPass.gff3 "${name}_codingquarry.gff3"
    mv out/Predicted_CDS.fa "${name}_codingquarry_cds.fna"
    mv out/Predicted_Proteins.faa "${name}_codingquarry_proteins.faa"
    mv out/fusions.txt "${name}_codingquarry_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarry_overlapreport.txt"

    rm -rf -- out
    """
}




codingQuarryPredictions.into {
    codingQuarryPredictions4SignalP;
    codingQuarryPredictions4PM;
}


/*
 * Predict signal peptides from CodingQuarry first pass predicted proteins.
 */
process getCodingQuarrySignalP {

    label "signalp"
    label "medium_task"

    tag "${name}"

    when:
    !params.notfungus && params.signalp

    input:
    set val(name), file("proteins.faa") from codingQuarryPredictions4SignalP
        .map { n, g, c, p, d, f, o -> [n, p] }

    output:
    set val(name), file("${name}.signalp5") into codingQuarryPredictionsSecreted

    script:
    """
    mkdir tmp
    signalp \
      -fasta "proteins.faa" \
      -prefix "${name}" \
      -org euk \
      -tmp tmp

    mv "${name}_summary.signalp5" "${name}.signalp5"

    rm -rf -- tmp
    """
}


/*
 * CQPM looks for genes that might not be predicted main set because of
 * genome compartmentalisation or differences with signal peptides.
 * NOTE: This fails if there are fewer than 500 secreted genes to train from.
 */
process runCodingQuarryPM {

    label "codingquarry"
    label "medium_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    when:
    !params.notfungus && params.signalp

    input:
    set val(name),
        file("transcripts.gtf"),
        file("genome.fasta"),
        file("genome.fasta.fai"),
        file("${name}_codingquarry.gff3"),
        file("${name}_codingquarry_cds.fna"),
        file("${name}_codingquarry_proteins.faa"),
        file("${name}_codingquarry_dubiousset.gff3"),
        file("${name}_codingquarry_fusions.txt"),
        file("${name}_codingquarry_overlapreport.txt"),
        file("secretome.signalp5") from stringtieMergedTranscripts4CodingQuarryPM
            .combine(genomes4RunCodingQuarryPM, by: 0)
            .combine(codingQuarryPredictions4PM, by: 0)
            .combine(codingQuarryPredictionsSecreted, by: 0)

    output:
    set val(name), file("cqpm") into codingQuarryPMPredictions

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
      -2 "${name}_codingquarry.gff3" \
      -p "${task.cpus}" \
      -g secretome.txt \
      -h

    mv out cqpm
    """
}


/*
 * Run the GeMoMa pipeline
 */
process extractGemomaRnaseqHints {

    label "gemoma"
    label "medium_task"

    tag "${name} - ${read_group}"

    input:
    set val(name),
        val(read_group),
        file(fasta),
        file(faidx),
        file(cram),
        val(strand) from crams4ExtractGemomaRnaseqHints

    output:
    set val(name),
        val(read_group),
        file("${name}_${read_group}_introns.gff"),
        file("${name}_${read_group}_forward.bedgraph"),
        file("${name}_${read_group}_reverse.bedgraph") into gemomaRnaseqHints

    script:
    // FR_UNSTRANDED also valid option
    def strand_flag = strand == "fr" ? "s=FR_SECOND_STRAND " : "s=FR_FIRST_STRAND "

    """
    # Convert cram to bam.
    # `-F 3328`  excludes these flags
    # not primary alignment (0x100)
    # read is PCR or optical duplicate (0x400)
    # supplementary alignment (0x800)
    samtools view \
        -b \
        -T "${fasta}" \
        -F 3328 \
        -q 25 \
        -@ "${task.cpus}" \
        -o "tmp.bam" \
        "${cram}"

    java -jar \${GEMOMA_JAR} CLI ERE ${strand_flag} m=tmp.bam c=true

    # m - mapped reads file (BAM/SAM files containing the mapped reads)	= null
    # u - use secondary alignments (allows to filter flags in the SAM or BAM, default = true)	= true

    mv introns.gff "${name}_${read_group}_introns.gff"
    mv coverage_forward.bedgraph "${name}_${read_group}_forward.bedgraph"
    mv coverage_reverse.bedgraph "${name}_${read_group}_reverse.bedgraph"

    rm -f *.bam protocol_ERE.txt
    rm -rf -- GeMoMa_temp
    """
}


process combineGemomaRnaseqHints {

    label "gemoma"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file("*i.gff"),
        file("*f.bedgraph"),
        file("*r.bedgraph") from gemomaRnaseqHints
            .map { n, rg, i, f, r -> [n, i, f, r] }
            .groupTuple(by: 0)

    output:
    set val(name),
        file("introns.gff"),
        file("forward_coverage.bedgraph"),
        file("reverse_coverage.bedgraph") into combinedGemomaRnaseqHints

    script:
    """
    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineIntronFiles \
      introns.gff \
      *i.gff

    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineCoverageFiles \
      forward_coverage.bedgraph \
      *f.bedgraph

    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineCoverageFiles \
      reverse_coverage.bedgraph \
      *r.bedgraph
    """
}


process extractGemomaCDSParts {

    label "gemoma"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file(gff) from genomes4ExtractGemomaCDSParts
            .filter { n, fa, fi, g -> g.name != 'WAS_NULL' }

    output:
    set val(name),
        file("cds-parts.fasta"),
        file("assignment.tabular"),
        file("proteins.fasta") into gemomaCDSParts

    script:
    """
    java -jar \${GEMOMA_JAR} CLI Extractor \
      a=${gff} \
      g=${fasta} \
      p=true \
      outdir=.

    rm -rf -- GeMoMa_temp
    rm protocol_Extractor.txt
    """
}

gemomaCDSParts.set { gemomaCDSParts4AlignGemomaCDSParts }


process alignGemomaCDSParts {

    label "mmseqs"
    label "medium_task"

    tag "${target_name} - ${ref_name}"

    input:
    set val(ref_name),
        file("cds-parts.fasta"),
        file("assignment.tabular"),
        file("proteins.fasta"),
        val(target_name),
        file(fasta),
        file(faidx) from gemomaCDSParts4AlignGemomaCDSParts
            .combine( genomes4AlignGemomaCDSParts )
            .filter { rn, c, a, p, tn, fa, fi -> rn != tn }

    output:
    set val(target_name),
        val(ref_name),
        file("matches.tsv"),
        file("assignment.tabular"),
        file("proteins.fasta") into alignedGemomaCDSParts

    script:
    // Todo add translation table option using gc option
    """
    mkdir -p genome
    # Stopping splitting by len is important. Otherwise scaffold names don't match.
    mmseqs createdb "${fasta}" genome/db --dont-split-seq-by-len

    mkdir -p proteins
    mmseqs createdb cds-parts.fasta proteins/db

    mkdir -p alignment tmp
    mmseqs search \
      proteins/db \
      genome/db \
      alignment/db \
      tmp \
      --threads ${task.cpus} \
      -e 100 \
      --min-length 10 \
      --comp-bias-corr 1 \
      --split-mode 1 \
      --realign \
      --max-seqs 100 \
      --mask 0 \
      --orf-start-mode 1 \
      --translation-table 1 \
      --use-all-table-starts

    mmseqs convertalis \
      proteins/db \
      genome/db \
      alignment/db \
      matches.tsv \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen'

    rm -rf -- genome proteins alignment tmp
    """
}


/*
*/
process runGemoma {

    label "gemoma"
    label "small_task"

    tag "${target_name} - ${ref_name}"

    input:
    set val(target_name),
        file(fasta),
        file(faidx),
        val(ref_name),
        file("matches.tsv"),
        file("assignment.tabular"),
        file("proteins.fasta"),
        file("introns.gff"),
        file("forward_coverage.bedgraph"),
        file("reverse_coverage.bedgraph") from genomes4RunGemoma
            .join(alignedGemomaCDSParts, by: 0)
            .combine(combinedGemomaRnaseqHints, by: 0)
            .filter { tn, fa, fi, rn, bl, a, p, i, fc, rc -> tn != rn }

    output:
    set val(target_name),
        val(ref_name),
        file("${name}_${ref_name}_preds.gff3") into indivGemomaPredictions

    script:
    // option g= allows genetic code to be provided as some kind of file.
    """
    mkdir -p out
    java -jar \${GEMOMA_JAR} CLI GeMoMa \
      s=matches.tsv \
      t=${fasta} \
      c=cds-parts.fasta \
      a=assignment.tabular \
      q=proteins.fasta \
      outdir=out \
      sort=true \
      i=introns.gff \
      r=2 \
      coverage=STRANDED \
      coverage_forward=coverage_forward.bedgraph \
      coverage_reverse=coverage_reverse.bedgraph

    mv out/predicted_annotation.gff "${name}_${ref_name}_preds.gff3"

    rm -rf -- GeMoMa_temp out
    """
}


process combineGemomaPredictions {

    label "gemoma"
    label "small_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag { name }

    input:
    set val(name),
        val(ref_names),
        file(pred_gffs) from indivGemomaPredictions

    output:
    set val(name), file("${name}_gemoma.gff3") into gemomaPredictions

    script:
    def preds = [ref_names, pred_gffs]
        .transpose()
        .collect { rn, pred -> "p=${rn} g=${pred}" }.join(' ')

    """
    mkdir -p gaf
    java -jar \${GEMOMA_JAR} CLI GAF \
      ${preds} \
      outdir=gaf

    mkdir -p finalised
    java -jar \${GEMOMA_JAR} CLI AnnotationFinalizer \
      g=${target} \
      a=gaf/filtered_predictions.gff \
      u=YES \
      i=introns.gff \
      c=STRANDED \
      coverage_forward=coverage_forward.bedgraph \
      coverage_reverse=coverage_reverse.bedgraph \
      outdir=finalised \
      rename=NO

    mv finalised/final_annotation.gff "${name}_gemoma.gff3"
    rm -rf -- gaf finalised GeMoMa_temp
    """
}

/*
 * Train Augustus models using PASA genes
process trainAugustus {

    script:
    """
    // To find "complete" cdss for training augustus...
    awk '\$0 ~ /^>.*type:complete/ {
      r = gensub(/^>[[:space:]]?([^[:space:]]+).*<REMOVETHIS>/, "\\1", "g", $0);
      print r;
    }' ${name}_cds.fna > complete.orfs
    
    grep -F -f complete.orfs "${name}.gff3" \
    | awk '$3 == "exon" || $3 == "CDS"' \
    | 's/cds\.//; s/\.exon[[:digit:]]*<REMOVETHIS>//' \
    | sort -s -n -k 1,1 -k 4 -k 9 \
    > training_set_complete.gff3
    """
}
 */



// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats
