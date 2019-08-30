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

    nextflow run -profile dev,nimbus main.nf -resume \
      --genomes "genomes/*.fasta" \
      --transcripts 'transcripts/*.fasta' \
      --proteins 'genomes/*.faa' \
      --crams config_crams.tsv \
      --known_sites config.tsv \
      --outdir run \
      --notrinity \
      --signalp \
      --augustus_species parastagonospora_nodorum_sn15 \
      --augustus_config ../augustus_training/07-optimise_utr/augustus_config \
      --genemark \
      --remote_proteins ./uniref-identity%3A0.9.fasta
    ```

    ## Parameters
    See the comments in main for now.
    I'll add documentation here when the interface becomes a bit more stable.

    ## Output
    
    You can specify the main output directory using the `--outdir` parameter.
    Here we use <outdir> inplace of that parameter and <name> inplace of the
    genome basenames.

    <outdir>/aligned/<name>/<name>_* -- Aligned reads, transcripts, and proteins to each genome.
    <outdir>/annotations/<name>/<name>_* -- Gene predictions for each genome, including output of individual methods.
    <outdir>/hints/<name>/<name>_*.gff3 -- GFF files from alignments and annotations configured to be provided as hints to augustus. 
    <outdir>/assembled/<read_group>_* -- Trinity assembled reads.
    <outdir>/qc/<name>_* -- Statistics and QC for each genome and step.

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}


// Target genomes to annotate as softmasked fasta files.
// Note the file basename (filename up to but excluding the
// last extension) is used to match genomes with other hints.
params.genomes = false

// A table specifying known containing the columns name and 
// known_sites specifying a set of high confidence genes to transfer
// across genomes.
params.known_sites = false

// Don't allow pre-hints yet, but we will.
// Probably need options for table style spec and directory
// (e.g. `hints/<name>/hint1.gff`).
// params.hints = false

// A glob of transcripts fasta files to map to the genomes.
// E.g. from prior trininty assembly.
params.transcripts = false

// A glob of protein fasta files from closely related organisms to align
// to genomes and use as hints.
params.proteins = false

// A single fasta file of many proteins from diverse taxonomic sources
// (E.g. uniref90). Used as weaker hints than proteins from more closely
// related organisms.
params.remote_proteins = false

// A MAF file including all species under consideration and may include
// additional ones. Create this using progressiveCactus or sibellia-z.
// If this is not provided, one will be created with sibellia-z.
params.genome_alignment = false

// The busco lineage to use to evaluate gene prediction completeness.
// If this is not provided, busco will not be run.
params.busco_lineage = false


// Augustus params

// The folder containing trained models for augustus [required].
params.augustus_config = false
// The species to use within the augustus_config [required].
params.augustus_species = false

// Run augustus denovo against all genomes.
// The denovo predictions don't form part of the output, this is only
// used to evaluate the performance of the augustus models.
params.augustus_denovo = false

// Run all augustus prediction methods with UTR prediction.
// Note that final combining prediction steps are always run
// with UTR predictions so augustus MUST be trained to use UTRs.
// Use this if the UTR model performs better than the non-UTR model.
params.augustus_utr = false

// The weighting config file for running augustus with hints.
params.augustus_hint_weights = "data/extrinsic_hints.cfg"

// The weighting config file for combining annotations from
// multiple sources.
params.augustus_pred_weights = "data/extrinsic_pred.cfg"


// RNAseq params

// A table specifying RNAseq fastq pairs.
// At the moment RNAseq must be stranded.
// Must contain a "read_group" and "read1" and "read2" column.
// When provided, will assemble read groups using Trinity, and align
// to all genomes using STAR.
params.fastq = false

// A table mapping pre-aligned RNAseq reads to the genomes.
// Should include a "name", and "cram" column.
params.crams = false

// RNAseq is FR stranded rather than RF
// (typical Illumina stranded configuration).
// This can be overwritten for individual fastq pairs or crams
// by adding a "strand" column specifying "fr" or "rf".
// This parameter sets a default for when that column is unspecified.
params.fr = false

// Misc parameters.

// Don't use parameters that are optimised for eukaryotes with
// high-gene density.
params.notfungus = false

// Run GeneMark. This must be explicitly specified because GeneMark requires a license.
params.genemark = false

// Use signalp5 rather than deepsig to train CodingQuarryPM models.
// This is not the default as signalp requires a license.
params.signalp = false

// Don't run trinity even if fastq files are provided.
// Useful if you have precomputed transcripts but not cram files.
params.notrinity = false

// Don't run STAR mapping even if fastq files are provided.
// Useful if you have precomputed cram files but not transcripts.
params.nostar = false

// Don't run codingquarry pathogen mode.
// Use if you're not working with genomes unlikely to have tricky-to-predict
// secreted genes, or if you don't have enough memory to run pathogen mode.
// Note, cqpm is also not run if the `--notfungus` option is used.
params.nocqpm = false

// Specify that we are running a training set to optimise the hints parameters.
// Should be the name of the "genome" that contains the training set.
params.training = false


/*
 * Sanitise the input
 */

def stranded = params.fr ? "fr" : "rf"

def is_null = { f -> (f == null || f == '') }

if ( !params.augustus_species ) {
    log.error "Please nominate one isolate as a reference or provide an Augustus species."
    exit 1
}

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


if ( params.remote_proteins ) {
    Channel
        .fromPath(params.remote_proteins, checkIfExists: true, type: "file")
        .first()
        .set { remoteProteins }

} else {
    remoteProteins = Channel.empty()
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


Channel
    .fromPath(params.augustus_pred_weights, checkIfExists: true, type: "file")
    .first()
    .set { augustusPredWeights }

Channel
    .fromPath(params.augustus_hint_weights, checkIfExists: true, type: "file")
    .first()
    .set { augustusHintWeights }


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

    tag "${name}"

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

genomesWithFaidx
    .tap { genomes4KnownSites }
    .filter  { n, f, i -> (!params.training || params.training == n) }
    .into {
        genomes4SpalnIndex;
        genomes4GmapIndex;
        genomes4MatchRemoteProteinsToGenome;
        genomes4ClusterRemoteProteinsToGenome;
        genomes4AlignRemoteProteinsToGenome;
        genomes4UserCrams;
        genomes4TidyBams;
        genomes4TidyFilteredBams;
        genomes4ExtractSpliceSites;
        genomes4Busco;
        genomes4RunGenemark;
        genomes4TidyGenemark;
        genomes4RunCodingQuarry;
        genomes4RunCodingQuarryPM;
        genomes4AlignGemomaCDSParts;
        genomes4RunGemoma;
        genomes4CombineGemoma;
        genomes4TidyGemoma;
        genomes4ChunkifyGenomes;
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

genomesWithKnownSites
    .tap { genomes4ExtractGemomaCDSParts }
    .filter  { n, f, i, g -> (!params.training || params.training == n) }
    .into {
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
    .into {
        fastq4TrinityAssemble;
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


/*
 * Index the genome for alignment with Spaln.
 */
process getSpalnIndex {

    label "spaln"
    label "small_task"

    tag "${name}"

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
    makeidx.pl -inp "${genome}"
    """
}

spalnIndices.into {
    spalnIndices4AlignSpalnTranscripts;
    spalnIndices4AlignSpalnProteins;
    spalnIndices4AlignSpalnRemoteProteins;
}


/*
 * Index the genomes for GMAP
 */
process getGmapIndex {

    label "gmap"
    label "medium_task"

    tag "${name}"

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
 * Preindex MMSeqs remote proteins
 */
process indexRemoteProteins {

    label "mmseqs"
    label "small_task"

    input:
    file fasta from remoteProteins

    output:
    file "proteins" into indexedRemoteProteins
    file "proteins.tsv" into remoteProteinsTSV

    script:
    """
    mkdir -p tmp proteins

    mmseqs createdb "${fasta}" proteins/db
    mmseqs createindex proteins/db tmp --threads "${task.cpus}"

    awk '
      /^>/ {
        b=gensub(/^>\\s*(\\S+).*\$/, "\\\\1", "g", \$0);
        printf("%s%s\\t", (N>0?"\\n":""), b);
        N++;
        next;
      }
      {
        printf("%s", \$0)
      }
      END {
        printf("\\n");
      }
    ' < "${fasta}" \
    > "proteins.tsv"

    rm -rf -- tmp
    """
}


/*
 * Index genome for STAR
 */
process getStarIndex {

    label "star"
    label "medium_task"

    tag "${name}"

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
 * This finds novel splice sites in the rnaseq.
 * We filter out SS with poor support later.
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

    def min_intron_len = 5
    def max_intron_len = 20000

    def sj_filters = "--outSJfilterReads Unique " +
                     "--outSJfilterOverhangMin 15 6 6 6 " + 
                     "--outSJfilterCountUniqueMin 10 5 5 5 " + 
                     "--outSJfilterDistToOtherSJmin 10 0 3 3 " +
                     "--outSJfilterIntronMaxVsReadN 10 100 500 1000 5000 10000"

    """
    STAR \
      --runThreadN "${task.cpus}" \
      --readFilesCommand zcat \
      --genomeDir "index" \
      --outSAMtype None \
      --outSAMmode None \
      ${sj_filters} \
      --alignIntronMin ${min_intron_len} \
      --alignIntronMax ${max_intron_len} \
      --alignSJoverhangMin 5 \
      --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds No \
      --outFileNamePrefix "${name}_${read_group}." \
      --readFilesIn "${r1_joined}" "${r2_joined}"
    """
}


/*
 * Perform second pass STAR alignment.
 * This uses the predicted splice sites from the previous step
 * but only includes higher confidence splice sites.
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

    def min_intron_len = 5
    def max_intron_len = 20000

    def sj_filters = "--outSJfilterReads All " +
                     "--outSJfilterCountUniqueMin 10 5 5 5 " +
                     "--outSJfilterIntronMaxVsReadN 5 500 5000"

    """
    STAR \
      --runThreadN ${task.cpus} \
      --readFilesCommand zcat \
      --genomeDir "index" \
      --sjdbFileChrStartEnd *SJ.out.tab \
      --outSAMtype BAM Unsorted \
      --outBAMcompression 1 \
      ${sj_filters} \
      --alignIntronMin ${min_intron_len} \
      --alignIntronMax ${max_intron_len} \
      --alignSJoverhangMin 10 \
      --alignSJDBoverhangMin 1 \
      --alignSoftClipAtReferenceEnds No \
      --outFilterType BySJout \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNmax 10 \
      --outFilterMismatchNoverLmax 0.2 \
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
      --outFilterIntronStrands RemoveInconsistentStrands \
      --outMultimapperOrder Random \
      --outSAMattributes All \
      --outSAMstrandField intronMotif \
      --outSAMattrIHstart 0 \
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
 * Assemble transcripts from RNAseq with Stringtie
 * Note that using the known sites option just completely ignores
 * rnaseq info at the loci so it may appear that No UTRs are included
 * if they aren't in the known sites file.
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
            .combine(
                genomes4AssembleStringtie.map { n, f, i, g -> [n, g] },
                by: 0,
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

    tag "${name}"

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
    def species = "phaenodo"
    def min_intron_len = 20
    def max_gene_size = 20000

    """
    spaln \
      -LS \
      -O0 \
      -Q7 \
      -S3 \
      -T${species} \
      -yX \
      -yS \
      -ya1 \
      -XG ${max_gene_size} \
      -yL${min_intron_len} \
      -t ${task.cpus} \
      -d "${name}" \
      "${fasta_clean}" \
    > "${name}_spaln_transcripts.gff3"
    """
}


/*
 * Mostly this is just to add intron features.
 */
process tidySpalnTranscripts {

    label "genometools"
    label "small_task"

    tag "${name}"

    input:
    set val(name), file("transcripts.gff3") from spalnAlignedTranscripts

    output:
    set val(name),
        file("tidied.gff3") into spalnTidiedTranscripts

    script:
    """
    gt gff3 \
      -tidy \
      -retainids \
      -addintrons <(grep -v "^#" transcripts.gff3) \
    > "tidied.gff3"
    """
}


/*
 * Get augustus hints from spaln results.
 */
process extractSpalnTranscriptHints {

    label "python3"
    label "small_task"

    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    input:
    set val(name), file("spaln.gff3") from spalnTidiedTranscripts

    output:
    set val(name),
        file("${name}_spaln_transcript_hints.gff3") into spalnTranscriptHints

    script:
    """
    gff2hints.py \
        --source E \
        --group-level mRNA \
        --priority 3 \
        --exon-trim 6 \
        --intron-trim 0 \
        spaln.gff3 \
    | awk '\$3 != "genicpart"' \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_spaln_transcripts_", \$9);
          print
        }
      ' \
    > "${name}_spaln_transcript_hints.gff3"
    """
}


/*
 * Align transcripts using gmap
 * We run this separately to PASA to control the number of threads PASA uses.
 * If you use BLAT + gmap in pasa it can run 2*cpus threads, which would screw
 * with our provisioning. Most of the PASA time is spent running transdecoder.
 */
process alignGmapTranscripts {

    label "gmap"
    label "medium_task"

    publishDir "${params.outdir}/aligned/${name}"

    tag "${name}"

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
    def min_intronlength = 20
    def max_intronlength_middle = 20000
    def max_intronlength_ends = 10000
    def trim_end_exons = 12
    def microexon_spliceprob = 0.95
    def canonical_mode = 1
    def cross_species = "" // "--cross-species "

    """
    gmap \
      --npaths=0 \
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

/*
 * Deduplicate identical user provided proteins and concat into single file.
 */
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
 * Align all proteins to genome with Spaln.
 * Spaln is good for closely related proteins, but not great for distant ones.
 */
process alignSpalnProteins {
    label "spaln"
    label "medium_task"

    publishDir "${params.outdir}/aligned/${name}"

    tag "${name}"

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
    set val(name), file("${name}_spaln_proteins.gff3") into spalnAlignedProteins

    script:
    def trans_table = 1
    def min_intron_len = 20
    def max_gene_size = 20000

    """
    spaln \
      -C${trans_table} \
      -KP \
      -LS \
      -M3 \
      -O0 \
      -Q7 \
      -ya1 \
      -yX \
      -yL${min_intron_len} \
      -XG${max_gene_size} \
      -t ${task.cpus} \
      -d "${name}" \
      "proteins.fasta" \
    > "${name}_spaln_proteins.gff3"
    """
}


/*
 * Get hints for augustus from spaln protein alignments.
 */
process extractSpalnProteinHints {

    label "braker"
    label "small_task"

    tag "${name}"

    publishDir "${params.outdir}/hints/${name}"

    input:
    set val(name), file("spaln.gff3") from spalnAlignedProteins

    output:
    set val(name), file("${name}_spaln_protein_hints.gff3") into spalnProteinHints

    script:
    def min_intron_len = 20

    """
    align2hints.pl \
      --in=spaln.gff3 \
      --out=hints.gff3 \
      --prg=spaln \
      --CDSpart_cutoff=12 \
      --minintronlen="${min_intron_len}" \
      --priority=3

    awk 'BEGIN {OFS="\\t"} {sub(/grp=/, "grp=${name}_spaln_proteins_", \$9); print}' \
      hints.gff3 \
    > "${name}_spaln_protein_hints.gff3"
    """
}


/*
 * Quickly find genomic regions with matches to remote proteins.
 * This is an approximate method which we refine later.
 */
process matchRemoteProteinsToGenome {

    label "mmseqs"
    label "big_task"

    tag "${name}"

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file("proteins") from genomes4MatchRemoteProteinsToGenome
            .combine(indexedRemoteProteins)

    output:
    set val(name), file("${name}_remote_proteins.tsv") into matchedRemoteProteinsToGenome

    script:
    """
    mkdir genome result tmp

    # the "dont-split" bit is important for keeping the ids correct.
    mmseqs createdb \
      "${fasta}" \
      genome/db \
      --dont-split-seq-by-len

    # Searching with genome as query is ~3X faster
    mmseqs search \
      genome/db \
      proteins/db \
      result/db \
      tmp \
      --threads "${task.cpus}" \
      -e 0.00001 \
      --min-length 10 \
      --comp-bias-corr 1 \
      --split-mode 1 \
      --max-seqs 50 \
      --mask 0 \
      --orf-start-mode 1 \
      --translation-table 1 \
      --use-all-table-starts

    # Extract match results.
    mmseqs convertalis \
      genome/db \
      proteins/db \
      result/db \
      results_unsorted.tsv \
      --threads "${task.cpus}" \
      --format-mode 0 \
      --format-output "query,target,qstart,qend,qlen,tstart,tend,tlen,alnlen,pident,mismatch,gapopen,evalue,bits"

    sort \
      -k1,1 \
      -k3,3n \
      -k4,4n \
      -k2,2 \
      --parallel="${task.cpus}" \
      --temporary-directory=tmp \
      results_unsorted.tsv \
    > "${name}_remote_proteins.tsv"

    sed -i '1i query\ttarget\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tpident\tmismatch\tgapopen\tevalue\tbitscore' "${name}_remote_proteins.tsv"

    rm -rf -- tmp genome result results_unsorted.tsv
    """
}


/*
 * Finds regions genome with high density of remote protein matches.
 * We align remote proteins to each of these regions individually.
 */
process clusterRemoteProteinsToGenome {

    label "bedtools"
    label "medium_task"

    tag "${name}"

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file("results.tsv") from genomes4ClusterRemoteProteinsToGenome
            .combine(matchedRemoteProteinsToGenome, by: 0)

    output:
    set val(name), file("clustered.bed") into clusteredRemoteProteinsToGenome

    script:
    // We extend the region to align against by this many basepairs.
    // This should be a bit longer than your absolute maximum expected gene length.
    def exonerate_buffer_region = 20000

    """
    mkdir -p tmp

      tail -n+2 results.tsv \
    | awk '
        BEGIN { OFS="\t" }
        \$3 > \$4 { print \$1, \$4, \$3, \$2 }
        \$3 < \$4 { print \$1, \$3, \$4, \$2 }
      ' \
    | sort \
        -k1,1 -k2,2n -k3,3n \
        --temporary-directory=tmp \
    | sed 's/,/%2C/g' \
    | bedtools merge -d 1000 -c 4 -o distinct -i - \
    | bedtools slop -g "${faidx}" -b "${exonerate_buffer_region}" -i - \
    > clustered.bed

    rm -rf -- tmp
    """
}


/*
 * This aligns the proteins identified in the "match" step to
 * the genomic regions that they matched. MMseqs is much faster at identifying
 * regions but is less accurate and can't model introns. 
 * Exonerate seems to be better for more remote proteins than spaln.
 * Spaln introduces lots of very short CDSs, i think they're frameshifts.
 */
process alignRemoteProteinsToGenome {

    label "exonerate"
    label "big_task"

    tag "${name}"

    input:
    set val(name),
        file(fasta),
        file(faidx),
        file("clustered.bed"),
        file("proteins.tsv") from genomes4AlignRemoteProteinsToGenome
            .combine(clusteredRemoteProteinsToGenome, by: 0)
            .combine(remoteProteinsTSV)

    output:
    set val(name),
        file("${name}_remote_proteins_exonerate.gff") into alignedRemoteProteinsToGenome

    script:
    """
    mkdir -p tmp

    exonerate_parallel.sh \
      -g "${fasta}" \
      -q "proteins.tsv" \
      -b "clustered.bed" \
      -n "${task.cpus}" \
      -t "tmp" \
      -o "${name}_remote_proteins_exonerate.gff"
    """
}


/*
 * Convert exonerate alignments to augustus hints
 */
process extractExonerateRemoteProteinHints {

    label "braker"
    label "small_task"

    tag "${name}"

    publishDir "${params.outdir}/hints/${name}"

    input:
    set val(name), file("exonerate.gff") from alignedRemoteProteinsToGenome

    output:
    set val(name),
        file("${name}_exonerate_remote_protein_hints.gff3") into exonerateRemoteProteinHints

    script:
    """
    align2hints.pl \
      --in=exonerate.gff \
      --out=hints.gff3 \
      --prg=exonerate \
      --CDSpart_cutoff=15 \
      --minintronlen=20 \
      --priority=2 \
      --source=T

    awk '
      BEGIN {OFS="\\t"}
      \$3 == "CDSpart" {
        sub(/grp=/, "grp=${name}_exonerate_remote_proteins_", \$9)
        print
      }
      ' \
      hints.gff3 \
    > "${name}_exonerate_remote_protein_hints.gff3"
    """
}


//
// 5 Run genemark, pasa, braker2, codingquarry on all genomes using previous steps.
//

/*
 * Add the stringtie assembled and gmap aligned transcripts to the pasa input channel.
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
 * Predict genes using pasa and transdecoder.
 */
process runPASA {

    label "pasa"
    label "medium_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

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
    // Stringtie and cufflinks tend to merge overlapping features,
    // which doesn't work well for organisms with high gene density.
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

pasaPredictions.into {
    pasaPredictions4Tidy;
    pasaPredictions4Busco;
}


/*
 * Add additional features to pasa prediction (e.g. start, stop, intron) etc.
 */
process tidyPasa {

    label "aegean"
    label "small_task"

    tag "${name}"

    input:
    set val(name), file("pasa.gff3") from pasaPredictions4Tidy
        .map { n, g, c, p -> [n, g] }

    output:
    set val(name), file("pasa_tidy.gff3") into tidiedPasa

    script:
    """
    gt gff3 -tidy -sort -retainids pasa.gff3 \
    | canon-gff3 -i - \
    > pasa_tidy.gff3
    """
}

tidiedPasa.into {
    pasaPredictions4Hints;
    pasaPredictions4Stats;
}


/*
 * Get hints for augustus from pasa predictions.
 * We fetch separate hints for PASA and transdecoder predictions.
 */
process extractPasaHints {

    label "python3"
    label "small_task"

    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    input:
    set val(name), file("pasa.gff3") from pasaPredictions4Hints

    output:
    set val(name),
        file("${name}_pasa_hints.gff3"),
        file("${name}_transdecoder_hints.gff3") into pasaHints

    script:
    """
    awk '\$3 == "exon" || \$3 == "intron" || \$3 == "mRNA"' pasa.gff3 \
    | gff2hints.py \
        --source PR \
        --group-level mRNA \
        --priority 4 \
        --exon-trim 9 \
        --intron-trim 0 \
        - \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_pasa_", \$9);
          print
        }
      ' \
    > "${name}_pasa_hints.gff3"

    awk '\$3 != "exon" && \$3 != "intron"' pasa.gff3 \
    | gff2hints.py \
        --source PR \
        --group-level mRNA \
        --priority 4 \
        --cds-trim 9 \
        --utr-trim 6 \
        - \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_transdecoder_", \$9);
          print
        }
      ' \
    > "${name}_transdecoder_hints.gff3"
    """
}


/*
 * Extract hints to be used for augustus and genemark.
 *
 * We only use the intron hints because we have spaln/pasa alignments,
 * and the coverage info causes lots of close genes to be merged or extended
 * even if you use the UTR model.
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
        file("${name}_${read_group}_intron_hints.gff3") into augustusRnaseqHints

    script:
    def min_coverage = 4
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
      --minintronlen=20 \
      --maxcoverage=100 \
      --priority 4 \
      --in="tmp.bam" \
      --out="tmp.gff3"

    filterIntronsFindStrand.pl \
        "${fasta}" \
        tmp.gff3 \
        --allowed=gtag,gcag,atac,ctac,gaag \
        --score \
    | awk '\$6 >= ${min_coverage}' \
    > "${name}_${read_group}_intron_hints.gff3"

    rm -f tmp.gff3
    """
}

augustusRnaseqHints.into {
    augustusRnaseqHints4RunGenemark;
    augustusRnaseqHints4JoinHints;
}


/*
 * Do denovo gene prediction with intron hint training.
 */
process runGenemark {

    label "genemarkes"
    label "medium_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    when:
    params.genemark

    input:
    set val(name),
        file(genome),
        file(faidx),
        file("*introns.gff3") from genomes4RunGenemark
            .join(
                augustusRnaseqHints4RunGenemark
                    .map {n, rg, introns -> [n, introns]}
                    .groupTuple(by: 0),
                by: 0
            )

    output:
    set val(name), file("${name}_genemark.gtf") into genemarkPredictions

    script:
    def use_fungus = params.notfungus ? '' : '--fungus '
    def is_training = params.training ? '--min_contig 300 ': ''

    """
    sort -k1,1V -k4,4n -k5,5rn -k3,3r *introns.gff3 > hints.gff3

    gmes_petap.pl \
      --cores "${task.cpus}" \
      --soft_mask 100 \
      --ET "hints.gff3" \
      ${use_fungus} \
      ${is_training} \
      --sequence "${genome}" \
      --evidence "hints.gff3" \

    mv genemark.gtf "${name}_genemark.gtf"
    """
}


/*
 * Convert genemark gtf to gff. Add introns etc. Extract protiens.
 */
process tidyGenemark {

    label "aegean"
    label "small_task"

    tag "${name}"

    when:
    params.genemark

    input:
    set val(name),
        file("genemark.gtf"),
        file(fasta),
        file(faidx) from genemarkPredictions
            .join(genomes4TidyGenemark, by: 0)

    output:
    set val(name), file("genemark.gff3") into tidiedGenemark
    set val(name),
        val("genemark"),
        file("genemark.faa") into genemarkPredictions4Busco

    script:
    def genetic_code = 1

    """
    gt gtf_to_gff3 -tidy genemark.gtf \
    | gt gff3 -tidy -sort -retainids \
    | canon-gff3 -i - \
    > genemark.gff3

    gt extractfeat \
      -type CDS \
      -join \
      -translate \
      -retainids \
      -gcode "${genetic_code}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      genemark.gff3 \
    > genemark.faa
    """
}

tidiedGenemark.into {
    genemarkPredictions4Hints;
    genemarkPredictions4Stats;
}


/*
 * Get hints for augustus.
 */
process extractGenemarkHints {

    label "python3"
    label "small_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    when:
    params.genemark

    input:
    set val(name), file("genemark.gff3") from genemarkPredictions4Hints

    output:
    set val(name), file("${name}_genemark_hints.gff3") into genemarkHints

    script:
    """
    gff2hints.py \
        --source PR \
        --group-level mRNA \
        --priority 3 \
        --cds-trim 9 \
        --exon-trim 9 \
        --intron-trim 0 \
        genemark.gff3 \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_genemark_", \$9);
          print
        }
      ' \
    > "${name}_genemark_hints.gff3"
    """
}


/*
 * Predict genes using codingquarry
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
    | sed 's/*\$//g' \
    > Predicted_Proteins.faa

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
    codingQuarryPredictions4SecretionPred;
    codingQuarryPredictions4PM;
    codingQuarryPredictions4Tidy;
    codingQuarryPredictions4Busco;
}


/*
 * Fix weird cq gffs and get introns.
 */
process tidyCodingQuarry {

    label "aegean"
    label "small_task"

    tag "${name}"

    when:
    !params.notfungus

    input:
    set val(name),
        file("codingquarry.gff3") from codingQuarryPredictions4Tidy
            .map { n, g, c, p, d, f, o -> [n, g] }

    output:
    set val(name), file("codingquarry_tidy.gff3") into tidiedCodingQuarry

    script:
    """
    awk 'BEGIN {OFS="\\t"} \$8 = "-1" {\$8="0"} {print}' codingquarry.gff3 \
    | awk 'BEGIN {OFS="\\t"} \$3 == "gene" {\$3="mRNA"} {print}' \
    | gt gff3 -tidy -sort -retainids \
    | canon-gff3 -i - \
    > codingquarry_tidy.gff3
    """
}

tidiedCodingQuarry.into {
    codingQuarryPredictions4Hints;
    codingQuarryPredictions4Stats;
}


/*
 * Get hints for augustus
 */
process extractCodingQuarryHints {

    label "python3"
    label "small_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    when:
    !params.notfungus

    input:
    set val(name), file("codingquarry.gff3") from codingQuarryPredictions4Hints

    output:
    set val(name), file("${name}_codingquarry_hints.gff3") into codingQuarryHints

    script:
    """
    gff2hints.py \
        --source PR \
        -g mRNA \
        --priority 4 \
        --cds-trim 6 \
        --exon-trim 6 \
        --intron-trim 0 \
        -- \
        codingquarry.gff3 \
    | awk 'BEGIN {OFS="\\t"} {sub(/group=/, "group=${name}_codingquarry_", \$9); print}' \
    > "${name}_codingquarry_hints.gff3"
    """
}


/*
 * Predict signal peptides from CodingQuarry first pass predicted proteins.
 */
if (params.signalp) {

    process getCodingQuarrySignalP {

        label "signalp"
        label "medium_task"

        tag "${name}"

        when:
        !params.notfungus

        input:
        set val(name), file("proteins.faa") from codingQuarryPredictions4SecretionPred
            .map { n, g, c, p, d, f, o -> [n, p] }

        output:
        set val(name), file("${name}_summary.signalp5") into codingQuarryPredictionsSecreted

        script:
        """
        mkdir tmp
        signalp \
          -fasta "proteins.faa" \
          -prefix "${name}" \
          -org euk \
          -tmp tmp

        rm -rf -- tmp
        """
    }

    process tidyCodingQuarrySignalp {

        label "posix"
        label "small_task"
        publishDir "${params.outdir}/annotations/${name}"

        tag "${name}"

        when:
        !params.notfungus || !params.nocqpm

        input:
        set val(name), file("secreted.txt") from codingQuarryPredictionsSecreted

        output:
        set val(name), file("${name}_codingquarry_proteins_secreted.txt") into codingQuarryPredictionsSecretedTidy

        script:
        """
        gawk '
          BEGIN {
            OFS=" "
          }
          \$2 ~ /^SP/ {
            match(\$0, /CS pos: ([0-9]*)-/, x)
            print \$1, x[1]
          }
        ' < "secreted.txt" > "${name}_codingquarry_proteins_secreted.txt"
        """
    }

} else {

    process getCodingQuarryDeepsig {

        label "deepsig"
        label "medium_task"
        publishDir "${params.outdir}/annotations/${name}"

        tag "${name}"

        when:
        !params.notfungus

        input:
        set val(name), file("proteins.faa") from codingQuarryPredictions4SecretionPred
            .map { n, g, c, p, d, f, o -> [n, p] }

        output:
        set val(name), file("secreted.txt") into codingQuarryPredictionsSecreted

        script:
        """
        deepsig.py \
          -f "proteins.faa" \
          -o secreted.txt \
          -k euk
        """
    }

    process tidyCodingQuarryDeepsig {

        label "posix"
        label "small_task"
        publishDir "${params.outdir}/annotations/${name}"

        tag "${name}"

        when:
        !params.notfungus || !params.nocqpm

        input:
        set val(name), file("secreted.txt") from codingQuarryPredictionsSecreted

        output:
        set val(name),
            file("${name}_codingquarry_proteins_secreted.txt") into codingQuarryPredictionsSecretedTidy

        script:
        """
        gawk '
          BEGIN {
            OFS=" "
          }
          \$2 == "SignalPeptide" {
            print \$1, \$4
          }
        ' < secreted.txt > "${name}_codingquarry_proteins_secreted.txt"
        """
    }
}


/*
 * CQPM looks for genes that might not be predicted main set because of
 * genome compartmentalisation or differences with signal peptides.
 * NOTE: This fails if there are fewer than 500 secreted genes to train from.
 *
 * CQPM also uses a lot of memory at a specific point, which can cause
 * random segfaults. We retry a few times, but if it keeps failing you
 * probably need to increase the available RAM.
 * Try setting `vm.overcommit_memory = 1` if you get segfaults before
 * reaching max memory.
 */
process runCodingQuarryPM {

    label "codingquarry"
    label "bigmem_task"

    errorStrategy "retry"
    maxRetries 10

    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    when:
    !params.notfungus || !params.nocqpm

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
        file("secretome.txt") from stringtieMergedTranscripts4CodingQuarryPM
            .combine(genomes4RunCodingQuarryPM, by: 0)
            .combine(codingQuarryPredictions4PM, by: 0)
            .combine(codingQuarryPredictionsSecretedTidy, by: 0)

    output:
    set val(name),
        file("${name}_codingquarrypm.gff3"),
        file("${name}_codingquarrypm_cds.fna"),
        file("${name}_codingquarrypm_proteins.faa") optional true into codingQuarryPMPredictions

    script:
    """
    mkdir -p ParameterFiles/RNA_secreted

    NSECRETED=\$(wc -l < secretome.txt)
    if [ \${NSECRETED} -lt 501 ]
    then
        exit 0
    fi

    grep -v "^#" transcripts.gtf > transcripts.tmp.gtf
    CufflinksGTF_to_CodingQuarryGFF3.py transcripts.tmp.gtf > transcripts.gff3

    CodingQuarry \
      -f genome.fasta \
      -t transcripts.gff3 \
      -2 "${name}_codingquarry.gff3" \
      -p "${task.cpus}" \
      -g secretome.txt \
      -h

    \${QUARRY_PATH}/scripts/fastaTranslate.py out/PGN_predicted_CDS.fa \
    | sed 's/*\$//g' \
    > PGN_predicted_Proteins.faa

    \${QUARRY_PATH}/scripts/gene_errors_Xs.py PGN_predicted_Proteins.faa out/PGN_predicted_Proteins.faa
    rm PGN_predicted_Proteins.faa

    mv out/PGN_predictedPass.gff3 "${name}_codingquarrypm.gff3"
    mv out/PGN_predicted_CDS.fa "${name}_codingquarrypm_cds.fna"
    mv out/PGN_predicted_Proteins.faa "${name}_codingquarrypm_proteins.faa"
    mv out/fusions.txt "${name}_codingquarrypm_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarrypm_overlapreport.txt"
    # rm -rf -- out
    """
}

codingQuarryPMPredictions.into {
    codingQuarryPMPredictions4Tidy;
    codingQuarryPMPredictions4Busco;
}


/*
 * Fix weird CodingQuarry gffs and add introns etc.
 */
process tidyCodingQuarryPM {

    label "aegean"
    label "small_task"

    tag "${name}"

    when:
    !params.notfungus

    input:
    set val(name),
        file("codingquarry.gff3") from codingQuarryPMPredictions4Tidy
            .map { n, g, c, p -> [n, g] }

    output:
    set val(name), file("codingquarrypm_tidy.gff3") into tidiedCodingQuarryPM

    script:
    """
      awk 'BEGIN {OFS="\\t"} \$8 == "-1" {\$8="0"} {print}' codingquarry.gff3 \
    | awk 'BEGIN {OFS="\\t"} \$3 == "gene" {\$3="mRNA"} {print}' \
    | gt gff3 -tidy -sort -retainids \
    | canon-gff3 -i - \
    > codingquarrypm_tidy.gff3
    """
}

tidiedCodingQuarryPM.into {
    codingQuarryPMPredictions4Hints;
    codingQuarryPMPredictions4Stats;
}


/*
 * Get hints for augustus
 */
process extractCodingQuarryPMHints {

    label "python3"
    label "small_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    when:
    !params.notfungus

    input:
    set val(name), file("codingquarry.gff3") from codingQuarryPMPredictions4Hints

    output:
    set val(name), file("${name}_codingquarrypm_hints.gff3") into codingQuarryPMHints

    script:
    """
      gff2hints.py \
        --source PR \
        -g mRNA \
        --priority 4 \
        --cds-trim 6 \
        --exon-trim 6 \
        --intron-trim 0 \
        codingquarry.gff3 \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_codingquarrypm_", \$9);
          print
        }
      ' \
    > "${name}_codingquarrypm_hints.gff3"
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


/*
 */
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

combinedGemomaRnaseqHints.into {
    combinedGemomaRnaseqHints4Run;
    combinedGemomaRnaseqHints4Combine;
}


/*
 * Gemoma needs to extract the splice site info as well as the proteins.
 */
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


/*
 * Align proteins to genomes for Gemoma.
 * Do this separately as the gemoma pipeline
 * currently crashes and we can control parallelism better.
 */
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
        file("cds-parts.fasta"),
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
 * Predict genes with gemoma for each known-sites set.
 * Merges adjacent MMseqs matches and checks intron-exon boundaries.
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
        file("cds-parts.fasta"),
        file("matches.tsv"),
        file("assignment.tabular"),
        file("proteins.fasta"),
        file("introns.gff"),
        file("coverage_forward.bedgraph"),
        file("coverage_reverse.bedgraph") from genomes4RunGemoma
            .combine(alignedGemomaCDSParts, by: 0)
            .combine(combinedGemomaRnaseqHints4Run, by: 0)
            .filter { tn, fa, fi, rn, cp, bl, a, p, i, fc, rc -> tn != rn }

    output:
    set val(target_name),
        val(ref_name),
        file("${target_name}_${ref_name}_preds.gff3") into indivGemomaPredictions

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

    mv out/predicted_annotation.gff "${target_name}_${ref_name}_preds.gff3"

    rm -rf -- GeMoMa_temp out
    """
}


/*
 * Combines gemoma predictions from multiple known-sites sets.
 * Also adds UTRs etc to gff based on RNAseq.
 */
process combineGemomaPredictions {

    label "gemoma"
    label "small_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    input:
    set val(name),
        val(ref_names),
        file(pred_gffs),
        file(fasta),
        file(faidx),
        file("introns.gff"),
        file("coverage_forward.bedgraph"),
        file("coverage_reverse.bedgraph") from indivGemomaPredictions
            .groupTuple(by: 0)
            .join(genomes4CombineGemoma, by: 0, remainder: false)
            .combine(combinedGemomaRnaseqHints4Combine, by: 0)

    output:
    set val(name), file("${name}_gemoma.gff3") into gemomaPredictions

    script:
    def ref_names_list = ref_names
    def pred_gffs_list = (pred_gffs instanceof List) ? pred_gffs : [pred_gffs]
    assert pred_gffs_list.size() == ref_names_list.size()

   // The format for GAF is a bit weird so do it here.
   // transpose is like zip() and collect is like map.
    def preds = [ref_names_list, pred_gffs_list]
        .transpose()
        .collect { rn, pred -> "p=${rn} g=${pred.name}" }
        .join(' ')

    """
    mkdir -p gaf
    java -jar \${GEMOMA_JAR} CLI GAF \
      ${preds} \
      outdir=gaf

    mkdir -p finalised
    java -jar \${GEMOMA_JAR} CLI AnnotationFinalizer \
      g=${fasta} \
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
 * Fix missing mRNA feature in output and add intron features etc.
 * extract proteins for busco.
 */
process tidyGemoma {

    label "aegean"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file("gemoma.gff3"),
        file(fasta),
        file(faidx) from gemomaPredictions
            .join(genomes4TidyGemoma, by: 0)

    output:
    set val(name), file("gemoma_tidy.gff3") into tidiedGemoma
    set val(name),
        val("gemoma"),
        file("gemoma.faa") into gemomaPredictions4Busco

    script:
    def genetic_code = 1
    """
      gt gff3 -tidy -sort -retainids gemoma.gff3 \
    | awk 'BEGIN {OFS="\\t"} \$3 == "prediction" {\$3="mRNA"} {print}' \
    | canon-gff3 -i - \
    > gemoma_tidy.gff3

    gt extractfeat \
      -type CDS \
      -join \
      -translate \
      -retainids \
      -gcode "${genetic_code}" \
      -matchdescstart \
      -seqfile "${fasta}" \
      gemoma_tidy.gff3 \
    > gemoma.faa
    """
}


/*
 * Get hints for augustus
 */
process extractGemomaHints {

    label "python3"
    label "small_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    input:
    set val(name), file("gemoma.gff3") from tidiedGemoma

    output:
    set val(name), file("${name}_gemoma_hints.gff3") into gemomaHints

    script:
    """
      gff2hints.py \
        --source PR \
        --group-level mRNA \
        --priority 4 \
        --cds-trim 6 \
        --exon-trim 6 \
        --utr-trim 9 \
        --intron-trim 0 \
        gemoma.gff3 \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_gemoma_", \$9);
          print
        }
      ' \
    > "${name}_gemoma_hints.gff3"
    """
}


/*
 *
 */
process chunkifyGenomes {

    label "python3"
    label "small_task"

    tag "${name}"

    input:
    set val(name),
        file("input.fasta"),
        file("input.fasta.fai") from genomes4ChunkifyGenomes

    output:
    set val(name), file("out_*.fasta") into chunkifiedGenomes

    script:
    """
    chunk_genomes.py -n 16 --prefix "out_" input.fasta
    """
}

chunkifiedGenomes
    .flatMap { n, fs -> fs.collect {f -> [n, f]} }
    .into {
        genomes4RunAugustusDenovo;
        genomes4RunAugustusDenovoUTR;
        genomes4RunAugustusHints;
        genomes4RunAugustusPreds;
    }


/*
 */
process runAugustusDenovo {

    label "augustus"
    label "small_task"

    tag "${name}"

    when:
    params.augustus_denovo

    input:
    set val(name), file(fasta) from genomes4RunAugustusDenovo
    file "augustus_config" from augustusConfig

    output:
    set val(name),
        val("augustus_denovo"),
        val("both"),
        file("out.gff") into augustusDenovoResults

    script:
    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    augustus \
      --species="${params.augustus_species}" \
      --softmasking=on \
      --singlestrand=true \
      --min_intron_len=20 \
      --start=on \
      --stop=on \
      --introns=on \
      --cds=on \
      --gff3=on \
      --UTR=off \
      --codingseq=on \
      --protein=on \
      --outfile="out.gff" \
      --errfile=augustus.err \
      "${fasta}"
    """
}


/*
 */
process runAugustusDenovoUTR {

    label "augustus"
    label "small_task"

    tag "${name} - ${strand}"

    when:
    params.augustus_denovo && params.augustus_utr

    input:
    set val(name),
        val(strand),
        file(fasta) from genomes4RunAugustusDenovoUTR
            .flatMap { n, f -> [[n, "forward", f], [n, "reverse", f]] }

    file "augustus_config" from augustusConfig

    output:
    set val(name),
        val("augustus_denovo_utr"),
        val(strand),
        file("out.gff") into augustusDenovoUTRResults

    script:
    if ( strand == "forward" ) {
        strand_param = "--strand=forward"
    } else if ( strand == "reverse" ) {
        strand_param = "--strand=backward"
    } else {
        log.error "This shouldn't happen"
        exit 1
    }
    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    augustus \
      --species="${params.augustus_species}" \
      --softmasking=on \
      ${strand_param} \
      --UTR=on \
      --min_intron_len=20 \
      --start=on \
      --stop=on \
      --introns=on \
      --cds=on \
      --gff3=on \
      --outfile="out.gff" \
      --errfile=augustus.err \
      "${fasta}"
    """
}

augustusRnaseqHints4JoinHints
    .map { n, rg, i -> [n, i] }
    .mix(spalnTranscriptHints, spalnProteinHints, exonerateRemoteProteinHints)
    .into {
        augustusExtrinsicHints4Hints;
        augustusExtrinsicHints4Preds;
    }


if ( params.augustus_utr ) {

    process filterHintStrandUTR {

        label "posix"
        label "small_task"
        publishDir "${params.outdir}/hints/${name}"

        tag "${name}"

        input:
        set val(name), file("*hints") from augustusExtrinsicHints4Hints
                .groupTuple(by: 0)

        output:
        set val(name),
            file("${name}_pos_hints_for_augustus_hints.gff"),
            file("${name}_neg_hints_for_augustus_hints.gff") into augustusFilteredHints4Hints

        script:
        """
        mkdir tmp

          cat *hints \
        | gawk '
            \$7 == "-" && (\$9 ~ /group/ || \$9 ~ /grp/) {
              b=gensub(/.*gr(ou)?p=([^;]+).*/, "\\\\2", "g", \$9);
              print b;
            }
          ' \
        | sort -u -T tmp \
        > neg_ids.txt

          cat *hints \
        | gawk '
            \$7 == "+" && (\$9 ~ /group/ || \$9 ~ /grp/) {
              b=gensub(/.*gr(ou)?p=([^;]+).*/, "\\\\2", "g", \$9);
              print b;
            }
          ' \
        | sort -u -T tmp \
        > pos_ids.txt
        rm -rf -- tmp

        cat *hints | grep -f neg_ids.txt -F > neg_groups.gff
        cat *hints | grep -f pos_ids.txt -F > pos_groups.gff

          cat *hints \
        | gawk '(\$7 == "-" || \$7 == ".") && !(\$9 ~ /group/ || \$9 ~ /grp/)' \
        > neg_single.gff

          cat *hints \
        | gawk '(\$7 == "+" || \$7 == ".") && !(\$9 ~ /group/ || \$9 ~ /grp/)' \
        > pos_single.gff

        cat neg_single.gff neg_groups.gff > "${name}_neg_hints_for_augustus_hints.gff"
        cat pos_single.gff pos_groups.gff > "${name}_pos_hints_for_augustus_hints.gff"
        """
    }

    augustusFilteredHints4Hints
        .flatMap { n, f, r -> [[n, "forward", f], [n, "reverse", r]]}
        .set { augustusFilteredHints4Hints }

} else {

    process filterHintStrand {

        label "posix"
        label "small_task"
        publishDir "${params.outdir}/hints/${name}"

        tag "${name}"

        input:
        set val(name),
            file("*hints") from augustusExtrinsicHints4Hints
                .groupTuple(by: 0)

        output:
        set val(name),
            val("both"),
            file("${name}_hints_for_augustus_hints.gff") into augustusFilteredHints4Hints

        script:
        """
          cat *hints \
        | awk '\$3 != "exonpart" && \$3 != "exon"' \
        > "${name}_hints_for_augustus_hints.gff"
        """
    }
}


/*
 */
process runAugustusHints {

    label "augustus"
    label "small_task"

    tag "${name} - ${strand}"

    input:
    set val(name),
        file(fasta),
        val(strand),
        file("hints.gff") from genomes4RunAugustusHints
            .combine(augustusFilteredHints4Hints, by: 0)

    file "augustus_config" from augustusConfig
    file "extrinsic.cfg" from augustusHintWeights

    output:
    set val(name),
        val("augustus_hints"),
        val(strand),
        file("out.gff") into augustusHintsResults

    script:
    if ( !params.augustus_utr ) {
        strand_param = "--singlestrand=true --UTR=off"
    } else if ( strand == "forward" ) {
        strand_param = "--strand=forward --UTR=on"
    } else if ( strand == "reverse" ) {
        strand_param = "--strand=backward --UTR=on"
    }

    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    perl -n -e'/>(\\S+)/ && print \$1."\\n"' < "${fasta}" > seqids.txt

    getLinesMatching.pl seqids.txt 1 < hints.gff > hints_filtered.gff

    augustus \
      --species="${params.augustus_species}" \
      --extrinsicCfgFile=extrinsic.cfg \
      --hintsfile=hints_filtered.gff \
      ${strand_param} \
      --allow_hinted_splicesites=atac \
      --softmasking=on \
      --alternatives-from-evidence=true \
      --min_intron_len=20 \
      --start=on \
      --stop=on \
      --introns=on \
      --cds=on \
      --gff3=on \
      --outfile="out.gff" \
      --errfile=augustus.err \
      "${fasta}"
    """
}

augustusDenovoResults
    .mix(
        augustusDenovoUTRResults,
        augustusHintsResults,
    )
    .map { n, p, s, g -> [n, p, g] }
    .groupTuple(by: [0, 1])
    .set {augustusChunks}


/*
 */
process joinAugustusChunks {

    label "genometools"
    label "small_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name} - ${paramset}"

    input:
    set val(name),
        val(paramset),
        file("*chunks.gff") from augustusChunks

    output:
    set val(name),
        val(paramset),
        file("${name}_${paramset}.gff3") into augustusJoinedChunks

    script:
    """
    for f in *chunks.gff
    do
        awk 'BEGIN {OFS="\t"} \$3 == "transcript" {\$3="mRNA"} {print}' \${f} \
      | gt gff3 -tidy -sort -o \${f}_tidied.gff3 -
    done

    gt merge -tidy -o "${name}_${paramset}.gff3" *_tidied.gff3
    """
}


/*
 */
process extractAugustusHintsHints {

    label "python3"
    label "small_task"
    publishDir "${params.outdir}/hints/${name}"

    tag "${name}"

    input:
    set val(name), file("augustus.gff3") from augustusJoinedChunks
        .filter { n, p, g -> p == "augustus_hints" }
        .map { n, p, g -> [n, g] }

    output:
    set val(name),
        file("${name}_augustus_hints_hints.gff3") into augustusHintsHints

    script:
    """
      gff2hints.py \
        --source PR \
        --group-level transcript \
        --priority 4 \
        --cds-trim 6 \
        --exon-trim 6 \
        --utr-trim 9 \
        --intron-trim 0 \
        augustus.gff3 \
    | awk 'BEGIN {OFS="\\t"} {sub(/group=/, "group=${name}_augustus_", \$9); print}' \
    > "${name}_augustus_hints_hints.gff3"
    """
}


/*
 */
process filterPredStrand {

    label "posix"
    label "small_task"

    tag "${name}"

    input:
    set val(name), file("*hints") from augustusExtrinsicHints4Preds
            .mix(
                genemarkHints,
                pasaHints.flatMap { n, p, t -> [[n, p], [n, t]] },
                codingQuarryHints,
                codingQuarryPMHints,
                gemomaHints,
                augustusHintsHints,
            )
            .groupTuple(by: 0)

    output:
    set val(name),
        file("${name}_pos_hints_for_augustus_preds.gff"),
        file("${name}_neg_hints_for_augustus_preds.gff") into augustusPredHints4PredTmp

    script:
    """
    mkdir tmp

      cat *hints \
    | gawk '
        \$7 == "-" && (\$9 ~ /group/ || \$9 ~ /grp/) {
          b=gensub(/.*gr(ou)?p=([^;]+).*/, "\\\\2", "g", \$9);
          print b;
        }
      ' \
    | sort -u -T tmp \
    > neg_ids.txt

      cat *hints \
    | gawk '
        \$7 == "+" && (\$9 ~ /group/ || \$9 ~ /grp/) {
          b=gensub(/.*gr(ou)?p=([^;]+).*/, "\\\\2", "g", \$9);
          print b;
        }
      ' \
    | sort -u -T tmp \
    > pos_ids.txt

    rm -rf -- tmp

    cat *hints | grep -f neg_ids.txt -F > neg_groups.gff
    cat *hints | grep -f pos_ids.txt -F > pos_groups.gff

      cat *hints \
    | gawk '(\$7 == "-" || \$7 == ".") && !(\$9 ~ /group/ || \$9 ~ /grp/)' \
    > neg_single.gff

      cat *hints \
    | gawk '(\$7 == "+" || \$7 == ".") && !(\$9 ~ /group/ || \$9 ~ /grp/)' \
    > pos_single.gff

    cat neg_single.gff neg_groups.gff > "${name}_neg_hints_for_augustus_preds.gff"
    cat pos_single.gff pos_groups.gff > "${name}_pos_hints_for_augustus_preds.gff"
    """
}

augustusPredHints4PredTmp
    .flatMap { n, f, r -> [[n, "forward", f], [n, "reverse", r]]}
    .set { augustusPredHints4Pred }

/*
 */
process runAugustusPreds {

    label "augustus"
    label "small_task"

    tag "${name} - ${strand}"

    input:
    set val(name),
        file(fasta),
        val(strand),
        file("hints.gff") from genomes4RunAugustusPreds
            .combine(augustusPredHints4Pred, by: 0)

    file "augustus_config" from augustusConfig
    file "extrinsic.cfg" from augustusPredWeights

    output:
    set val(name),
        val(strand),
        file("out.gff") into augustusPredsResults

    script:
    if ( strand == "forward" ) {
        strand_param = "--strand=forward"
    } else if ( strand == "reverse" ) {
        strand_param = "--strand=backward"
    } else {
        log.error "This shouldn't happen"
        exit 1
    }

    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    perl -n -e'/>(\\S+)/ && print \$1."\\n"' < "${fasta}" > seqids.txt

    getLinesMatching.pl seqids.txt 1 < hints.gff > hints_filtered.gff

    augustus \
      --species="${params.augustus_species}" \
      --extrinsicCfgFile=extrinsic.cfg \
      --hintsfile=hints_filtered.gff \
      ${strand_param} \
      --UTR=on \
      --allow_hinted_splicesites=atac \
      --softmasking=on \
      --alternatives-from-evidence=true \
      --min_intron_len=20 \
      --start=on \
      --stop=on \
      --introns=on \
      --cds=on \
      --gff3=on \
      --outfile="out.gff" \
      --errfile=augustus.err \
      "${fasta}"
    """
}


/*
 */
process joinAugustusPredsChunks {

    label "genometools"
    label "small_task"
    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    input:
    set val(name), file("*chunks.gff") from augustusPredsResults
        .map { n, s, f -> [n, f] }
        .groupTuple(by: 0)

    output:
    set val(name),
        file("${name}_preds.gff3") into augustusJoinedPreds

    script:
    """
    for f in *chunks.gff
    do
        awk 'BEGIN {OFS="\t"} \$3 == "transcript" {\$3="mRNA"} {print}' \${f} \
      | grep -v "^#" \
      | gt gff3 -tidy -sort -o \${f}_tidied.gff3 -
    done

    gt merge -tidy -o "${name}_preds.gff3" *_tidied.gff3
    """
}

// 6 If no genome alignment, run sibelliaz

// 7 Combine estimates using augustus.

// 8 Screen proteins using database of TEs

// 9 stats


/*
 * Evaluate genome completeness with BUSCO on the genomes.
 * Later we evaluate each gene prediction set too.
 * Could compare this number with that one.
 */
process runBusco {
    label "busco"
    label "medium_task"

    tag "${name}"

    publishDir "${params.outdir}/qc/${name}"

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


pasaPredictions4Busco.map { n, g, c, p -> [n, "transdecoder", p] }
    .mix(
        genemarkPredictions4Busco,
        codingQuarryPredictions4Busco
            .map { n, g, c, p, d, f, o -> [n, "codingquarry", p] },
        codingQuarryPMPredictions4Busco
            .map { n, g, c, p -> [n, "codingquarrypm", p] },
    )
    .set { proteins4Busco }


/*
 */
process runBuscoProteins {
    label "busco"
    label "medium_task"

    tag "${name}"

    publishDir "${params.outdir}/qc/${name}"

    when:
    params.busco_lineage

    input:
    set val(name), val(analysis), file(fasta) from proteins4Busco
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
