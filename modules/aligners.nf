#!/usr/bin/env nextflow

/*
 * Index the genome for alignment with Spaln.
 */
process get_spaln_index {

    label "spaln"
    label "small_task"
    time '3h'

    tag "${name}"

    input:
    tuple val(name), path("in.fasta")

    output:
    tuple val(name),
        path("${name}_spaln_index.bkn"),
        path("${name}_spaln_index.ent"),
        path("${name}_spaln_index.idx"),
        path("${name}_spaln_index.bkp"),
        path("${name}_spaln_index.grp"),
        path("${name}_spaln_index.seq")

    script:
    """
    makeidx.pl -inp in.fasta

    mv "in.bkn" "${name}_spaln_index.bkn"
    mv "in.ent" "${name}_spaln_index.ent"
    mv "in.idx" "${name}_spaln_index.idx"
    mv "in.bkp" "${name}_spaln_index.bkp"
    mv "in.grp" "${name}_spaln_index.grp"
    mv "in.seq" "${name}_spaln_index.seq"
    """
}


/*
 * Align transcripts using Spaln
 * This performs pretty well, especially if you have a
 * good set of consensus splice sites in their package.
 */
process spaln_align_transcripts {

    label "spaln"
    label "medium_task"
    time '5h'

    tag "${name}"

    input:
    val species
    val max_gene_hard
    val min_intron_soft
    tuple val(name),
        path("db.bkn"),
        path("db.ent"),
        path("db.idx"),
        path("db.bkp"),
        path("db.grp"),
        path("db.seq"),
        path("transcripts.fasta")

    output:
    tuple val(name), path("${name}_spaln_transcripts.gff3")

    script:
    def species_params = species ? "-T${species} -yS " : ""

    """
    spaln \
      -LS \
      -O0 \
      -Q7 \
      -S3 \
      -yX \
      -ya1 \
      ${species_params} \
      -XG ${max_gene_hard} \
      -yL${min_intron_soft} \
      -t ${task.cpus} \
      -d db \
      "transcripts.fasta" \
    > "${name}_spaln_transcripts.gff3"
    """
}


/*
 * Align all proteins to genome with Spaln.
 * Spaln is good for closely related proteins, but not great for distant ones.
 */
process spaln_align_proteins {

    label "spaln"
    label "medium_task"
    time '6h'

    tag "${name}"

    input:
    val trans_table
    val min_intron_soft
    val max_gene_hard
    tuple val(name),
        path("db.bkn"),
        path("db.ent"),
        path("db.idx"),
        path("db.bkp"),
        path("db.grp"),
        path("db.seq"),
        path("proteins.fasta")

    output:
    tuple val(name), path("${name}_spaln_proteins.gff3")

    script:
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
      -yL${min_intron_soft} \
      -XG${max_gene_hard} \
      -t ${task.cpus} \
      -d "${name}" \
      "proteins.fasta" \
    > "${name}_spaln_proteins.gff3"
    """
}


/*
 * Spaln CDS output doesn't include the stop codon.
 * Fix this.
 * the type is also "cds" instead of "CDS"
 */
process fix_spaln_proteins_stop {

    label "gffpal"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    tuple val(name), path("spaln.gff3")

    output:
    tuple val(name), path("${name}_spaln_proteins_fixed_stop.gff3")

    script:
    """
      awk -F '\\t' 'BEGIN {OFS="\\t"} \$3 == "cds" {\$3="CDS"} {print}' "spaln.gff3" \
    | gffpal expandcds -o "${name}_spaln_proteins_fixed_stop.gff3" --cds-type "CDS" -
    """
}


/*
 * Index the genomes for GMAP
 */
process get_gmap_index {

    label "gmap"
    label "medium_task"
    time '3h'

    tag "${name}"

    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path("${name}_gmap_index")

    script:
    """
    gmap_build \
      -k 13 \
      -D "${name}_gmap_index" \
      -d "${name}" \
      ${genome}
    """
}


/*
 * Align transcripts using map
 * We run this separately to PASA to control the number of threads PASA uses.
 * If you use BLAT + gmap in pasa it can run 2*cpus threads, which would screw
 * with our provisioning. Most of the PASA time is spent running transdecoder.
 */
process gmap_align_transcripts {

    label "gmap"
    label "medium_task"
    time '6h'

    tag "${name}"

    input:
    val min_intron_hard
    val max_intron_hard
    tuple val(name),
        path("db"),
        path("transcripts.fasta")

    output:
    tuple val(name), path("${name}_gmap_transcripts.gff3")

    script:
    def trim_end_exons = 12
    def microexon_spliceprob = 0.95
    def canonical_mode = 1
    def cross_species = "" // "--cross-species "

    """
    gmap \
      --npaths=0 \
      --chimera-margin=50 \
      --min-intronlength="${min_intron_hard}" \
      --max-intronlength-middle="${max_intron_hard}" \
      --max-intronlength-ends="${max_intron_hard}" \
      --trim-end-exons="${trim_end_exons}" \
      --microexon-spliceprob="${microexon_spliceprob}" \
      --canonical-mode=1 \
      ${cross_species} \
      --format=gff3_match_cdna \
      --nthreads "${task.cpus}" \
      -D db \
      -d "${name}" \
      transcripts.fasta \
    > ${name}_gmap_transcripts.gff3
    """
}


/*
 * Preindex MMSeqs remote proteins
 */
process get_mmseqs_protein_db {

    label "mmseqs"
    label "small_task"
    time '3h'

    tag "${name}"

    input:
    tuple val(name), path("seqs.fasta")

    output:
    tuple val(name), path("${name}_mmseqs_protein_db")

    script:
    """
    mkdir -p tmp "${name}_mmseqs_protein_db"

    mmseqs createdb seqs.fasta "${name}_mmseqs_protein_db/db"
    # mmseqs createindex "${name}_mmseqs_protein_db/db" tmp --threads "${task.cpus}"

    rm -rf -- tmp
    """
}


/*
 * Preindex MMSeqs remote proteins
 */
process get_mmseqs_genome_db {

    label "mmseqs"
    label "small_task"
    time '3h'

    tag "${name}"

    input:
    tuple val(name), path("seqs.fasta")

    output:
    tuple val(name), path("${name}_mmseqs_genome_db")

    script:
    """
    mkdir -p "${name}_mmseqs_genome_db"
    mmseqs createdb seqs.fasta "${name}_mmseqs_genome_db/db" --dont-split-seqs-by-len
    """
}


/*
 * Quickly find genomic regions with matches to remote proteins.
 * This is an approximate method which we refine later with exonerate.
 */
process mmseqs_search_genome_against_proteins {

    label "mmseqs"
    label "big_task"
    time '6h'

    tag "${genome_name} - ${protein_name}"

    input:
    val trans_table
    tuple val(genome_name),
        path("genome"),
        val(protein_name),
        path("proteins")

    output:
    tuple val(genome_name),
        val(protein_name),
        path("${genome_name}_${protein_name}_genome_v_protein_mmseqs_matches.tsv")

    script:
    """
    mkdir result tmp
    cp -rL proteins proteins_tmp

    # Searching with genome as query is ~3X faster
    mmseqs search \
      genome/db \
      proteins_tmp/db \
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
      --translation-table "${trans_table}" \
      --use-all-table-starts

    # Extract match results.
    mmseqs convertalis \
      genome/db \
      proteins_tmp/db \
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
    > "${genome_name}_${protein_name}_genome_v_protein_mmseqs_matches.tsv"

    sed -i '1i query\ttarget\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tpident\tmismatch\tgapopen\tevalue\tbitscore' "${genome_name}_${protein_name}_genome_v_protein_mmseqs_matches.tsv"

    rm -rf -- tmp genome result results_unsorted.tsv proteins_tmp
    """
}


process cluster_genome_vs_protein_matches {

    label "bedtools"
    label "small_task"
    time "2h"

    tag "${name}"

    input:
    val pad_size // We extend the region to align against by this many basepairs. Suggest max_gene_hard.
    val merge_distance // 1000 is good
    tuple val(genome_name),
        path(faidx),
        path(protein_name),
        path("matches.tsv")

    output:
    tuple val(genome_name),
        val(protein_name),
        path("${genome_name}_${protein_name}_genome_v_proteins_clustered_matches.bed")

    script:
    """
    mkdir -p tmp

      tail -n+2 matches.tsv \
    | awk '
        BEGIN { OFS="\t" }
        \$3 > \$4 { print \$1, \$4, \$3, \$2 }
        \$3 < \$4 { print \$1, \$3, \$4, \$2 }
      ' \
    | sort \
        -k1,1 -k2,2n -k3,3n \
        --temporary-directory=tmp \
    | sed 's/,/%2C/g' \
    | bedtools merge -d "${merge_distance}" -c 4 -o distinct -i - \
    | bedtools slop -g "${faidx}" -b "${pad_size}" -i - \
    > "${genome_name}_${protein_name}_genome_v_proteins_clustered_matches.bed"

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
process exonerate_regions {

    label "exonerate"
    label "big_task"
    time '6h'

    tag "${name}"

    input:
    val trans_table
    val min_intron_hard
    val max_intron_hard
    tuple val(genome_name),
        fasta(fasta),
        val(protein_name),
        fasta("clustered.bed"),
        fasta("proteins.tsv")

    output:
    set val(genome_name),
        val(protein_name),
        file("${genome_name}_${protein_name}_genome_v_proteins_exonerate.gff")

    script:
    """
    mkdir -p tmp
    # Sometimes this gets touched, so the checkpointing goes a bit skiwiff.
    cp -L proteins.tsv proteins_tmp.tsv

    exonerate_parallel.sh \
      -g "${fasta}" \
      -q "proteins_tmp.tsv" \
      -b "clustered.bed" \
      -n "${task.cpus}" \
      -t "tmp" \
      -m "${min_intron_hard}" \
      -x "${max_intron_hard}" \
      -r "${trans_table}" \
      -o "${genome_name}_${protein_name}_genome_v_proteins_exonerate.gff"

    rm -rf -- tmp proteins_tmp.tsv
    """
}



/*
 * Index genome for STAR
 */
process get_star_index {

    label "star"
    label "medium_task"
    time '4h'

    tag "${name}"

    input:
    tuple val(name),
        path(fasta),
        path(gff)

    output:
    tuple val(name), file("${name}_star_index")

    script:
    // If no gff was in known_sites, use it, otherwise dont
    def sjdb = gff.name != 'WAS_NULL' ? "--sjdbGTFfile ${gff} --sjdbOverhang 149 " : ''
    // Possibly need option to set this to CDS, incase user input doesn't have exons?
    def exon_feature = "exon"

    """
    mkdir -p "${name}_star_index"
    STAR \
      --runThreadN ${task.cpus} \
      --runMode genomeGenerate \
      --genomeDir "${name}_star_index" \
      --genomeFastaFiles "${fasta}" \
      --genomeSAindexNbases 11 \
      --sjdbGTFtagExonParentTranscript Parent \
      --sjdbGTFfeatureExon "${exon_feature}" \
      ${sjdb}
    """
}


/*
 * Perform first pass for STAR.
 *
 * This finds novel splice sites in the rnaseq.
 * We filter out SS with poor support later.
 */
process star_find_splicesites {

    label "star"
    label "medium_task"
    time '1d'

    tag "${name} - ${read_group}"

    input:
    val min_intron_len
    val max_intron_len
    val extra_params
    tuple val(name),
        path("index"),
        val(read_group),
        path(r1s),
        path(r2s)

    output:
    tuple val(name),
        val(read_group),
        path("${name}_${read_group}.SJ.out.tab")

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
      ${extra_params} \
      --alignIntronMin ${min_intron_len} \
      --alignIntronMax ${max_intron_len} \
      --alignSJoverhangMin 10 \
      --alignSJDBoverhangMin 3 \
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
process star_align_reads {

    label "star"
    label "medium_task"
    time '1d'

    tag "${name} - ${read_group}"

    input:
    val min_intron_len
    val max_intron_len
    val extra_params
    tuple val(name),
        val(read_group),
        path("genome.fasta"),
        path("index"),
        path(r1s),
        path(r2s),
        path("*SJ.out.tab")

    output:
    tuple val(name),
        val(read_group),
        path("${name}_${read_group}.cram")

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
      ${extra_params} \
      --alignIntronMin ${min_intron_len} \
      --alignIntronMax ${max_intron_len} \
      --alignSJoverhangMin 10 \
      --alignSJDBoverhangMin 3 \
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
        -T "genome.fasta" \
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
