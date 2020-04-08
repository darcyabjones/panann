#!/usr/bin/env nextflow


/*
 * Get augustus hints from spaln results.
 */
process extract_augustus_hints {

    label "gffpal"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    val analysis
    val source
    val hint_source // E.g. E, P, T, M
    val priority
    val exon_trim
    val cds_trim
    val utr_trim
    val gene_trim
    val is_final
    tuple val(name), path("in.gff3")

    output:
    tuple val(name),
        val(analysis),
        path("${name}_${analysis}_hints.gff3")

    script:

    def use_part = is_final ? "": "part"
    """
    gffpal hints \
        --source "${hint_source}" \
        --group-level mRNA \
        --priority "${priority}" \
        --exon "exon${use_part}" \
        --cds "CDS${use_part}" \
        --utr "UTR${use_part}" \
        --utr3 "UTR${use_part}" \
        --utr5 "UTR${use_part}" \
        --intron "intron" \
        --exon-trim "${exon_trim}" \
        --cds-trim "${cds_trim}" \
        --utr-trim "${utr_trim}" \
        --gene-trim "${gene_trim}" \
        --intron-trim 0 \
        in.gff3 \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_${analysis}_", \$9);
          \$2 = "${source}";
          print
        }
      ' \
    > "${name}_${analysis}_hints.gff3"
    """
}


/*
 * Get augustus hints from spaln results.
 */
process extract_augustus_split_hints {

    label "gffpal"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    val analysis
    val exon_source
    val cds_source
    val exon_hint_source  // e.g. E, I
    val cds_hint_source  // e.g. P, T
    val exon_priority
    val cds_priority
    val exon_trim
    val cds_trim
    val utr_trim
    val gene_trim
    val is_final
    tuple val(name), path("in.gff3")

    output:
    tuple val(name),
        val(analysis),
        path("${name}_${analysis}_hints.gff3")

    script:

    def use_part = is_final ? "": "part"
    """
    awk -F '\t' '\$3 == "exon" || \$3 == "intron" || \$3 == "mRNA"' in.gff3 \
    | gffpal hints \
        --source "${exon_hint_source}" \
        --group-level mRNA \
        --priority "${exon_priority}" \
        --exon "exon${use_part}" \
        --intron "intron" \
        --exon-trim "${exon_trim}" \
        --gene-trim "${gene_trim}" \
        --intron-trim 0 \
        - \
    | awk '
        BEGIN {OFS="\\t"}
        {
          sub(/group=/, "group=${name}_${exon_source}_", \$9);
          \$2 = "${exon_source}";
          print
        }
      ' \
    > exon_hints.gff3

    awk -F '\t' '\$3 != "exon" && \$3 != "intron"' in.gff3 \
    | gffpal hints \
        --source "${cds_hint_source}" \
        --group-level mRNA \
        --priority "${cds_priority}" \
        --cds "CDS${use_part}" \
        --utr "UTR${use_part}" \
        --utr3 "UTR${use_part}" \
        --utr5 "UTR${use_part}" \
        --cds-trim "${cds_trim}" \
        --utr-trim "${utr_trim}" \
        - \
    | awk -F '\t' '
        BEGIN {OFS="\\t"}
        \$3 != "genicpart" {
          sub(/group=/, "group=${name}_${cds_source}_", \$9);
          \$2 = "${cds_source}";
          print
        }
      ' \
    > "cds_hints.gff3"

    cat exon_hints.gff3 cds_hints.gff3 > "${name}_${analysis}_hints.gff3"
    """
}


/*
 * Convert exonerate alignments to augustus hints
 * This is a bit of a special case because it's gff2
 */
process extract_exonerate_hints {

    label "braker"
    label "small_task"
    time '2h'

    tag "${genome_name} - ${protein_name}"

    input:
    val min_intron_hard
    val max_intron_hard
    tuple val(genome_name),
          val(protein_name),
          path("exonerate.gff")

    output:
    tuple val(genome_name),
          val(protein_name),
          path("${genome_name}_${protein_name}_exonerate_hints.gff3")

    script:
    """
    align2hints.pl \
      --in=exonerate.gff \
      --out=hints.gff3 \
      --prg=exonerate \
      --CDSpart_cutoff=15 \
      --minintronlen="${min_intron_hard}" \
      --maxintronlen="${max_intron_hard}" \
      --priority=2 \
      --source=T

    awk -F '\\t' '
      BEGIN {OFS="\\t"}
      \$3 == "CDSpart" {
        sub(/grp=/, "grp=${genome_name}_${protein_name}_exonerate_", \$9)
        \$2 = "exonerate";
        print
      }
      ' \
      hints.gff3 \
    > "${genome_name}_${protein_name}_exonerate_hints.gff3"
    """
}


/*
 * Convert exonerate alignments to augustus hints
 * This is a bit of a special case because it's gff2
 */
process extract_exonerate_evm_hints {

    label "braker"
    label "small_task"
    time '2h'

    tag "${genome_name} - ${protein_name}"

    input:
    val min_intron_hard
    val max_intron_hard
    set val(genome_name),
        val(protein_name),
        file("exonerate.gff")

    output:
    set val(genome_name),
        val(protein_name),
        file("${genome_name}_${protein_name}_exonerate_evm_hints.gff3")

    script:
    """
    align2hints.pl \
      --in=exonerate.gff \
      --out=evm.gff3 \
      --prg=exonerate \
      --CDSpart_cutoff=0 \
      --minintronlen="${min_intron_hard}" \
      --maxintronlen="${max_intron_hard}" \
      --priority=2 \
      --source=T

    awk -F '\\t' '
      BEGIN {OFS="\\t"}
      \$3 == "CDSpart" {
        id=gensub(/.*grp=([^;]+).*/, "\\\\1", "g", \$9);
        \$9="ID=${genome_name}_${protein_name}_exonerate_"id;
        \$2 = "exonerate";
        \$3 = "nucleotide_to_protein_match";
        print
      }
      ' \
      evm.gff3 \
    > "${genome_name}_${protein_name}_exonerate_evm_hints.gff3"
    """
}


/*
 * Extract hints to be used for augustus and genemark.
 *
 * We only use the intron hints because we have spaln/pasa alignments,
 * and the coverage info causes lots of close genes to be merged or extended
 * even if you use the UTR model.
 */
process extract_augustus_rnaseq_hints {

    label "braker"
    label "medium_task"
    time '4h'

    tag "${name} - ${read_group}"

    input:
    val min_intron_hard
    val max_intron_hard
    val valid_splicesites
    tuple val(name),
        val(read_group),
        path(fasta),
        path(cram)

    output:
    tuple val(name),
        val(read_group),
        path("${name}_${read_group}_intron_hints.gff3")

    script:
    max_gap_len = min_intron_hard - 1

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
      --maxgaplen="${max_gap_len}" \
      --minintronlen="${min_intron_hard}" \
      --maxintronlen="${max_intron_hard}" \
      --maxcoverage=1000 \
      --priority=4 \
      --ssOn \
      --source="I" \
      --in="tmp.bam" \
      --out="tmp.gff3"

    filterIntronsFindStrand.pl \
        "${fasta}" \
        tmp.gff3 \
        --allowed="${valid_splicesites}" \
        --score \
    > "${name}_${read_group}_intron_hints.gff3"

    rm -f tmp.gff3
    """
}


/*
 * Run the GeMoMa pipeline
 */
process extract_gemoma_rnaseq_hints {

    label "gemoma"
    label "medium_task"
    time '2h'

    tag "${name} - ${read_group}"

    input:
    tuple val(name),
        val(read_group),
        path(fasta),
        path(cram),
        val(strand)

    output:
    tuple val(name),
        val(read_group),
        path("${name}_${read_group}_gemoma_introns.gff"),
        path("${name}_${read_group}_gemoma_forward.bedgraph"),
        path("${name}_${read_group}_gemoma_reverse.bedgraph")

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

    mv introns.gff "${name}_${read_group}_gemoma_introns.gff"
    mv coverage_forward.bedgraph "${name}_${read_group}_gemoma_forward.bedgraph"
    mv coverage_reverse.bedgraph "${name}_${read_group}_gemoma_reverse.bedgraph"

    rm -f *.bam protocol_ERE.txt
    rm -rf -- GeMoMa_temp
    """
}


/*
 * Gemoma joins the rnaseq hints into one.
 * I think it's basically just a cat | sort.
 */
process combine_gemoma_rnaseq_hints {

    label "gemoma"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    tuple val(name),
        tuple("*i.gff"),
        tuple("*f.bedgraph"),
        tuple("*r.bedgraph")

    output:
    tuple val(name),
        path("${name}_gemoma_introns.gff"),
        path("${name}_gemoma_forward.bedgraph"),
        path("${name}_gemoma_reverse.bedgraph")

    script:
    """
    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineIntronFiles \
      ${name}_gemoma_introns.gff \
      *i.gff

    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineCoverageFiles \
      "${name}_gemoma_forward.bedgraph" \
      *f.bedgraph

    java -cp \${GEMOMA_JAR} \
      projects.gemoma.CombineCoverageFiles \
      "${name}_gemoma_reverse.bedgraph" \
      *r.bedgraph
    """
}


process extract_gmap_evm_hints {

    label "posix"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    tuple val(name),
          path("in.gff3")

    output:
    tuple val(name),
          path("${name}_gmap_evm_hints.gff3")

    script:
    """
    awk -F'\t' '
      BEGIN { OFS="\t" }
      \$3 == "cDNA_match" {
        \$2="gmap";
        print
      }
    ' in.gff3 \
    >> "${name}_gmap_evm_hints.gff3"
    """
}


process extract_spaln_transcript_evm_hints {

    label "posix"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    tuple val(name),
          path("in.gff3")

    output:
    tuple val(name),
          path("${name}_spaln_transcript_evm_hints.gff3")

    script:
    """
    awk -F'\t' '
      BEGIN { OFS="\t" }
      \$3 == "exon" {
        parent=gensub(/.*Parent=([^;]+).*/, "\\\\1", "g", \$9);
        target=gensub(/.*Target=([^;]+).*/, "\\\\1", "g", \$9);
        \$9="ID=${name}_spaln_transcript" parent ";Target=" target;
        \$2="spaln_transcript";
        \$3="cDNA_match";
        print
      }
    ' in.gff3 \
    >> "${name}_spaln_transcript_evm_hints.gff3"
    """
}


process extract_spaln_protein_evm_hints {

    label "posix"
    label "small_task"
    time "1h"

    tag "${name}"

    input:
    tuple val(name),
          path("in.gff3")

    output:
    tuple val(name),
          path("${name}_spaln_protein_evm_hints.gff3")

    script:
    """
    awk -F'\t' '
      BEGIN { OFS="\t" }
      \$3 == "CDS" || \$3 == "cds" {
        parent=gensub(/.*Parent=([^;]+).*/, "\\\\1", "g", \$9);
        target=gensub(/.*Target=([^;]+).*/, "\\\\1", "g", \$9);
        \$9="ID=${name}_spaln_protein" parent ";Target=" target;
        \$2="spaln_protein";
        \$3="nucleotide_to_protein_match";
        print
      }
    ' in.gff3 \
    >> "${name}_spaln_protein_evm_hints.gff3"
    """
}
