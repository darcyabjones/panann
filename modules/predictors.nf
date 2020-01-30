#!/usr/bin/env nextflow

/*
 * Do denovo gene prediction with intron hint training.
 */
process genemark {

    label "genemarkes"
    label "medium_task"
    time '12h'

    tag "${name}"

    input:
    val not_fungus
    val training
    tuple val(name),
        path(genome),
        path(faidx),
        path("*introns.gff3")

    output:
    tuple val(name), path("${name}_genemark.gtf")

    script:
    def use_fungus = not_fungus ? '' : '--fungus '
    def is_training = training ? '--min_contig 300 ': ''

    """
    sort -k1,1V -k4,4n -k5,5rn -k3,3r *introns.gff3 > hints.gff3

    gmes_petap.pl \
      --cores "${task.cpus}" \
      --soft_mask 100 \
      --ET "hints.gff3" \
      ${use_fungus} \
      ${is_training} \
      --sequence "${genome}"

    mv genemark.gtf "${name}_genemark.gtf"
    """
}

/*
 * Predict genes using pasa and transdecoder.
 */
process pasa {

    label "pasa"
    label "small_task"
    time '1d'

    tag "${name}"

    input:
    val not_fungus
    val max_intron_hard
    tuple val(name),
        path(genome_fasta),
        path(known_sites),
        path(stringtie_gtf),
        path(gmap_aligned),
        path(transcripts_fasta),
        path(transcripts_fasta_cln),
        path(transcripts_fasta_clean)

    output:
    tuple val(name), path("${name}_pasa.gff3")

    script:
    def use_stringent = params.notfungus ? '' : "--stringent_alignment_overlap 30.0 "

    // Don't use stringtie if it is fungus
    // Stringtie and cufflinks tend to merge overlapping features,
    // which doesn't work well for organisms with high gene density.
    def use_stringtie = (stringtie_gtf.name == "SGT_WAS_NULL" || !not_fungus) ? '' : "--trans_gtf ${stringtie_gtf} "
    def use_known = known_sites.name == "WAS_NULL" ? '' : "-L --annots ${known_sites} "

    // Transdecoder doesn't support standard integer based table access.
    def gen_code = "Universal"

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
      --MAX_INTRON_LENGTH "${max_intron_hard}" \
      --ALIGNERS blat \
      --CPU "${task.cpus}" \
      --transcribed_is_aligned_orient \
      --TRANSDECODER \
      ${use_stringent} \
      ${use_stringtie} \
      ${use_known}

    pasa_asmbls_to_training_set.dbi \
      -G "${gen_code}" \
      --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta \
      --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3

    mv \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 \${PWD}/${name}_pasa.gff3
    """
}


/*
 * Predict genes using codingquarry
 * Note that we still need the proteins extracted from here for codingquarrypm.
 */
process codingquarry {

    label "codingquarry"
    label "big_task"
    time '1d'

    errorStrategy "retry"
    maxRetries 10

    tag "${name}"

    input:
    tuple val(name),
        path("transcripts.gtf"),
        path("genome.fasta")

    output:
    tuple val(name), path("${name}_codingquarry_fixed.gff3")
    tuple val(name), path("${name}_codingquarry.gff3") // We need to hold onto this because cqpm doesn't like the fixed versions.
    tuple val(name), path("${name}_codingquarry.faa")
    path "${name}_codingquarry.fna"
    path "${name}_codingquarry_dubiousset.gff3"
    path "${name}_codingquarry_fusions.txt"
    path "${name}_codingquarry_overlapreport.txt"

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

    # Tidy some of the CQ weirdness.
    # Sometimes you can get a -ve cds phase for unknown reasons.
    # CQ puts the CDS directly on the gene, which breaks some tools.
    awk -F '\t' 'BEGIN {OFS="\\t"} \$8 = "-1" {\$8="0"} {print}' out/PredictedPass.gff3 \
    | awk -F '\t' \
      '
        BEGIN {OFS="\\t"}
        \$3 == "gene" {
            print;
            \$9=gensub(/.*ID=([^;]+).*/, "ID=mRNA:\\\\1;Parent=\\\\1;", "1", \$9);
            \$3="mRNA";
            print
        }
        \$3 == "CDS" {
            \$9=gensub(/Parent=/, "Parent=mRNA:", "1", \$9);
            print
        }
      ' \
    > "${name}_codingquarry_fixed.gff3"

    mv out/PredictedPass.gff3 "${name}_codingquarry.gff3"

    mv out/DubiousSet.gff3 "${name}_codingquarry_dubiousset.gff3"
    mv out/Predicted_CDS.fa "${name}_codingquarry.fna"
    mv out/Predicted_Proteins.faa "${name}_codingquarry.faa"
    mv out/fusions.txt "${name}_codingquarry_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarry_overlapreport.txt"

    rm -rf -- out
    """
}


process signalp {

    label "signalp"
    label "medium_task"
    time '12h'

    tag "${name}"

    input:
    tuple val(name),
        path("proteins.faa")

    output:
    tuple val(name),
        path("${name}_secreted.txt")

    script:
    """
    mkdir tmp
    signalp \
      -fasta "proteins.faa" \
      -prefix "cq" \
      -org euk \
      -tmp tmp

    rm -rf -- tmp

    gawk '
      BEGIN {
        OFS=" "
      }
      \$2 ~ /^SP/ {
        match(\$0, /CS pos: ([0-9]*)-/, x)
        print \$1, x[1]
      }
    ' < "cq_summary.signalp5" \
    > "${name}_secreted.txt"
    """
}


process deepsig {

    label "deepsig"
    label "medium_task"
    time '12h'

    tag "${name}"

    input:
    tuple val(name),
        path("proteins.faa")

    output:
    tuple val(name),
        path("${name}_secreted.txt")

    script:
    """
    deepsig.py \
      -f "proteins.faa" \
      -o secreted.txt \
      -k euk

    gawk '
      BEGIN {
        OFS=" "
      }
      \$2 == "SignalPeptide" {
        print \$1, \$4
      }
    ' < secreted.txt > "${name}_secreted.txt"
    """
}


/*
 * CQPM looks for genes that might not be predicted main set because of
 * genome compartmentalisation or differences with signal peptides.
 * NOTE: This fails if there are fewer than 500 secreted genes to train from.
 *
 * CQPM also uses a lot of memory at a specific point, which can cause
 * random segfaults. We retry a few times, but if it keeps failing you
 * probably need to increase the available RAM.
 *
 * Try setting `vm.overcommit_memory = 1` if you get segfaults before
 * reaching max memory.
 *
 * We also set the output to be optional because it requires > 500 proteins to train from.
 */
process codingquarrypm {

    label "codingquarry"
    label "big_task"
    time '1d'

    errorStrategy "retry"
    maxRetries 10

    tag "${name}"

    input:
    tuple val(name),
        path("transcripts.gtf"),
        path("genome.fasta"),
        path("codingquarry.gff3"),
        path("secretome.txt")

    output:
    tuple val(name),
        path("${name}_codingquarrypm.gff3") optional true
    path "${name}_codingquarrypm_fusions.txt" optional true
    path "${name}_codingquarrypm_overlapreport.txt" optional true

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
      -2 "codingquarry.gff3" \
      -p "${task.cpus}" \
      -g secretome.txt \
      -h

    # Tidy some of the CQ weirdness.
    # Sometimes you can get a -ve cds phase for unknown reasons.
    # CQ puts the CDS directly on the gene, which breaks some tools.
    awk -F '\t' 'BEGIN {OFS="\\t"} \$8 = "-1" {\$8="0"} {print}' out/PGN_predictedPass.gff3 \
    | awk -F '\t' \
      '
        BEGIN {OFS="\\t"}
        \$3 == "gene" {
            print;
            \$9=gensub(/.*ID=([^;]+).*/, "ID=mRNA:\\\\1;Parent=\\\\1;", "1", \$9);
            \$3="mRNA";
            print
        }
        \$3 == "CDS" {
            \$9=gensub(/Parent=/, "Parent=mRNA:", "1", \$9);
            print
        }
      ' \
    > "${name}_codingquarrypm.gff3"

    mv out/fusions.txt "${name}_codingquarrypm_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarrypm_overlapreport.txt"
    rm -rf -- out
    """
}


/*
 * Gemoma needs to extract the splice site info as well as the proteins.
 * Can probably replace this with the comparative version below?
 */
process extract_gemoma_cds_parts {

    label "gemoma"
    label "small_task"
    time '2h'

    tag "${name}"

    input:
    tuple val(name),
        path(fasta),
        path(gff)

    output:
    tuple val(name),
        path("${name}_gemoma_cds_parts.fasta"),
        path("${name}_gemoma_assignment.tabular"),
        path("${name}_gemoma_proteins.fasta")

    script:
    """
    java -jar \${GEMOMA_JAR} CLI Extractor \
      a=${gff} \
      g=${fasta} \
      p=true \
      outdir=.

    mv cds-parts.fasta "${name}_gemoma_cds_parts.fasta"
    mv assignment.tabular "${name}_gemoma_assignment.tabular"
    mv proteins.fasta "${name}_gemoma_proteins.fasta"

    rm -rf -- GeMoMa_temp
    rm protocol_Extractor.txt
    """
}


/*
 * Extract CDS parts for each individual prediction set.
 */
process extract_gemoma_comparative_cds_parts {

    label "gemoma"
    label "small_task"
    time '3h'

    tag "${name} - ${analysis}"

    input:
    set val(name),
        val(analysis),
        file(fasta),
        file(gff)

    output:
    set val(name),
        val(analysis),
        file("${name}_${analysis}_cdsparts.fasta"),
        file("${name}_${analysis}_assignment.tsv"),
        file("${name}_${analysis}_proteins.fasta")

    script:
    """
    awk -F'\\t' -v name="${name}" -v analysis="${analysis}" '
      BEGIN { OFS="\\t" }
      !/^#/ {
        \$1=name"."\$1;
        \$9=gensub(/Parent=([^;]+)/, "Parent="name"."analysis".\\\\1", "g", \$9);
        \$9=gensub(/ID=([^;]+)/, "ID="name"."analysis".\\\\1", "g", \$9);
        print
      }
    ' < "${gff}" \
    > renamed.gff3

    sed "/^>/s/^>/>${name}./" < "${fasta}" > renamed.fasta

    java -jar \${GEMOMA_JAR} CLI Extractor \
      a=renamed.gff3 \
      g=renamed.fasta \
      p=true \
      outdir=.

    mv cds-parts.fasta "${name}_${analysis}_cdsparts.fasta"
    mv assignment.tabular "${name}_${analysis}_assignment.tsv"
    mv proteins.fasta "${name}_${analysis}_proteins.fasta"

    rm -rf -- GeMoMa_temp
    rm protocol_Extractor.txt
    """
}


/*
 * Concatenate the gemoma parts for all isolates and prediction
 * methods into single files.
 */
/*
 * To reduce the number of alignments/prediction steps, we
 * do a pretty strict clustering to remove redundancy, and
 * use the seed sequence as the representative.
 * The coverage requirements and high identity means that this
 * should mostly just remove nearly identical matches.
 */
process cluster_gemoma_cds_parts {

    label "mmseqs"
    label "medium_task"

    time '1d'

    input:
    tuple path("*c.fasta"), path("*a.tsv"), path("*p.fasta")

    output:
    tuple path("cdsparts.fasta"),
          path("assignment.tsv"),
          path("proteins.fasta")

    script:
    """
    mkdir -p proteins protein_clusters tmp

    cat *c.fasta > "old_cdsparts.fasta"
    cat *p.fasta > "old_proteins.fasta"

    FIRST_FILE=1
    for f in *a.tsv;
    do
      if [ \${FIRST_FILE} == 1 ];
      then
        cat "\${f}" > old_assignment.tsv
        FIRST_FILE=0
      else
        tail -n+2 "\${f}" >> old_assignment.tsv
      fi
    done

    mmseqs createdb old_proteins.fasta proteins/db
    mmseqs cluster \
      proteins/db \
      protein_clusters/db \
      tmp \
      --threads "${task.cpus}" \
      --min-seq-id 0.9 \
      -c 0.98 \
      --cov-mode 0 \
      --cluster-mode 0

    mmseqs createtsv proteins/db proteins/db protein_clusters/db protein_clusters.tsv

    rm -rf -- proteins protein_clusters tmp

    # The script outputs are hardcoded.
    select_comparative_proteins.py \
      --clusters protein_clusters.tsv \
      --assignments old_assignment.tsv \
      --cdsparts old_cdsparts.fasta \
      --proteins old_proteins.fasta
    """
}


/*
 * Align proteins to genomes for Gemoma.
 * Do this separately as the gemoma pipeline
 * currently crashes and we can control parallelism better.
 * The output format order is important.
 * I had to reverse engineer it from the GeMoMa jar binaries!
 */
process mmseqs_search_gemoma_cds_parts {

    label "mmseqs"
    label "medium_task"
    time '6h'

    tag "${target_name} - ${ref_name}"

    input:
    val trans_table
    tuple val(ref_name),
        path("cds-parts.fasta"),
        path("assignment.tabular"),
        path("proteins.fasta"),
        val(target_name),
        path("genome") // Should be mmseqs db

    output:
    set val(target_name),
        val(ref_name),
        file("${target_name}_${ref_name}_mmseqs_search_gemoma_cds_parts_matches.tsv")

    script:
    """
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
      --translation-table "${trans_table}" \
      --use-all-table-starts

    mmseqs convertalis \
      proteins/db \
      genome/db \
      alignment/db \
      "${target_name}_${ref_name}_mmseqs_search_gemoma_cds_parts_matches.tsv" \
      --threads ${task.cpus} \
      --format-mode 0 \
      --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen'

    rm -rf -- proteins alignment tmp
    """
}


/*
 * Predict genes with gemoma for each known-sites set.
 * Merges adjacent MMseqs matches and checks intron-exon boundaries.
 */
process gemoma {

    label "gemoma"
    label "small_task"
    time '6h'

    tag "${target_name} - ${ref_name}"

    input:
    tuple val(target_name),
        val(ref_name),
        path(fasta),
        path("cds-parts.fasta"),
        path("assignment.tabular"),
        path("proteins.fasta"),
        path("matches.tsv"),
        path("introns.gff")

    output:
    tuple val(target_name),
        val(ref_name),
        path("${target_name}_${ref_name}_gemoma.gff3")

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
      r=2

    mv out/predicted_annotation.gff "${target_name}_${ref_name}_gemoma.gff3"

    rm -rf -- GeMoMa_temp out
    """
}


/*
 * Combines gemoma predictions from multiple known-sites sets.
 * Also adds UTRs etc to gff based on RNAseq.
 */
process gemoma_combine {

    label "gemoma"
    label "small_task"
    time '6h'

    tag "${name}"

    input:
    tuple val(name),
        val(ref_names),
        path(pred_gffs),
        path(fasta),
        path("introns.gff"),

    output:
    tuple val(name),
        path("${name}_gemoma_combined.gff3")

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

    # if ${get_utr}
    # then
    #   mkdir -p finalised
    #   java -jar \${GEMOMA_JAR} CLI AnnotationFinalizer \
    #     g=${fasta} \
    #     a=gaf/filtered_predictions.gff \
    #     i=introns.gff \
    #     u=YES \
    #     c=STRANDED \
    #     coverage_forward=coverage_forward.bedgraph \
    #     coverage_reverse=coverage_reverse.bedgraph \
    #     outdir=finalised \
    #     rename=NO

    #   mv finalised/final_annotation.gff gemoma_tmp.gff3
    # else
    #   mv gaf/filtered_predictions.gff gemoma_tmp.gff3
    # fi

    mv gaf/filtered_predictions.gff gemoma_tmp.gff3

    awk 'BEGIN {OFS="\\t"} \$3 == "prediction" {\$3="mRNA"} {print}' \
      gemoma_tmp.gff3 > "${name}_gemoma_combined.gff3"

    rm -rf -- gaf finalised GeMoMa_temp gemoma_tmp.gff3
    """
}


/*
 * Augustus denovo is mosly just run for QC to see what it looks like
 * without hints.
 * We don't actually use the output for the final predictions.
 */
process augustus_denovo {

    label "augustus"
    label "small_task"
    time '1d'

    tag "${name} - ${strand}"

    input:
    val augustus_species
    val augustus_utr
    val not_fungus
    val min_intron_hard
    tuple val(name), val(strand), path(fasta)
    path "augustus_config"

    output:
    tuple val(name),
        val(strand),
        path("${name}_${strand}_augustus_denovo.gff")

    script:
    if ( augustus_utr && not_fungus ) {
        strand_param = "--singlestrand=false --UTR=on"
    } else if ( !augustus_utr && not_fungus ) {
        strand_param = "--singlestrand=true --UTR=off"
    } else if ( !augustus_utr && !not_fungus ) {
        strand_param = "--singlestrand=true --UTR=off"
    } else if ( strand == "forward" ) {
        strand_param = "--strand=forward --UTR=on"
    } else if ( strand == "reverse" ) {
        strand_param = "--strand=backward --UTR=on"
    }

    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"

    augustus \
      --species="${augustus_species}" \
      --softmasking=on \
      ${strand_param} \
      --min_intron_len=${min_intron_hard} \
      --start=on \
      --stop=on \
      --introns=on \
      --cds=on \
      --gff3=on \
      --codingseq=on \
      --protein=on \
      --outfile="${name}_${strand}_augustus_denovo.gff" \
      --errfile=augustus.err \
      "${fasta}"
    """
}


/*
 * Run augustus with hints for each chunked genome.
 * If the genome is expected to have overlapping genes
 * and we are predicting UTRs, then we run each strand separately.
 * Without utrs, we can use the singlestrand parameter.
 */
process augustus_hints {

    label "augustus"
    label "small_task"
    time '1d'

    tag "${name} - ${strand}"

    input:
    val augustus_species
    val augustus_utr
    val not_fungus
    val min_intron_hard
    val valid_splicesites
    tuple val(name),
        val(strand),
        path(fasta),
        path("*hints")

    path "augustus_config"
    path "extrinsic.cfg"

    output:
    tuple val(name),
        val(strand),
        path("${name}_${strand}_augustus_hints.gff3")

    script:
    if ( augustus_utr && not_fungus ) {
        strand_param = "--singlestrand=false --UTR=on"
    } else if ( !augustus_utr && not_fungus ) {
        strand_param = "--singlestrand=true --UTR=off"
    } else if ( !augustus_utr && !not_fungus ) {
        strand_param = "--singlestrand=true --UTR=off"
    } else if ( strand == "forward" ) {
        strand_param = "--strand=forward --UTR=on"
    } else if ( strand == "reverse" ) {
        strand_param = "--strand=backward --UTR=on"
    }

    is_utr = augustus_utr ? "true" : "false"
    is_fungus = not_fungus ? "false" : "true"

    """
    export AUGUSTUS_CONFIG_PATH="\${PWD}/augustus_config"
    perl -n -e'/>(\\S+)/ && print \$1."\\n"' < "${fasta}" > seqids.txt

    if ${is_utr} && ${is_fungus}
    then
      # Gemoma doesn't do fungal utrs well.
      cat *hints \
      | awk -F '\t' '! (\$2 == "GeMoMa" && (\$3 == "exon" || \$3 == "UTRpart"))' \
      > hints.gff
    elif ${is_utr} && ! ${is_fungus}
    then
      cat *hints > hints.gff
    else
      cat *hints \
      | awk -F '\t' '\$3 != "exonpart" && \$3 != "exon" && \$3 != "UTRpart"' \
      > hints.gff
    fi

    getLinesMatching.pl seqids.txt 1 < hints.gff > hints_filtered.gff

    augustus \
      --species="${augustus_species}" \
      --extrinsicCfgFile=extrinsic.cfg \
      --hintsfile=hints_filtered.gff \
      ${strand_param} \
      --allow_hinted_splicesites="${valid_splicesites}" \
      --softmasking=on \
      --alternatives-from-evidence=true \
      --min_intron_len="${min_intron_hard}" \
      --start=off \
      --stop=off \
      --introns=off \
      --gff3=on \
      --outfile="out.gff3" \
      --errfile=augustus.err \
      "${fasta}"


    awk -F '\\t' '
      BEGIN {OFS="\\t"}
      \$3 == "transcript" {\$3="mRNA"}
      \$0 !~ /^#/ {print}
    ' out.gff3 \
    > "${name}_${strand}_augustus_hints.gff3"
    """
}


process evm {

    label "evm"
    label "big_task"
    time '1d'

    tag "${name}"

    input:
    val min_intron_hard
    path "weights.txt"
    tuple val(name),
        path(fasta),
        path(faidx),
        path("denovo.gff3"),
        path("transcripts.gff3"),
        path("proteins.gff3"),
        path("other.gff3")

    output:
    tuple val(name),
        path("${name}_evm.gff3")

    script:
    """
    partition_EVM_inputs.pl \
      --genome "${fasta}" \
      --gene_predictions <(cat denovo.gff3 other.gff3) \
      --protein_alignments proteins.gff3 \
      --transcript_alignments transcripts.gff3 \
      --segmentSize 500000 \
      --overlapSize 10000 \
      --partition partitions_list.out

    write_EVM_commands.pl \
      --genome "${fasta}" \
      --weights "\${PWD}/weights.txt" \
      --gene_predictions <(cat denovo.gff3 other.gff3) \
      --min_intron_length "${min_intron_hard}" \
      --protein_alignments proteins.gff3 \
      --transcript_alignments transcripts.gff3 \
      --output_file_name evm.out \
      --partitions partitions_list.out \
    > commands.list

    xargs -I{} -P${task.cpus} -- sh -c "{};" < commands.list

    recombine_EVM_partial_outputs.pl \
      --partitions partitions_list.out \
      --output_file_name evm.out

    convert_EVM_outputs_to_GFF3.pl \
      --partitions partitions_list.out \
      --output evm.out \
      --genome "${fasta}"

    find . -regex ".*evm.out.gff3" -exec cat {} \\; > "${name}_evm.gff3"
    """
}


/*
 * Select regions to re-predict genes in because EVM threw them out.
 *
 * If EVM found another gene at that locus we keep that one.
 * We filter out any hints that overlap CDS sequence of an EVM prediction
 * in the same strand, then we merge hints for these loci and pad the
 * region out a bit to give some context.
 */
process find_missing_evm_predictions {

    label "bedtools"
    label "small_task"
    time '3h'

    tag "${name}"

    input:
    tuple val(name),
        path("evm.gff3"),
        path("denovo.gff3"),
        path("transcripts.gff3"),
        path("proteins.gff3"),
        path("other.gff3"),
        path(faidx)

    output:
    tuple val(name), path("clustered.bed")

    script:
    """
    awk -F'\\t' '
      BEGIN { OFS="\\t" }
      \$3 == "mRNA" {
        print \$1, \$4, \$5, ".", ".", \$7
      }' other.gff3 \
    | sort \
      -k1,1 -k2,2n -k3,3n \
      --temporary-directory=tmp \
    | bedtools subtract \
        -a - \
        -b <(awk '\$3 == "CDS"' "evm.gff3") \
        -s \
        -A \
    | bedtools merge -s -c 6 -o distinct -i - \
    | bedtools slop -g "${faidx}" -b 5 -i - \
    > clustered.bed
    """
}


process augustus_gap_filler {

    label "augustus"
    label "big_task"
    time '1d'

    tag "${name}"

    input:
    val augustus_species
    val not_fungus
    val augustus_utr
    val valid_splicesites
    val min_intron_hard
    tuple val(name),
        path(fasta),
        path("toredo.bed"),
        path("*hints")

    path "augustus_config" from augustusConfig
    path "extrinsic.cfg" from augustusGapFillerWeights

    output:
    tuple val(name), path("augustus_gaps/*.gff3")

    script:
    is_utr = augustus_utr ? "true" : "false"
    is_fungus = not_fungus ? "false" : "true"
    utr_flag = augustus_utr ? "-u" : ""

    """
    if ${is_utr} && ${is_fungus}
    then
      # Gemoma doesn't do fungal utrs well.
      cat *hints \
      | awk -F '\t' '! ((\$2 == "GeMoMa" || \$2 == "ComparativeGeMoMa") &&
                        (\$3 == "exon" || \$3 == "UTRpart"))' \
      > hints.gff
    elif ${is_utr} && ! ${is_fungus}
    then
      cat *hints > hints.gff
    else
      cat *hints \
      | awk -F '\t' '\$3 != "exonpart" && \$3 != "exon" && \$3 != "UTRpart"' \
      > hints.gff
    fi

    augustus_region.sh \
      -f "${fasta}" \
      -b "toredo.bed" \
      -g "hints.gff" \
      -s "${augustus_species}" \
      -c "extrinsic.cfg" \
      -a "\${PWD}/augustus_config" \
      -p "${valid_splicesites}" \
      ${utr_flag} \
      -n "${task.cpus}" \
      -m "${min_intron_hard}" \
      -o "augustus_gaps"
    """
}

