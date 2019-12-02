#!/usr/bin/env nextflow

params.genomes = false
params.stringtie = false

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

genomes.into {
    genomes4RunCodingQuarry;
    genomes4RunCodingQuarryPM;
}


if ( params.stringtie ) {
    Channel
        .fromPath(params.stringtie, checkIfExists: true, type: 'file')
        .splitCsv(by: 1, sep: '\t', header: true)
        .filter { (!is_null(it.name) && !is_null(it.stringtie) }
        .map {[it.name, file(it.stringtie, checkIfExists: true]}
        .unique()
        .into { stringtie4CodingQuarry; stringtie4CodingQuarryPM }
} else {
    log.error "Running codingquarry requires stringtie input."
    exit 1
}


/*
 * Predict genes using codingquarry
 * Note that we still need the proteins extracted from here for codingquarrypm.
 */
process runCodingQuarry {

    label "codingquarry"
    label "bigmem_task"
    time '1d'
    errorStrategy "retry"
    maxRetries 10

    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    input:
    set val(name),
        file("transcripts.gtf"),
        file("genome.fasta") from stringtie4CodingQuarry
            .combine(genomes4RunCodingQuarry, by: 0)

    output:
    set val(name), file("${name}_codingquarry.gff3") into codingQuarryPredictions
    set val(name), file("${name}_codingquarry.faa") into codingQuarryPredictionsProteins
    file "${name}_codingquarry.fna"
    file "${name}_codingquarry_dubiousset.gff3"
    file "${name}_codingquarry_fusions.txt"
    file "${name}_codingquarry_overlapreport.txt"

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
    mv out/Predicted_CDS.fa "${name}_codingquarry.fna"
    mv out/Predicted_Proteins.faa "${name}_codingquarry.faa"
    mv out/fusions.txt "${name}_codingquarry_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarry_overlapreport.txt"

    rm -rf -- out
    """
}

codingQuarryPredictions.into {
    codingQuarryPredictions4SecretionPred;
    codingQuarryPredictions4PM;
    codingQuarryPredictions4Tidy;
}


/*
 * Fix weird cq gffs and get introns.
 */
process tidyCodingQuarry {

    label "aegean"
    label "small_task"
    time '2h'

    tag "${name}"

    publishDir "${params.outdir}/annotations/${name}"

    input:
    set val(name), file("codingquarry.gff3") from codingQuarryPredictions4Tidy

    output:
    set val(name),
        val("codingquarry"),
        file("${name}_codingquarry_tidy.gff3") into tidiedCodingQuarry

    script:
    """
    awk -F '\t' 'BEGIN {OFS="\\t"} \$8 = "-1" {\$8="0"} {print}' codingquarry.gff3 \
    | awk -F '\t' 'BEGIN {OFS="\\t"} \$3 == "gene" {\$3="mRNA"} {print}' \
    | gt gff3 -tidy -sort -retainids -setsource 'CodingQuarry' \
    | canon-gff3 -i - \
    > "${name}_codingquarry_tidy.gff3"
    """
}


/*
 * Predict signal peptides from CodingQuarry first pass predicted proteins.
 */
if (params.signalp) {

    process getCodingQuarrySignalP {

        label "signalp"
        label "medium_task"
        time '12h'

        tag "${name}"

        when:
        !params.nocqpm

        input:
        set val(name),
            file("proteins.faa") from codingQuarryPredictionsProteins

        output:
        set val(name),
            file("${name}_codingquarry_secreted.txt") into codingQuarryPredictionsSecreted

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
        > "${name}_codingquarry_secreted.txt"
        """
    }

} else {

    process getCodingQuarryDeepsig {

        label "deepsig"
        label "medium_task"
        time '12h'

        publishDir "${params.outdir}/annotations/${name}"

        tag "${name}"

        when:
        !params.nocqpm

        input:
        set val(name),
            file("proteins.faa") from codingQuarryPredictionsProteins

        output:
        set val(name),
            file("${name}_codingquarry_secreted.txt") into codingQuarryPredictionsSecreted

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
        ' < secreted.txt > "${name}_codingquarry_secreted.txt"
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
 *
 * Try setting `vm.overcommit_memory = 1` if you get segfaults before
 * reaching max memory.
 *
 * We also set the output to be optional because it requires > 500 proteins to train from.
 */
process runCodingQuarryPM {

    label "codingquarry"
    label "bigmem_task"
    time '1d'

    errorStrategy "retry"
    maxRetries 10

    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    when:
    !params.nocqpm

    input:
    set val(name),
        file("transcripts.gtf"),
        file("genome.fasta"),
        file("codingquarry.gff3"),
        file("secretome.txt") from stringtie4CodingQuarryPM
            .combine(genomes4RunCodingQuarryPM, by: 0)
            .combine(codingQuarryPredictions4PM, by: 0)
            .combine(codingQuarryPredictionsSecreted, by: 0)

    output:
    set val(name),
        file("${name}_codingquarrypm.gff3") optional true into codingQuarryPMPredictions

    file "${name}_codingquarrypm_fusions.txt" optional true
    file "${name}_codingquarrypm_overlapreport.txt" optional true

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

    mv out/PGN_predictedPass.gff3 "${name}_codingquarrypm.gff3"
    mv out/fusions.txt "${name}_codingquarrypm_fusions.txt"
    mv out/overlapReport.txt "${name}_codingquarrypm_overlapreport.txt"
    rm -rf -- out
    """
}


/*
 * Fix weird CodingQuarry gffs and add introns etc.
 */
process tidyCodingQuarryPM {

    label "aegean"
    label "small_task"
    time '2h'

    publishDir "${params.outdir}/annotations/${name}"

    tag "${name}"

    input:
    set val(name),
        file("codingquarry.gff3") from codingQuarryPMPredictions

    output:
    set val(name),
        val("codingquarrypm"),
        file("${name}_codingquarrypm_tidy.gff3") into tidiedCodingQuarryPM

    script:
    """
      awk -F '\t' 'BEGIN {OFS="\\t"} \$8 == "-1" {\$8="0"} {print}' codingquarry.gff3 \
    | awk -F '\t' 'BEGIN {OFS="\\t"} \$3 == "gene" {\$3="mRNA"} {print}' \
    | gt gff3 -tidy -sort -retainids -setsource "CodingQuarryPM" \
    | canon-gff3 -i - \
    > "${name}_codingquarrypm_tidy.gff3"
    """
}
