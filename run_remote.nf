#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/


def helpMessage() {
    log.info"""
    # panann


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
params.min_contig_length = 500
params.min_intron_soft = 20
params.min_intron_hard = 5
params.max_intron_hard = 15000
params.max_gene_hard = 20000
params.genomes = false

// Remote proteins aligned to each genome with exonerate.
// As a tsv file, with name and remote columns
params.remote_alignments = false

params.trans_table = 1

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


if ( params.remote_proteins ) {
    Channel
        .fromPath(params.remote_proteins, checkIfExists: true, type: "file")
        .first()
        .set { remoteProteins }

} else {
    log.error "Please provide remote proteins"
    exit 1
}




process getFaidx {

    label "samtools"
    label "small_task"
    time '1h'

    tag "${name}"

    input:
    set val(name), file("orig.fa") from genomes

    output:
    set val(name), file("${name}.fasta"), file("${name}.fasta.fai") into genomesWithFaidx

    script:
    """
    # braker panics if the genome has descriptions
    sed -r 's/^(>[^[:space:]]*).*\$/\\1/' orig.fa \
    | fasta_to_tsv.sh \
    | awk 'length(\$2) >= ${params.min_contig_length}' \
    | tsv_to_fasta.sh \
    > "${name}.fasta"


    samtools faidx "${name}.fasta"
    """
}

genomesWithFaidx
    .into {
        genomes4MatchRemoteProteinsToGenome;
        genomes4ClusterRemoteProteinsToGenome;
        genomes4AlignRemoteProteinsToGenome;
    }


/*
 * Preindex MMSeqs remote proteins
 */
process indexRemoteProteins {

    label "mmseqs"
    label "small_task"
    time '3h'

    input:
    file fasta from remoteProteins.collectFile(newLine: true)

    output:
    file "proteins" into indexedRemoteProteins
    file "proteins.tsv" into remoteProteinsTSV

    script:
    """
    mkdir -p tmp proteins

    mmseqs createdb "${fasta}" proteins/db
    # mmseqs createindex proteins/db tmp --threads "${task.cpus}"

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
 * Quickly find genomic regions with matches to remote proteins.
 * This is an approximate method which we refine later with exonerate.
 */
process matchRemoteProteinsToGenome {

    label "mmseqs"
    label "big_task"
    time '6h'

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
    cp -rL proteins proteins_tmp

    # the "dont-split" bit is important for keeping the ids correct.
    mmseqs createdb \
      "${fasta}" \
      genome/db \
      --dont-split-seq-by-len

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
      --translation-table "${params.trans_table}" \
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
    > "${name}_remote_proteins.tsv"

    sed -i '1i query\ttarget\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tpident\tmismatch\tgapopen\tevalue\tbitscore' "${name}_remote_proteins.tsv"

    rm -rf -- tmp genome result results_unsorted.tsv proteins_tmp
    """
}


/*
 * Finds regions genome with high density of remote protein matches.
 * We align remote proteins to each of these regions individually.
 */
process clusterRemoteProteinsToGenome {

    label "bedtools"
    label "medium_task"
    time '3h'

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
    def exonerate_buffer_region = params.max_gene_hard

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
    time '6h'

    publishDir "${params.outdir}/remote_hints"
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
        file("${name}_remote_proteins_exonerate.gff") into tmpAlignedRemoteProteins

    script:
    """
    mkdir -p tmp
    cp -L proteins.tsv proteins_tmp.tsv

    exonerate_parallel.sh \
      -g "${fasta}" \
      -q "proteins_tmp.tsv" \
      -b "clustered.bed" \
      -n "${task.cpus}" \
      -t "tmp" \
      -m "${params.min_intron_hard}" \
      -x "${params.max_intron_hard}" \
      -r "${params.trans_table}" \
      -o "${name}_remote_proteins_exonerate.gff"

    rm -rf -- tmp proteins_tmp.tsv
    """
}


