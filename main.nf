#!/usr/bin/env nextflow


def help_message() {
    log.info"""
    # panann

    A pipeline to predict genes in a multiple fungal genomes.

    ## Usage

    Say you have several genomes to annotate and an already well annotated genome
    "genome1". I also assume that you have genemark installed (See INSTALL).

    ```bash
    cat <<EOF > known_sites.tsv
    name	gff3
    genome1	genome1.gff3
    EOF

    cat <<EOF > rnaseq.tsv
    read_group	r1	r2
    inplanta	fastq/ip_r1.fq.gz	fastq/ip_r2.fq.gz
    inplanta	fastq/ip2_r1.fq.gz	fastq/ip2_r2.fq.gz
    invitro	fastq/iv_r1.fq.gz	fastq/iv_r2.fq.gz
    EOF

    nextflow run -profile singularity main.nf -resume \
      --genomes "genomes/*.fasta" \
      --transcripts 'transcripts/*.fasta' \
      --proteins 'proteins/*.faa' \
      --remote_proteins ./uniref-identity%3A0.9.fasta \
      --fastq rnaseq.tsv \
      --known_sites known_sites.tsv \
      --augustus_species "parastagonospora_nodorum_sn15" \
      --augustus_config ./my_augustus_config \
      --genemark
    ```

    ## Parameters

    | parameter | Default | Description |
    | :---      | :---    | :---        |
    | `--genomes` | Required | A glob of the fasta genomes to search for genes in. The genome should already be soft-masked. The basename of the file is used as the genome name. |
    | `--known_sites` | Optional | Specify existing high-confidence genes to use as manual hints with augustus, and for annotation transfer with GeMoMa. A tab-separated table with a header containing the columns `name` and `gff3`. The gff3 column should contain the full path to the gff3 file, and the name should match the basename of one of the fasta files. |
    | `--transcripts` | Optional | A glob to fasta files containing ESTs/transcripts to align to the genome and use with PASA. If assembling fastq reads with trinity, both will be passed to aligners as a combined file. |
    | `--proteins` | Optional | A glob of fasta files containing Proteins to use as hints for augustus and EVM. This is intended for aligning proteins from very closely related genomes (i.e. same species). |
    | `--remote_proteins` | Optional | A glob to fasta file containing proteins from diverse sources to align to the genome and use as hints for augustus and EVM. Intended for using uniref90 or similar. |
    | `--remote_alignments` | Optional | A tab separated file containing precomputed alignments of remote proteins to the genomes. The file should have a header, and the columns `name` and `remote_alignments`. The name should match the basename of one of the `genome` fasta files. `remote_alignments` should be the full path to an Exonerate GFF2 results file. If this is specified, anything in `--remote_proteins` will not be aligned and this will be used instead even if no alignments exist for a particular genome.|
    | `--busco_lineage` | Optional | The path to the directory containing busco lineage information. If provided, different prediction methods' completeness will be evaluated. The folder can be downloaded from <https://busco.ezlab.org/> under datasets. Unzip the folder. |
    | `--busco_genomes` | false | A flag specifying that you want to run busco on the genome sequence as well as the protein completeness tests. Requires that `--busco_lineage` is set. |
    | `--augustus_config` | environment variable `AUGUSTUS_CONFIG_PATH` in the augustus environment | A path to the directory containing augustus config information. This is the `config` directory in the git repository or where-ever you installed augustus <https://github.com/Gaius-Augustus/Augustus>. Use this to provide custom trained augustus species (recommended). |
    | `--augustus_species` | Required | The name of the augustus model to use for gene prediction. This should correspond to one of the folders in `\${AUGUSTUS_CONFIG_PATH}/species`. |
    | `--augustus_utr` | false | A flag specifying that augustus should be run with the `--UTR` flag. Requires that the species model has been trained with UTR parameters. Not recommended for organisms with high-gene density like fungi, as it tends to exclude genes on the same strand that are very close together. |
    | `--augustus_hint_weights` | data/extrinsic_hints.cfg | The extrinsic config hints file to use for Augustus prediction. |
    | `--augustus_gapfiller_weights` | data/extrinsic_gapfiller.cfg | The extrinsic config file to use for filling in gaps in EVM predictions with Augustus. |
    | `--evm_config` | data/evm.cfg | The evidence modeller weights that you would like to use. |
    | `--fastq` | Optional if `--crams` is provided. | A tab-separated table containing the fastq reads to align to the genome and use as hints. The table should contain a header with the columns `read_group`, `r1`, and `r2`. `r1` and `r2` are the paths to gzipped fastq files. Only illumina paired-end stranded reads are supported. The `read_group` should be specified for each treatment, and each `read_group` may have multiple fastq pairs on separate lines." |
    | `--crams` | Optional | A tab-separated table containing pre-aligned RNAseq reads to use instead of `--fastq`. The table should contain the columns `name`, `read_group`, and `cram`. The name should correspond to the basename of one of the fasta genomes. The read group is as for `--fastq`. The `cram` column should be the path to a cram file <http://samtools.github.io/hts-specs/> that were created using the genome fasta with `name` as a reference. |
    | `--fr` | false | Flag telling tools that the stranded RNAseq library was created with FR rather than the usual illumina protocols RF. |
    | `--valid_splicesites` | gtag,gcag,atac,ctac | A comma separated list of valid splice sites to allow. Any pairs not in this list will be filtered from the augustus predictions and from training sets for genemark. Including lots of sites can result in spurious introns, so evaluate your genome(s) first. |
    | `--min_intron_soft` | 20 | The minimum intron length to consider generally. Tools that use this parameter still allow smaller introns, as this is more of a shape parameter than a hard limit. |
    | `--min_intron_hard` | 5 | The absolute minimum intron length to consider. Anything smaller than this will be filtered out. |
    | `--max_intron_hard` | 15000 | The absolute maximum intron length to consider. Any bigger introns than this will be filtered out. |
    | `--star_novel_params` | See code (sorry) | Parameters to pass to the STAR aligner for the first stage of RNAseq read alignment (finding novel splice sites). These should generally be splice junction filter options. |
    | `--star_align_params` | See code (sorry) | Parameters to pass to the STAR aligner for the second stage of RNAseq read alignment (aligning around the splice-sites). These should generally be splice junction filter options. Reads that don't pass these filters will be excluded from output. |
    | `--max_gene_hard` | 20000 | The absolute maximum length in basepairs that a gene can be (including introns). |
    | `--spaln_species | phaenodo | The species specific parameters to use with the spaln transcript aligner. Options are available here <https://github.com/ogotoh/spaln/tree/master/table> (one of the directory names). Selecting a good parameterset dramatically improves alignment quality, so trial some if you have the time. |
    | `--trans_table` | 1 | The ncbi translation table number to use. Note that not all tools currently respect this parameter, so organisms with significantly different tables to standard might not work well.|
    | `--notfungus` | false | Disable some of the tuning that is done to support organisms with high gene density. Also disable codingquarry analysis. |
    | `--genemark` | false | Specify that you want to run genemark and have it installed. Because genemark has licensing issues for non-academic institutions, you may not be able to use it. |
    | `--notrinity` | false | Don't assemble the RNAseq reads even if you provided `--fastq`. |
    | `--nostar` | false | Don't align RNAseq reads even if you provided `--fastq`. |
    | `--nocqpm` | false | Don't run codingquarry pathogen mode. Use this if you are unlikely to have secreted genes that look significantly different to the rest of the genes. |
    | `--training` | Optional | Only predict genes in the genome with this name. Other genomes provided with known sites etc will still be transferred across. Use this if you want to train augustus or generate a set of hints to optimise your extrinsic parameters config. |
    | `--outdir` | results | Where to store the results. |


    ## Output

    You can specify the main output directory using the `--outdir` parameter.
    Here we use <outdir> inplace of that parameter and <name> inplace of the
    genome basenames.

    <outdir>/aligned/<name>/<name>_* -- Aligned reads, transcripts, and proteins to each genome.
    <outdir>/annotations/<name>/<name>_* -- Gene predictions for each genome, including output of individual methods.
    <outdir>/hints/<name>/<name>_*.gff3 -- GFF files from alignments and annotations configured to be provided as hints to augustus.
    <outdir>/assembled/<read_group>_* -- Trinity assembled reads.
    <outdir>/qc/<name>_* -- Statistics and QC for each genome and step.
    <outdir>/trace -- Information about the nextflow run (e.g. memory use, runtime etc).

    ## Exit codes

    - 0: All ok.
    - 1: Incomplete parameter inputs.
    """.stripIndent()
}


workflow {

    main:
    if ( params.help ) {
        helpMessage()
        exit 0
    }

    def should_run(validate_params)
}
