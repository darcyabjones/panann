# panann

A pipeline for gene annotation of pangenomes.
We have a particular focus on Fungal pathogen genomes, but the principles could be adapted to other systems.

This is very-much a work in progress and mostly developed for personal use.
I try to make things work well for general use, but some parameters may be highly specific to the organisms and datasets I am actually using.

See the [wiki](https://github.com/darcyabjones/panann/wiki) for details, research, and updates.

NOTE: As of 2019-10-10 i have not yet put the containers on singularity or docker hub so some of the install instructions wont work.
I'll add these in the coming days, in the meantime you can still follow the instructions to build the containers yourself.

## What does this do?

Panann combines comparative, de-novo, and evidence based gene prediction methods in a single pipeline.
It is indended for annotating large numbers of genomees from a single organism.
We focus on sensitivity rather than specificity, because filtering by comparative genomics is possible (e.g. Dn/Ds, frequency across the population etc) and because I work with fungal pathogens and our genes of interest frequently get missed by gene predictions.

To run panann you'll need:

- Stranded paired end illumina RNA-seq for one of your isolates.
  Pre-aligned/assembled Iso-seq transcripts may also work.
- A well trained augustus model for your species of interest (or one that performs well from another species).
  We give strong weights to augustus predictions, so this is important.
- A set of proteins from closely related organisms (e.g. same species or genus).
- A number of already soft-masked fasta genomes.
  You can use repeat-masker to soft-mask your genomes.
  Ideally, at least one of these genomes should already have a high-quality set of annotations,
  though this is not required.


The pipeline follows these main steps:

1. Align Illumina RNA-seq reads, and full length transcripts to the genome, optionally assembling the transcripts from fastq reads using Trinity.
2. Align proteins from closely related organisms to the genome, and align remotely related proteins (e.g. uniref90) to the genomes.
3. Predict genes and transfer known annotations using Augustus, CodingQuarry, GeMoMa, PASA/TransDecoder, and GeneMark-ES.
4. Collect proteins from all predicted annotations in all genomes, cluster them, and transfer representative members of each cluster to each genome using GeMoMa.
5. Combine the comparative, de-novo, and evidence based predictions using EvidenceModeler.
6. Find loci that were excluded from EVM output (e.g. because it doesn't like non-standard splicesites), and use all predictions as hints for augustus to predict genes in these regions.
7. Combine EVM and the gap-filling augustus output into a single final prediction set.
8. Collect statistics about each annotation set, including BUSCO completeness, gene-counts, gene-lengths, splice-site frequencies, and concordance with known annotation sets.


Not implemented but on the wish-list:

- Filtering false positives via Dn/Ds, Frequency, RNAseq coverage, and matches to databases.
- Tidying up GFFs and adding information to the attributes column.
- Summarising statistics across the population.


## How not to use this software

There are many pipelines for doing genome annotation.
Even fungal genome annotation!

Everyones needs are slightly different.

Here are some honourable mentions:

- [Funannotate](https://funannotate.readthedocs.io/en/latest/) great for single genomes. Heavy reliance on PASA/Transdecoder, which performed poorly on our datasets.
- [MAKER](https://www.yandell-lab.org/software/maker.html) seems to be the top choice for people that just want to get it done quickly.
- [LoReAn](https://github.com/lfaino/LoReAn) uses long reads for annotations.
- [Comparative annotation toolkit](https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit) is a neat pan-genome annotation pipeline. It is heavily based around augustus and uses the [CGP](https://github.com/Gaius-Augustus/Augustus/blob/master/README-cgp.md) mode for transferring annotations. The run-time of the [Cactus genome aligner](https://github.com/ComparativeGenomicsToolkit/cactus) may be restrictive for you if you have many genomes or very large genomes.


## Quick start

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

nextflow run -profile singularity darcyabjones/panann -resume \
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



## Install

The pipeline is written in [Nextflow](https://www.nextflow.io/), which you will need to install on the executing computer.
The pipeline itself has many dependencies, so I recommend that you use the containers.
To run the containers, you'll need to install either [Singularity](https://sylabs.io/docs/) (recommended) or [Docker](https://www.docker.com/).

The pipeline itself will pull the containers for you.

On an ubuntu server, the process to install nextflow and singularity might look like this.

```bash
set -eu

sudo apt-get update
sudo apt-get install -y \
    default-jre-headless \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    wget \
    pkg-config \
    git

VERSION=1.12
OS=linux
ARCH=amd64
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
rm go$VERSION.$OS-$ARCH.tar.gz

echo 'export PATH=/usr/local/go/bin:$PATH' >> ~/.bashrc
source ~/.bashrc

VERSION=3.4.0
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz
tar -xzf singularity-${VERSION}.tar.gz
cd singularity
./mconfig
make -C builddir
sudo make -C builddir install

cd ..
rm -rf -- singularity

curl -s https://get.nextflow.io | bash

./nextflow run darcyabjones/panann --help
```

Because [GeneMark-ES](http://exon.gatech.edu/GeneMark/gmes_instructions.html) and [SignalP](http://www.cbs.dtu.dk/services/SignalP/) both have restricted licenses, they are not available in the regular containers.
You can instead build the containers yourself if you provide the source files.
See the [containers](containers/README.md) directory for more.

To tell panann to use these tools, provide the `--genemark` or `--signalp` flags.
If you specify those and the software isn't installed or a container doesn't exist, the pipeline will unceremoniously fail.
Without `--signalp`, panann will instead use [DeepSig](https://github.com/BolognaBiocomp/deepsig) to train codingquarry pathogen-mode, which is a perfectly fine alternative.
You can also disable codingquarry pathogen mode with the `--notfungus` or `--nocqpm` options and neither signalp or deepsig will be run.
Without `--genemark` panann just won't run genemark and the other predictors will probably pick up the slack.


## Profiles

The strength of pipeline engines like nextflow is that you can run it on different compute systems
simply by switching some configuration files.

Some preset config files are included in this repo.
You can view these config files in the `conf` directory.

The configuration to use at runtime is controlled by the `-profile` parameter.

Multiple profiles can be specified by separating them with a comma e.g. `-profile laptop,singularity`.
Panann generally has a separate config file for a compute environment (e.g. cloud, HPC, laptop), and for a software environment (e.g. singularity, docker, local).
It's likely that you'll have to tailor the compute configuration, but you shouldn't need to change the software config so this allows you to mix-and-match.

For more info on configuration see the [nextflow documentation](https://www.nextflow.io/docs/latest/config.html).
You can also raise an issue on the github repository and I'll try to help.


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
| `--spaln_species` | phaenodo | The species specific parameters to use with the spaln transcript aligner. Options are available here <https://github.com/ogotoh/spaln/tree/master/table> (one of the directory names). Selecting a good parameterset dramatically improves alignment quality, so trial some if you have the time. |
| `--trans_table` | 1 | The ncbi translation table number to use. Note that not all tools currently respect this parameter, so organisms with significantly different tables to standard might not work well.|
| `--notfungus` | false | Disable some of the tuning that is done to support organisms with high gene density. Also disable codingquarry analysis. |
| `--genemark` | false | Specify that you want to run genemark and have it installed. Because genemark has licensing issues for non-academic institutions, you may not be able to use it. |
| `--notrinity` | false | Don't assemble the RNAseq reads even if you provided `--fastq`. |
| `--nostar` | false | Don't align RNAseq reads even if you provided `--fastq`. |
| `--nocqpm` | false | Don't run codingquarry pathogen mode. Use this if you are unlikely to have secreted genes that look significantly different to the rest of the genes. |
| `--training` | Optional | Only predict genes in the genome with this name. Other genomes provided with known sites etc will still be transferred across. Use this if you want to train augustus or generate a set of hints to optimise your extrinsic parameters config. |
| `--outdir` | results | Where to store the results. |


# Output

You can specify the main output directory using the `--outdir` parameter.
Here we use <outdir> inplace of that parameter and <name> inplace of the
genome basenames.

<outdir>/aligned/<name>/<name>_* -- Aligned reads, transcripts, and proteins to each genome.
<outdir>/annotations/<name>/<name>_* -- Gene predictions for each genome, including output of individual methods.
<outdir>/hints/<name>/<name>_*.gff3 -- GFF files from alignments and annotations configured to be provided as hints to augustus.
<outdir>/assembled/<read_group>_* -- Trinity assembled reads.
<outdir>/qc/<name>_* -- Statistics and QC for each genome and step.
<outdir>/trace -- Information about the nextflow run (e.g. memory use, runtime etc).

# Exit codes

I'm hoping to add some better error-handling in the future to provide more useful/nextflow-agnostic tips to users.
In the meantime, it's just input parameter validation that is handled elegantly.

- 0: All ok.
- 1: Incomplete parameter inputs.
