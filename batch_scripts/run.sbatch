#!/bin/bash --login
 
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00
#SBATCH --account=y95
#SBATCH --export=NONE

module load singularity/3.3.0
module load nextflow/19.01.0.5050-bin


nextflow run -profile pawsey_zeus -resume -with-singularity ./panann-plus.sif ./main.nf \
  --genomes "input/stago/*.fasta" \
  --transcripts "input/transcripts/*.fasta" \
  --proteins "input/proteins/*.faa" \
  --remote_proteins "input/uniref0.5fungi_20191029.fasta" \
  --augustus_config "input/augustus_config" \
  --augustus_species "parastagonospora_nodorum_sn15" \
  --busco_lineage "input/pezizomycotina_odb9" \
  --crams "input/cram_config.tsv" \
  --known_sites "input/known_sites_config.tsv" \
  --spaln_species "phaenodo" \
  --genemark \
  --signalp \
  --notrinity \
  --nostar \
  --outdir "results"
