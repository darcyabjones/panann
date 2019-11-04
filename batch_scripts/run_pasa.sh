#!/usr/bin/env bash

GENOMES=()

INDEX=$(( ${SLURM_ARRAY_TASK_ID:-0} + ${SLURM_PROCID:-0} ))

OUTDIR="${PWD}/pasa_sep"
TMPDIR=".pasa_tmp_$$"
mkdir -p "${TMPDIR}"

echo "DATABASE=${PWD}/pasa.sqlite" > align_assembly.config
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
  --MAX_INTRON_LENGTH "${params.max_intron_hard}" \
  --ALIGNERS blat \
  --CPU 1 \
  --transcribed_is_aligned_orient \
  --TRANSDECODER \
  --stringent_alignment_overlap 30.0 \ 
  ${use_known}

    pasa_asmbls_to_training_set.dbi \
      -G "${gen_code}" \
      --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta \
      --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3

    ln -s \${PWD}/pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 \${PWD}/${name}_pasa.gff3
