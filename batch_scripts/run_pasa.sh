#!/usr/bin/env bash

MAX_INTRON_HARD=15000
GENOMES=()

INDEX=$(( ${SLURM_ARRAY_TASK_ID:-0} + ${SLURM_PROCID:-0} ))

NAME="${GENOMES[${INDEX}]}"

RESULTS_DIR="${PWD}/results"
GENOME_DIR="${PWD}/input/stago"

GENOME="${GENOME_DIR}/${NAME}.fasta"
GMAP_RESULTS="${RESULTS_DIR}/aligned/${NAME}/${NAME}_gmap_transcripts.gff3"
TRANSCRIPTS="${RESULTS_DIR}/transcripts/transcripts/transcripts.fasta"
TRANSCRIPTS_CLEAN="${TRANSCRIPTS}.clean"

OUTDIR="${PWD}/pasa_sep"
TMPDIR="${PWD}/.pasa_tmp_$$"

mkdir -p "${TMPDIR}"
cd "${TMPDIR}"

echo "DATABASE=${PWD}/pasa.sqlite" > align_assembly.config
echo "validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80" >> align_assembly.config
echo "validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=90" >> align_assembly.config
echo "subcluster_builder.dbi:-m=50" >> align_assembly.config

Launch_PASA_pipeline.pl \
  --config align_assembly.config \
  --create \
  --run \
  --genome "${GENOME}" \
  --transcripts "${TRANSCRIPTS_CLEAN}" \
  --IMPORT_CUSTOM_ALIGNMENTS_GFF3 "${GMAP_RESULTS}" \
  -T -u "${TRANSCRIPTS}" \
  --MAX_INTRON_LENGTH "${MAX_INTRON_HARD}" \
  --ALIGNERS blat \
  --CPU 1 \
  --transcribed_is_aligned_orient \
  --TRANSDECODER \
  --stringent_alignment_overlap 30.0

pasa_asmbls_to_training_set.dbi \
  -G "Universal" \
  --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta \
  --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3

cp -L pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 "${OUTDIR}/${NAME}_pasa.gff3"
