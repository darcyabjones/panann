#!/usr/bin/env bash

set -euo pipefail

MAX_INTRON_HARD=15000
GENOMES=( 15FG04 15FG101 15FG103 15FG104 15FG105 15FG107 15FG109 15FG110 15FG112 15FG113 15FG114 15FG116 15FG117 15FG118 15FG119 15FG120 15FG226 15FG229 15FG28 15FG37 15FG47 15FG49 15FG99 16FG06 16FG158 16FG160 16FG162 16FG164 16FG165 16FG166 16FG167 16FG168 16FG169 16FG170 201FG208 201FG211 201FG218 201FG219 201FG49 202FG212 202FG414 203FG213 203FG58 204FG214 204FG221 204FG222 204FG223 205FG142 205FG215_1 205FG215_2 205FG216 205FG225 205FG410 206FG226 206FG227 206FG66 206FG67 206FG68 53FG143_1 53FG143_2 903FG214 FG106 FG107 FG108 Gerald1 Gerald4 Meck1 Meck3 Meck6 Meck8 Nor_RAC2182_1 Northam_Emu1 Northam_Mace1 Northam_Mace2 Northam_Magenta Northam_WGT RSID01 RSID02 RSID03 RSID04 RSID25 RSID28 RSID30 RSID31 RSID33 RSID35 RSID36 RSID37 RSID39 RSID42 S1FT3B S4FT3A SN2000 SN4 SN79 WAC13068 WAC13069 WAC13070 WAC13071 WAC13072 WAC13073 WAC13074 WAC13075 WAC13076 WAC13077 WAC13402 WAC13403 WAC13405 WAC13443 WAC13447 WAC13524 WAC13525 WAC13526 WAC13527 WAC13528 WAC13530 WAC13532 WAC13616 WAC13617 WAC13630 WAC13631 WAC13632 WAC13690 WAC13955 WAC2285 WAC2810 WAC2813 WAC4303 WAC4321 WAC4648 WAC4808 WAC739 WAC740 WAC741 WAC8384 WAC8390 WAC8410 WAC8635 WAC9178 )

INDEX=$(( (${SLURM_ARRAY_TASK_ID:-0} * 28) + ${SLURM_PROCID:-0} ))

if [ ${INDEX} -gt ${#GENOMES[@]} ]
then
   echo 0
fi

NAME="${GENOMES[${INDEX}]}"

echo "STARTING: ${NAME}"

RESULTS_DIR="${PWD}/results"
GENOME_DIR="${PWD}/input/stago"

GENOME="${GENOME_DIR}/${NAME}.fasta"
GMAP_RESULTS="${RESULTS_DIR}/aligned/${NAME}/${NAME}_gmap_transcripts.gff3"
TRANSCRIPTS="${RESULTS_DIR}/transcripts/transcripts.fasta"
TRANSCRIPTS_CLEAN="${TRANSCRIPTS}.clean"

OUTDIR="${PWD}/pasa_sep"
TMPDIR="${PWD}/.pasa_tmp_${NAME}"

mkdir -p "${TMPDIR}"
cd "${TMPDIR}"

#echo "DATABASE=${PWD}/pasa.sqlite" > align_assembly.config
#echo "validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80" >> align_assembly.config
#echo "validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=90" >> align_assembly.config
#echo "subcluster_builder.dbi:-m=50" >> align_assembly.config

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
  --stringent_alignment_overlap 30.0 \
> pasa_pipeline.log


echo "PASA: ${NAME}"

pasa_asmbls_to_training_set.dbi \
  -G "Universal" \
  --pasa_transcripts_fasta pasa.sqlite.assemblies.fasta \
  --pasa_transcripts_gff3 pasa.sqlite.pasa_assemblies.gff3 \
> transdecoder.log

cp -L pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 "${OUTDIR}/${NAME}_pasa.gff3"

echo "FINISHED: ${NAME}"
