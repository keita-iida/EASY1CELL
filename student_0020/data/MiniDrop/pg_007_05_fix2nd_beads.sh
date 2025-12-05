#!/bin/sh

DIR1=007_bam
INCOMING=${DIR1}/pg_007_04_tiny1.clean_substitution.bam
OUTGOING=${DIR1}/pg_007_05_tiny1.clean.bam
#-----------------------------------------------------------------------------80
# Fixed bead error (composite).
#-----------------------------------------------------------------------------80
DetectBeadSynthesisErrors \
  I=${INCOMING} \
  O=${OUTGOING} \
  REPORT=${OUTGOING}.clean.indel_report.txt \
  OUTPUT_STATS=${OUTGOING}.synthesis_stats.txt \
  SUMMARY=${OUTGOING}.synthesis_stats.summary.txt \
  PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
  TMP_DIR=${DIR1}/tmp

samtools view -h ${OUTGOING} > ${OUTGOING}.sam
samtools index ${OUTGOING}

