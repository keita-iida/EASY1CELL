#!/bin/sh

DIR1=005_bam
INCOMING=${DIR1}/pg_005_04_tiny1.cb.umi.qc.bam
OUTGOING=${DIR1}/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam
#-----------------------------------------------------------------------------80
# Remove potential SMART adapter sequence (start site trimming).
#-----------------------------------------------------------------------------80
TrimStartingSequence \
  INPUT=${INCOMING} \
  OUTPUT=${OUTGOING} \
  OUTPUT_SUMMARY=${OUTGOING}.report.txt \
  SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG MISMATCHES=0 NUM_BASES=5

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

