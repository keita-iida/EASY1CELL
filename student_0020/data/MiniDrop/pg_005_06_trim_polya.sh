#!/bin/sh

DIR1=005_bam
INCOMING=${DIR1}/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam
OUTGOING=${DIR1}/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam
#-----------------------------------------------------------------------------80
# Remove potential SMART adapter sequence (start site trimming).
#-----------------------------------------------------------------------------80
PolyATrimmer \
  INPUT=${INCOMING} \
  OUTPUT=${OUTGOING} \
  OUTPUT_SUMMARY=${OUTGOING}.report.txt \
  MISMATCHES=0 NUM_BASES=6 USE_NEW_TRIMMER=true \
  ADAPTER=

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

