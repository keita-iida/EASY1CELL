#!/bin/sh

DIR1=005_bam
INCOMING=${DIR1}/pg_005_02_tiny1.cb.bam
OUTGOING=${DIR1}/pg_005_03_tiny1.cb.umi.bam
#-----------------------------------------------------------------------------80
# Extract UMIs and tag them.
# BASE_RANGE information is written in the original paper.
#-----------------------------------------------------------------------------80
TagBamWithReadSequenceExtended \
  INPUT=${INCOMING} \
  OUTPUT=${OUTGOING} \
  SUMMARY=${OUTGOING}.summary.txt \
  BASE_RANGE=13-20 BASE_QUALITY=10 BARCODED_READ=1 DISCARD_READ=TRUE \
  TAG_NAME=XM NUM_BASES_BELOW_QUALITY=1

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

