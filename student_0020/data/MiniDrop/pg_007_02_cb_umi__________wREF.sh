#!/bin/sh

DIR1=005_bam
DIR2=007_bam
REF=001_data/001_reference/chr8.fa
UNMAPPED=${DIR1}/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam
INCOMING=${DIR2}/pg_007_01_tiny1.aligned.sorted.bam
OUTGOING=${DIR2}/pg_007_02_tiny1.merged.bam
#-----------------------------------------------------------------------------80
# This program requires 001_data/001_reference/MYC_region.dict.
#-----------------------------------------------------------------------------80
picard MergeBamAlignment ALIGNED_BAM=${INCOMING} \
  UNMAPPED_BAM=${UNMAPPED} \
  OUTPUT=${OUTGOING} \
  REFERENCE_SEQUENCE=${REF} \
  INCLUDE_SECONDARY_ALIGNMENTS=false PAIRED_RUN=false

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

