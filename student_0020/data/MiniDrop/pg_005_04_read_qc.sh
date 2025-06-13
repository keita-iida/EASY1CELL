#!/bin/sh

DIR1=005_bam
INCOMING=${DIR1}/pg_005_03_tiny1.cb.umi.bam
OUTGOING=${DIR1}/pg_005_04_tiny1.cb.umi.qc.bam
#-----------------------------------------------------------------------------80
# Filtering low quality reads.
#-----------------------------------------------------------------------------80
FilterBam TAG_REJECT=XQ INPUT=${INCOMING} OUTPUT=${OUTGOING}

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

