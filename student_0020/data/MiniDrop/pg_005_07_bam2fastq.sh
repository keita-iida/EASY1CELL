#!/bin/sh

DIR1=005_bam
DIR2=006_fastq
INCOMING=${DIR1}/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam
OUTGOING=${DIR2}/pg_005_07_tiny1.cb.umi.qc.woAdapter.trimPolyA.fastq
#-----------------------------------------------------------------------------80
# Make directory.
#-----------------------------------------------------------------------------80
mkdir -p ${DIR2}

#-----------------------------------------------------------------------------80
# Change the data format from BAM to FASTQ.
#-----------------------------------------------------------------------------80
picard SamToFastq INPUT=${INCOMING} FASTQ=${OUTGOING}

