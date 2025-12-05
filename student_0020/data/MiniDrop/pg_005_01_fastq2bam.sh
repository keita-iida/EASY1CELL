#!/bin/sh

DIR1=004_fastq
DIR2=005_bam
OUTGOING=${DIR2}/pg_005_01_tiny1.bam
#-----------------------------------------------------------------------------80
# Change the data format from FASTQ to BAM for using Drop-seq tools.
#-----------------------------------------------------------------------------80
mkdir -p ${DIR2}

FILE1=${DIR1}/pg_004_01_tiny1_R1.fastq
FILE2=${DIR1}/pg_004_01_tiny1_R2.fastq
picard FastqToSam F1=${FILE1} F2=${FILE2} O=${OUTGOING} SM=tiny1

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

