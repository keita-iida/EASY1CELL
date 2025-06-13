#!/bin/sh

DIR1=006_fastq
DIR2=007_bam
INCOMING=${DIR1}/pg_006_01_tiny1.star_Aligned.out.sam
OUTGOING=${DIR2}/pg_007_01_tiny1.aligned.sorted.bam
#-----------------------------------------------------------------------------80
# Make directory.
#-----------------------------------------------------------------------------80
mkdir -p ${DIR2}

#-----------------------------------------------------------------------------80
# Sort the alignment result.
#-----------------------------------------------------------------------------80
picard SortSam I=${INCOMING} O=${OUTGOING} SO=queryname

samtools view -h ${OUTGOING} > ${OUTGOING}.sam
