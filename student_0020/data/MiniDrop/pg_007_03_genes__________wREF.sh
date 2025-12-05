#!/bin/sh

DIR1=005_bam
DIR2=007_bam
REF=001_data/001_reference/genes.refFlat
UNMAPPED=${DIR1}/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam
INCOMING=${DIR2}/pg_007_02_tiny1.merged.bam
OUTGOING=${DIR2}/pg_007_03_tiny1.star_gene_exon_tagged.bam
#-----------------------------------------------------------------------------80
# Add gene information.
#-----------------------------------------------------------------------------80
TagReadWithGeneFunction I=${INCOMING} O=${OUTGOING} ANNOTATIONS_FILE=${REF}

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

