#!/bin/sh

DIR1=006_fastq
REF=001_data/001_reference/STAR_index_chr8_MYC_masked
INCOMING=${DIR1}/pg_005_07_tiny1.cb.umi.qc.woAdapter.trimPolyA.fastq
#-----------------------------------------------------------------------------80
# Number of threads for STAR
#-----------------------------------------------------------------------------80
NCORES=1

#-----------------------------------------------------------------------------80
# STAR alignment: It is assumed that the Drop-seq cell barcode is assumed to be
# the first 12 bases of Read 1, and the UMI is assumed to be 13-20 bases.
#-----------------------------------------------------------------------------80
STAR --runThreadN ${NCORES} \
  --genomeDir ${REF} \
  --readFilesIn ${INCOMING} \
  --outFileNamePrefix ${DIR1}/pg_006_01_tiny1.star_ 

