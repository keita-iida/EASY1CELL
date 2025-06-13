#!/bin/sh

DIR1=007_bam
INCOMING=${DIR1}/pg_007_03_tiny1.star_gene_exon_tagged.bam
OUTGOING=${DIR1}/pg_007_04_tiny1.clean_substitution.bam
#-----------------------------------------------------------------------------80
# Make directory.
#-----------------------------------------------------------------------------80
mkdir -p ${DIR1}/tmp

#-----------------------------------------------------------------------------80
# Fix bead errors (substitution).
# Caution!!
# DetectBeadSubstitutionErrors automatically changes the tags of
# corrected cell barcodes and UMIs into XC and XM, respectively.
#-----------------------------------------------------------------------------80
DetectBeadSubstitutionErrors \
  I=${INCOMING} \
  O=${OUTGOING} \
  OUTPUT_REPORT=${OUTGOING}.report.txt \
  TMP_DIR=${DIR1}/tmp \
  MIN_UMIS_PER_CELL=1

samtools view -h ${OUTGOING} > ${OUTGOING}.sam

