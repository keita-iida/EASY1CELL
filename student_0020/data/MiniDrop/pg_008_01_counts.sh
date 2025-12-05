#!/bin/sh

DIR1=007_bam
DIR2=008_count
INCOMING=${DIR1}/pg_007_05_tiny1.clean.bam
OUTGOING=${DIR2}/pg_008_01_tiny1.txt.gz
#-----------------------------------------------------------------------------80
# Make directory.
#-----------------------------------------------------------------------------80
mkdir -p ${DIR2}

#-----------------------------------------------------------------------------80
# NUM_CORE_BARCODES is manually set.
#-----------------------------------------------------------------------------80

NUM_CORE_BARCODES=1

DigitalExpression \
  I=${INCOMING} \
  O=${OUTGOING} \
  SUMMARY=${OUTGOING}.summary.txt \
  NUM_CORE_BARCODES=${NUM_CORE_BARCODES}

