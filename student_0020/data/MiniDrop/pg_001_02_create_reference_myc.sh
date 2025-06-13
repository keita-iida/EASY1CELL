#!/bin/sh

DIR1=001_data/001_reference
INCOMING=${DIR1}/myc.fa
OUTGOING=${DIR1}/myc.dict
#-----------------------------------------------------------------------------80
# Sort the alignment result.
#-----------------------------------------------------------------------------80
STAR --runThreadN 1 \
  --runMode genomeGenerate \
  --genomeDir ${DIR1}/STAR_index_myc \
  --genomeFastaFiles ${DIR1}/myc.fa \
  --sjdbGTFfile ${DIR1}/myc.gtf \
  --sjdbOverhang 50 \
  --genomeSAindexNbases 10

picard CreateSequenceDictionary R=${INCOMING} O=${OUTGOING}

