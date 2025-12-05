#!/bin/sh

OUTDIR="004_fastq"
R1=${OUTDIR}/pg_004_01_tiny1_R1.fastq
R2=${OUTDIR}/pg_004_01_tiny1_R2.fastq
#-----------------------------------------------------------------------------80
# Make directory.
#-----------------------------------------------------------------------------80
mkdir -p ${OUTDIR}
rm -f ${OUTDIR}/pg_004_01_tiny1_*
#---------------------------------------------------------60
# First Read (Normal read)
#---------------------------------------------------------60
# R1
CELLBC1="AATACAACGGTA"               # 12 bases
UMI1="GGAGCAAA"                      # 8 bases
POLYT1="T"                           # 1 bases
SEQ1="${CELLBC1}${UMI1}${POLYT1}"
QUAL1="IIIIFFFFF0BBBB<F<<<<<"        # III...II doesn't work for Drop-seq tools.

# R2
SEQ2="GAATGTCAAGAGGCGAACACACAACGTCTTGGAGCGCCAGAGGAGGAACGA"     # R2 (transcript)
QUAL2=$(printf 'I%.0s' $(seq 1 ${#SEQ2}))
#---------------------------------------------------------60
# Second Read (Normal read with poly A tail)
#---------------------------------------------------------60
# R1
CELLBC3="AATACAACGGTA"               # 12 bases
UMI3="CCACCAAA"                      # 8 bases
POLYT3="T"                           # 1 bases
SEQ3="${CELLBC3}${UMI3}${POLYT3}"
QUAL3="IIIIIIIIIIIIBB<F<<<<<"        # III...II doesn't work for Drop-seq tools.

# R2
POLYA4="AAAAAAA"
SEQ4="TTTAATGTAACCTTGCTAAAGGAGTGATTTCTATTTCCTTTCTT${POLYA4}"
QUAL4=$(printf 'I%.0s' $(seq 1 ${#SEQ4}))
#---------------------------------------------------------60
# Third Read (Artifact including SMART adapter sequence)
#---------------------------------------------------------60
# R1
CELLBC5="AATACAACGGTA"               # 12 bases
UMI5="AAAAAAAA"                      # 8 bases
POLYT5="T"                           # 1 bases
SEQ5="${CELLBC5}${UMI5}${POLYT5}"
QUAL5="IIIIIIFFFIIIIIIIIFFF#"        # III...II doesn't work for Drop-seq tools.

# R2
ADAPTER6="AAGCAGTGGTATCAACGCAGAGTGAATGGG"
SEQ6="${ADAPTER6}TTTTTCCCCCAAAAAATTTTT"              # Dummy transcript sequence
QUAL6=$(printf 'I%.0s' $(seq 1 ${#SEQ6}))
#-----------------------------------------------------------------------------80
# Output
#-----------------------------------------------------------------------------80
echo "@TEST.1" > $R1
echo "$SEQ1" >> $R1
echo "+" >> $R1
echo "$QUAL1" >> $R1

echo "@TEST.1" > $R2
echo "$SEQ2" >> $R2
echo "+" >> $R2
echo "$QUAL2" >> $R2

echo "@TEST.2" >> $R1
echo "$SEQ3" >> $R1
echo "+" >> $R1
echo "$QUAL3" >> $R1

echo "@TEST.2" >> $R2
echo "$SEQ4" >> $R2
echo "+" >> $R2
echo "$QUAL4" >> $R2

echo "@TEST.3" >> $R1
echo "$SEQ5" >> $R1
echo "+" >> $R1
echo "$QUAL5" >> $R1

echo "@TEST.3" >> $R2
echo "$SEQ6" >> $R2
echo "+" >> $R2
echo "$QUAL6" >> $R2

echo ""
echo "FASTQ files written to ${R1} and ${R2}"
echo ""
