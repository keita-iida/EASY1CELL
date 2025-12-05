#=============================================================================80
#=============================================================================80
#
# Tutorial for Creating a Minimal scRNA-seq Dataset Featuring a Single Cell
# and the Gene MYC
#
# Author: Keita Iida and Kiwamu Arakane
#
#=============================================================================80
#=============================================================================80

This README provides a step-by-step tutorial for generating a minimal Drop-seq
single-cell RNA sequencing (scRNA-seq) dataset, intended for educational and
demonstrative purposes.
The example illustrates how to create a highly simplified dataset using real
paired-end Drop-seq data. The final output will contain only one cell barcode
and one gene - MYC (ENSG00000136997).

#=============================================================================80
# Environmental setting
#=============================================================================80
This tutorial requires a computational environment with the following tools
installed:

* samtools
* Picard
* Dropseq Tools
* STAR

It is recommended to use miniconda to manage the installation of these tools:

$ conda install samtools
$ conda install picard
$ conda install dropseq_tools
$ conda install star

You will also need to prepare the following input files:

* 002_input/cell_barcodes.txt - a list of desired cell barcodes
* 002_input/gene_ensemblids.txt - a list of gene IDs to retain (e.g., MYC)

Note: The last line in both text files must end with a newline (\n).
For example:

  AAA\\n
  BBB\\n
  CCC\\n
         <--- This final newline is required.

#=============================================================================80
# pg_004_01_create_fastq.sh
#=============================================================================80
This script manually creates small FASTQ files based on the sequences from
003_bam/tiny1.bam. Output files:

* 004_fastq/pg_004_01_tiny1_R1.fastq
* 004_fastq/pg_004_01_tiny1_R2.fastq

These files contain three reads associated with a single cell barcode.
Drop-seq assumes the following structure:

Read 1 (R1): Cell barcode (12 bp) + UMI (8 bp) + (polyT, spacer, linker, etc.)
Read 2 (R2): cDNA sequence

#-----------------------------------------------50
# Example: 004_fastq/pg_004_01_tiny1_R1.fastq
#-----------------------------------------------50
@TEST.1
AATACAACGGTAGGAGCAAAT
+
IIIIFFFFF0BBBB<F<<<<<
@TEST.2
AATACAACGGTACCACCAAAT
+
IIIIIIIIIIIIBB<F<<<<<
@TEST.3
AATACAACGGTAAAAAAAAAT
+
IIIIIIFFFIIIIIIIIFFF#
#-----------------------------------------------50

#-----------------------------------------------50
# Example: 004_fastq/pg_004_01_tiny1_R2.fastq
#-----------------------------------------------50
@TEST.1
GAATGTCAAGAGGCGAACACACAACGTCTTGGAGCGCCAGAGGAGGAACGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@TEST.2
TTTAATGTAACCTTGCTAAAGGAGTGATTTCTATTTCCTTTCTTAAAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@TEST.3
AAGCAGTGGTATCAACGCAGAGTGAATGGGTTTTTCCCCCGGGGGGTTTTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
#-----------------------------------------------50

Note: Read TEST.2 ends in AAAAAAA, which will later be interpreted as a polyA
tail.

#=============================================================================80
# pg_005_01_fastq2bam.sh
#=============================================================================80
This script uses Drop-seq tools to convert the FASTQ files (R1 and R2) into a
single BAM file:

* 005_bam/pg_005_01_tiny1.bam
* 005_bam/pg_005_01_tiny1.bam.sam

Header Explanation
#-----------------------------------------------50
@HD     VN:1.5  SO:queryname
@RG     ID:A    SM:tiny1
#-----------------------------------------------50
* @HD: Header line
  - VN:1.5: SAM version 1.5
  - SO:queryname: Sorted by query name

* @RG: Read group
  - ID:A: Group ID "A"
  - SM:tiny1: Sample name

Alignment Lines
#-----------------------------------------------50
TEST.1  77   *  0  0  *  *  0  0  AAT...    III...    RG:Z:
TEST.1  141  *  0  0  *  *  0  0  GAA...... III...... RG:Z:A
...
TEST.3  77   *  0  0  *  *  0  0  AAT...    III...    RG:Z:A
TEST.3  141  *  0  0  *  *  0  0  AAG...... III...... RG:Z:A
#-----------------------------------------------50
* 77: Binary flag indicating Read 1, paired, and unmapped
      77 = 64 (Read1) + 13 (unmapped) = Read1, paired, unmapped
* 141: Binary flag indicating Read 2, paired, and unmapped
       141 = 128（Read2）+ 13 = Read2, paired, unmapped
* RG:Z:A: Read group A.
* Use Picard's Explain SAM Flags to interpret these values:
  https://broadinstitute.github.io/picard/explain-flags.html

#=============================================================================80
# pg_005_02_cellbarcode_tag.sh
#=============================================================================80
This script adds cell barcode tags to the BAM file using Drop-seq tools:

* 005_bam/pg_005_02_tiny1.cb.bam
* 005_bam/pg_005_02_tiny1.cb.bam.sam

Example of Cell Barcode Tags
#-----------------------------------------------50
TEST.1  77   *  0  0  *  *  0  0  AAT...    III...                       RG:Z:A
TEST.1  141  *  0  0  *  *  0  0  GAA...... III...... XC:Z:AATACAACGGTA  RG:Z:A
...
TEST.3  77   *  0  0  *  *  0  0  AAT...    III...                       RG:Z:A
TEST.3  141  *  0  0  *  *  0  0  TTT...... III...... XC:Z:AATACAACGGTA  RG:Z:A
#-----------------------------------------------50
* XC:Z: indicates the cell barcode (12 bp). This is added only to Read 2 (R2).
* Drop-seq tools assume the tag name is XC. You can rename it, but it's not
  recommended, as some functions rely on this default.

Also note the addition of a new header line @PG, which records the program used
to modify the BAM file.

#=============================================================================80
# pg_005_03_umi_tag.sh
#=============================================================================80
This script tags the reads with UMI sequences (8 bp) using Drop-seq tools:

* 005_bam/pg_005_03_tiny1.cb.umi.bam
* 005_bam/pg_005_03_tiny1.cb.umi.bam.sam

Example with UMI Tags
#-----------------------------------------------50
TEST.1  4  *  0  0  *  *  0  0  GAA... III... XC:Z:AATACAACGGTA  RG:Z:A XM:Z:GGAGCAAA
TEST.2  4  *  0  0  *  *  0  0  TTT... III... XC:Z:AATACAACGGTA  RG:Z:A XM:Z:CCACCAAA
TEST.3  4  *  0  0  *  *  0  0  AAG... III... XC:Z:AATACAACGGTA  RG:Z:A XM:Z:AAAAAAAA
#-----------------------------------------------50
* XM:Z: is the tag used for UMI.
* Flag 4 indicates the read is unmapped (merged from Read 1 and Read 2).
* As with XC, changing the tag name is possible but not recommended.

#=============================================================================80
# pg_005_04_read_qc.sh
#=============================================================================80
This script removes low-quality reads tagged by XC using Drop-seq tools:

* 005_bam/pg_005_04_tiny1.cb.umi.qc.bam
* 005_bam/pg_005_04_tiny1.cb.umi.qc.bam.sam

In this example, the FASTQ files do not contain low-quality reads, so the output
is unchanged.

Practice Tip: Try simulating FASTQ files with low-quality reads to see how they
are filtered out in the resulting BAM file.

#=============================================================================80
# pg_005_05_remove_adapter.sh
#=============================================================================80
This script removes the SMART adapter sequence:
AAGCAGTGGTATCAACGCAGAGTGAATGGG
This sequence can appear at the 5′ end of cDNA reads and needs to be trimmed.

* 005_bam/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam
* 005_bam/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam.report.txt
* 005_bam/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam.sam

#-----------------------------------------------50
# 005_bam/pg_005_05_tiny1.cb.umi.qc.woAdapter.bam.sam
#-----------------------------------------------50
TEST.1  4  *  0  0  *  *  0  0  GAA... III... XC:Z:AAT... RG:Z:A XM:Z:GGA...
TEST.2  4  *  0  0  *  *  0  0  TTT... III... XC:Z:AAT... RG:Z:A XM:Z:CCA...
TEST.3  4  *  0  0  *  *  0  0  TTT... III... XC:Z:AAT... RG:Z:A XM:Z:AAA... ZS:i:30
#-----------------------------------------------50

For TEST.3, a new tag "ZS:i:30" appeared at the end, which means that the
first 30 bases were trimmed. One can see the following change:

From AAGCAGTGGTATCAACGCAGAGTGAATGGGTTTTTCCCCCGGGGGGTTTTT
to                                 TTTTTCCCCCGGGGGGTTTTT

#=============================================================================80
# pg_005_06_trim_polya.sh
#=============================================================================80
This script trims polyA tails from the 3' end of reads. The tail is detected
when NUM_BASES or more consecutive A bases are found (default: 6).

Output files:
* 005_bam/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam
* 005_bam/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam.report.txt
* 005_bam/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam.sam

Example:
#-----------------------------------------------50
# 005_bam/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam.sam
#-----------------------------------------------50
TEST.1  4  *  0  0  *  *  0  0  GAA... III... XC:Z:AAT... RG:Z:A XM:Z:GGA...
TEST.2  4  *  0  0  *  *  0  0  TTT... III... XC:Z:AAT... RG:Z:A XM:Z:CCA... ZP:i:45
TEST.3  4  *  0  0  *  *  0  0  TTT... III... XC:Z:AAT... RG:Z:A XM:Z:AAA... ZS:i:30
#-----------------------------------------------50

For TEST.2, a new tag "ZP:i:45" appeared at the end, which means that polyA
was detected from 45th base, which were trimmed here:

Before  CAAGAGGCGAACACACAACGTCTTGGAGCGCCAGAGGAGGAACGAAAAAAA
After   CAAGAGGCGAACACACAACGTCTTGGAGCGCCAGAGGAGGAACG

Reads like TEST.3 with internal A-rich regions are not trimmed:

#=============================================================================80
# pg_005_07_bam2fastq.sh
#=============================================================================80
This script converts BAM to FASTQ using Drop-seq tools:

* 006_fastq/pg_005_07_tiny1.cb.umi.qc.woAdapter.trimPolyA.fastq

#-----------------------------------------------50
# 006_fastq/pg_005_07_tiny1.cb.umi.qc.woAdapter.trimPolyA.fastq
#-----------------------------------------------50
@TEST.1
GAATGTCAAGAGGCGAACACACAACGTCTTGGAGCGCCAGAGGAGGAACGA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@TEST.2
TTTAATGTAACCTTGCTAAAGGAGTGATTTCTATTTCCTTTCTT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@TEST.3
TTTTTCCCCCGGGGGGTTTTT
+
IIIIIIIIIIIIIIIIIIIII
#-----------------------------------------------50

#=============================================================================80
# pg_006_01_star__________wREF.sh
#=============================================================================80
This script runs STAR to align reads to the reference genome:

* 006_fastq/pg_006_01_tiny1.star_Aligned.out.sam
* 006_fastq/pg_006_01_tiny1.star_Log.final.out
* 006_fastq/pg_006_01_tiny1.star_Log.out
* 006_fastq/pg_006_01_tiny1.star_Log.progress.out
* 006_fastq/pg_006_01_tiny1.star_SJ.out.tab

Note: Cell barcode (XC) and UMI (XM) tags are not included in STAR's SAM output.

Header Explanation (006_fastq/pg_006_01_tiny1.star_Aligned.out.sam)
#-----------------------------------------------50
* @SQ: Sequence Dictionary, describing information about each chromosome and
       contig in the reference genome.
  - SN:8: Sequence Name. Name of the reference sequence (e.g., chromosome name).
  - LN:145138636: Length. The length (in bases) of the corresponding reference
                  sequence (in this case chromosome 8), which is 145,138,636 bp
                  in the above example.
#-----------------------------------------------50

Alignment Lines (006_fastq/pg_006_01_tiny1.star_Aligned.out.sam)
#-----------------------------------------------50
TEST.1  0    8  127740694  255  51M    *  0  0  GAATGTCAAGAGGCGAACACA... III... NH:i:1 HI:i:1 AS:i:50 nM:i:0
TEST.2  0    8  127740348  255  44M    *  0  0  TTTAATGTAACCTTGCTAAAG... III... NH:i:1 HI:i:1 AS:i:43 nM:i:0
TEST.3  256  8  139003813  0    13M8S  *  0  0  TTTTTCCCCCGGGGGGTTTTT    III... NH:i:7 HI:i:1 AS:i:12 nM:i:0
TEST.3  256  8  38719248   0    13M8S  *  0  0  TTTTTCCCCCGGGGGGTTTTT    III... NH:i:7 HI:i:2 AS:i:12 nM:i:0
TEST.3  16   8  70608363   0    5S16M  *  0  0  AAAAACCCCCCGGGGGAAAAA    III... NH:i:7 HI:i:3 AS:i:13 nM:i:1
TEST.3  256  8  132058839  0    13M8S  *  0  0  TTTTTCCCCCGGGGGGTTTTT    III... NH:i:7 HI:i:4 AS:i:12 nM:i:0
TEST.3  272  8  72390115   0    13M8S  *  0  0  AAAAACCCCCCGGGGGAAAAA    III... NH:i:7 HI:i:5 AS:i:12 nM:i:0
TEST.3  272  8  80279758   0    20M1S  *  0  0  AAAAACCCCCCGGGGGAAAAA    III... NH:i:7 HI:i:6 AS:i:13 nM:i:3
TEST.3  256  8  7714463    0    6S15M  *  0  0  TTTTTCCCCCGGGGGGTTTTT    III... NH:i:7 HI:i:7 AS:i:12 nM:i:1
#-----------------------------------------------50
Interestingly, STAR automatically duplicated TEST.3 to test several candidates.

For BAM flags, reference the following:
  https://broadinstitute.github.io/picard/explain-flags.html

@CO,                        This is a comment header describing the used command.
TEST.1  0,                  Single-end reads that are uniquely aligned in the
                            sense direction. Although the original data is
                            paired-end, the created BAM file appears single-end.
TEST.1  ...  8  127740694,  This read was aligned onto this position in
                            chromosome 8.
TEST.1  ...  255,           This alignment is unique, but sometimes means
                            "invalid."
TEST.1  ...  51M,           CIGAR string, describing how many aligned reads have
                            mappings, insertions, deletions, etc.
TEST.1  ...  *  0  0,       RNEXT, PNEXT, TLEN, where
                            RNEXT: read name of the partner, "*,"
                            PNEXT: read position of the partner, "0,"
                            TLEN: total length of the template, "0."
NH:i:1,                     The total number of locations this read was mapped to.
HI:i:1,                     Hit index (to distinguish between multiple hits).
AS:i:50,                    Alignment score (higher is better). It equals
                              −10 log10 Pr{mapping position is wrong},
                            rounded to the nearest integer.
nM:i:0,                     Number of nucleotide mismatches.

#=============================================================================80
# pg_007_01_sorting.sh
#=============================================================================80
This script sorts BAM by query name using Picard's SortSam, preparing for later
merging steps:

* 007_bam/pg_007_01_tiny1.aligned.sorted.bam
* 007_bam/pg_007_01_tiny1.aligned.sorted.bam.sam

Note: Sorting only affects tag order; sequence content remains unchanged.

#=============================================================================80
# pg_007_02_cb_umi__________wREF.sh
#=============================================================================80
This script merges STAR-aligned BAM (with no cell barcode/UMI) with the
barcode-tagged original BAM using Picard’s MergeBamAlignment.

Input:
* Aligned BAM: 007_bam/pg_007_01_tiny1.aligned.sorted.bam
* Tagged BAM: 005_bam/pg_005_06_tiny1.cb.umi.qc.woAdapter.trimPolyA.bam

Output:
* 007_bam/pg_007_02_tiny1.merged.bam
* 007_bam/pg_007_02_tiny1.merged.bam.sam

#-----------------------------------------------50
# 007_bam/pg_007_02_tiny1.merged.bam.sam
#-----------------------------------------------50
@HD     VN:1.5  SO:coordinate
@SQ     SN:8    LN:145138636    M5:c67955b5f...
@RG     ID:A    SM:tiny1
@PG     ID:STAR PN:STAR VN:2.7.11b      CL:...
TEST.2 ... 44M   ... TTT...... III... XC:Z:AAT... MD:Z:44   ... XM:Z:CCA... ZP:i:45 UQ:i:0  AS:i:43
TEST.1 ... 51M   ... GAA...... III... XC:Z:AAT... MD:Z:51   ... XM:Z:GGA... UQ:i:0  AS:i:50
TEST.3 ...           TTT...    III... XC:Z:AAT... PG:Z:STAR ... XM:Z:AAA... ZS:i:30
#-----------------------------------------------50

Explanation of fields:
* MD:Z:2G13
  Mismatch descriptor: match 2 bases, then mismatch (ref has G), then match 13 bases.
* UQ:i:0
  Uncertainty Quality: the sum of mismatch quality scores. A higher value
  indicates lower alignment confidence.
* ZS:i:30
  Supplementary alignment score from STAR. Occasionally used for filtering or
  scoring.
* ZP:i:45
  A value estimating the uniqueness and reliability of the alignment.
  Higher values indicate greater mapping confidence.

#=============================================================================80
# pg_007_03_genes__________wREF.sh
#=============================================================================80
This script assigns gene-level annotations to each aligned read, tagging them
with gene names (GE), IDs (GN), and transcript regions (XF: exonic/intronic/intergenic):

* 007_bam/pg_007_03_tiny1.star_gene_exon_tagged.bam
* 007_bam/pg_007_03_tiny1.star_gene_exon_tagged.bam.sam

#-----------------------------------------------50
# 007_bam/pg_007_03_tiny1.star_gene_exon_tagged.bam.sam
#-----------------------------------------------50
TEST.2 ... XC:Z:AAT... XF:Z:INTRONIC ... XM:Z:CCA... gf:Z:INTRONIC gn:Z:ENSG00000136997 gs:Z:+
TEST.1 ... XC:Z:AAT... XF:Z:CODING   ... XM:Z:GGA... gf:Z:CODING   gn:Z:ENSG00000136997 gs:Z:+
TEST.3 ... XC:Z:AAT...                   XM:Z:AAA... ZS:i:30
#-----------------------------------------------50
Explanation of tags:
* XF:Z:INTRONIC
  Region classification tag: indicates whether the read maps to a
  coding exon (CODING), intron (INTRONIC), untranslated region (UTR) or
  intergenic region.
* gf:Z:INTRONIC
  Gene function tag: same information as XF, added redundantly for tool
  compatibility.
* gn:Z:ENSG00000175806
  Gene ID (Ensembl format) corresponding to the mapped gene.
* gs:Z:+
  Strand information: indicates the gene's transcriptional strand (+ or -),
  useful for interpreting read direction and restoring transcriptional
  orientation.

#=============================================================================80
# pg_007_04_fix_beads.sh
#=============================================================================80
This script corrects point mutations in cell barcodes—typically single
nucleotide substitutions such as GTTAC → GTTTC.

Notes
* The output file will mostly look the same, except for an added @PG ID:1 tag
  in the header.
* If multiple similar barcodes share UMIs, they may be merged and classified
  as a barcode "subtype."

Practice Tip: Try simulating FASTQ files with point-mutated cell barcodes to
see how they are fixed in the resulting BAM file.

#=============================================================================80
# pg_007_05_fix2nd_beads.sh
#=============================================================================80
This script corrects synthesis errors in barcodes, particularly insertions and
deletions (indels), targeting experimental errors that produce multiple
near-identical barcodes.

#=============================================================================80
# pg_008_01_counts.sh
#=============================================================================80
This script counts the number of reads mapped to each gene for each cell.

* 008_count/pg_008_01_tiny1.txt.gz
* 008_count/pg_008_01_tiny1.txt.gz.summary.txt

#-----------------------------------------------50
# 008_count/pg_008_01_tiny1.txt.gz
#-----------------------------------------------50
GENE    AATACAACGGTA
ENSG00000136997 1

Notes
* Each entry represents the number of reads per gene (row) per cell (column).

