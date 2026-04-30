# DECODE
Gene-specific reads were identified by matching the first 10 nucleotides of each reference sequence (allowing one mismatch), and their relative abundance was calculated as the proportion of reads assigned to each gene among all aligned reads.
## Dependencies

- Python 3.6 or higher
- Standard library only (no pip install needed)

Built-in modules used:
- os, sys, argparse, re, typing, collections
## Directory Structure
workdir/
├── index/ # BWA reference index
│ └── ref.fa # Reference genome
├── 01.data/ # Input FASTQ files
│ ├── example.R1.fastq
│ └── example.R2.fastq
├── 02.bwa/ # Output directory (auto-created)
│ ├── {sample}.sam # SAM alignment file
│ ├── {sample}.filtered.bam # Filtered BAM file
│ ├── {sample}.fastq # Extracted FASTQ
│ ├── {sample}.extracted_seqs.txt # Extracted 19bp sequences
│ ├── {sample}.total_mapped.txt # Mapped read count
│ ├── {sample}.unknown.tsv # Unmatched sequences
│ ├── all_stats.tsv # Final gene statistics
│ └── summary.txt # Summary report
├── gene_NC.txt # NC group gene reference
├── gene_treat.txt # Treatment group gene reference
└── list # Sample list file
