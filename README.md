# DECODE
Gene-specific reads were identified by matching the first 10 nucleotides of each reference sequence (allowing one mismatch), and their relative abundance was calculated as the proportion of reads assigned to each gene among all aligned reads.
### Software Dependencies

| Software | Version Tested | Minimum Version |
|----------|---------------|-----------------|
| Python | 3.6, 3.7, 3.8, 3.9, 3.10 | 3.6 |
| BWA | 0.7.17 | 0.7.12 |
| SAMtools | 1.9, 1.10, 1.11 | 1.9 |

### Python Dependencies (Standard Library Only)
No external pip packages required. Uses built-in modules:
- `os`, `sys`, `argparse`, `re`, `typing`, `collections`

### Hardware Requirements
- **No non-standard hardware required**

### Tested On
- Ubuntu 20.04 LTS (Jammy Jellyfish)


### Typical Install Time

**2-5 minutes** on a normal desktop computer (Intel Core i5 or equivalent, 8GB RAM, SSD storage, 100 Mbps internet connection).

**Detailed breakdown:**

| Step | Time |
|------|------|
| Install BWA, SAMtools, GNU Parallel | 60-90 seconds |
| Download DECODE scripts | 5-10 seconds |
| Set executable permissions | < 5 seconds |
| Build BWA index (demo reference) | 30-60 seconds |
| **Total** | **~2-5 minutes** |

**For users with dependencies pre-installed:** Less than 1 minute.

**Slow internet connection (< 10 Mbps):** Allow 10-15 minutes for dependency downloads.

### Run with Example Data

# Create demo directory
-mkdir -p demo_workdir/{01.data,02.bwa,index}
-cd demo_workdir

# Create reference genome
index/ref.fa 

# Build BWA index
bwa index index/ref.fa

# Create gene reference files
-gene_NC.txt 
-gene_treat.txt 

# Create sample list
cat > list << 'EOF'
example
EOF

# Create demo FASTQ files
-01.data/example.R1.fastq 
-01.data/example.R2.fastq 

# Run gene extraction
../gene_extract.sh

# Run gene counting
python3 ../gene_count.py --workdir . --mismatch 1 --include-mismatch

## Directory Structure

**workdir/**
- **index/** - BWA reference index
  - `ref.fa` - Reference genome
- **01.data/** - Input FASTQ files
  - `example.R1.fastq`
  - `example.R2.fastq`
- **02.bwa/** - Output directory (auto-created)
  - `{sample}.sam` - SAM alignment file
  - `{sample}.filtered.bam` - Filtered BAM file
  - `{sample}.fastq` - Extracted FASTQ
  - `{sample}.extracted_seqs.txt` - Extracted 19bp sequences
  - `{sample}.total_mapped.txt` - Mapped read count
  - `{sample}.unknown.tsv` - Unmatched sequences
  - `all_stats.tsv` - Final gene statistics
  - `summary.txt` - Summary report
- `gene_NC.txt` - NC group gene reference
- `gene_treat.txt` - Treatment group gene reference
- `list` - Sample list file
