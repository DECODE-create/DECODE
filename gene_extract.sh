#!/bin/bash
Workdir=""

process_sample() {
    sample=$1
    
    echo "=== Processing $sample ==="
    
    # Step 1: align and save raw sam
    echo "1. Aligning..."
    bwa mem -t 4 -M $Workdir/index/ref.fa \
        $Workdir/${sample}.R1.fastq.gz \
        $Workdir/${sample}.R2.fastq.gz \
        > $Workdir/bwa/${sample}.sam 2> $Workdir/bwa/${sample}.align.log
    
    # Step 2: convert to bam and drop unmapped reads
    echo "2. Converting to BAM..."
    samtools view -bS -F 4 $Workdir/02.bwa/${sample}.sam \
        > $Workdir/bwa/${sample}.filtered.bam 2> $Workdir/bwa/${sample}.view.log
    
    # Step 3: count mapped reads
    total=$(samtools view -c $Workdir/02.bwa/${sample}.filtered.bam)
    echo "3. Total mapped reads: $total"
    echo "$total" > $Workdir/bwa/${sample}.total_mapped.txt
    
    if [ "$total" -eq 0 ]; then
        echo "WARNING: No mapped reads for $sample"
        return 1
    fi
    
    # Step 4: dump back to fastq
    echo "4. Extracting FASTQ..."
    samtools fastq $Workdir/bwa/${sample}.filtered.bam \
        > $Workdir/bwa/${sample}.fastq 2> $Workdir/bwa/${sample}.fastq.log
    
    # Step 5: pull out target sequences (the key step you wanted)
    echo "5. Extracting target sequences..."
    grep -o 'AAACACCG...................' $Workdir/02.bwa/${sample}.fastq | \
    awk '{print substr($0, 11, 19)}' \
    > $Workdir/02.bwa/${sample}.extracted_seqs.txt
    
    # Show how many we got and a peek
    seq_count=$(wc -l < $Workdir/02.bwa/${sample}.extracted_seqs.txt)
    echo "6. Extracted $seq_count sequences"
    
    if [ "$seq_count" -eq 0 ]; then
        echo "WARNING: No target sequences extracted for $sample"
        echo "Checking first few lines of FASTQ:"
        head -4 $Workdir/02.bwa/${sample}.fastq
    else
        echo "First 5 extracted sequences:"
        head -5 $Workdir/02.bwa/${sample}.extracted_seqs.txt
    fi
    
    echo "Done with $sample"
    echo ""
}

export -f process_sample
export Workdir

cat list | parallel -j 4 process_sample
