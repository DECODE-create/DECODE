#!/usr/bin/env python3
"""
Gene sequence statistics script (supports mismatch matching)
Function: Count gene frequencies from extracted sequence files, supports mismatch matching
Usage: python gene_count.py --workdir /path/to/workdir --mismatch 1
workdir/
├── 02.bwa/                    
│   ├── {sample}.total_mapped.txt    
│   ├── {sample}.extracted_seqs.txt  
│   └── {sample}.unknown.tsv         
├── gene_NC.txt                
├── gene_treat.txt            
└── list                      
"""

import os
import sys
import argparse
import re
from typing import Dict, List, Tuple
from collections import defaultdict, Counter

def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate Hamming distance (sequences must be same length)"""
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequence length mismatch: {seq1} ({len(seq1)}) vs {seq2} ({len(seq2)})")
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def load_gene_dict(gene_list_file: str) -> Tuple[Dict[str, str], List[str]]:
    """Load gene dictionary, returns {sequence: gene_name} and list of gene names"""
    gene_dict = {}
    gene_names = set()
    
    if not os.path.exists(gene_list_file):
        raise FileNotFoundError(f"Gene list file does not exist: {gene_list_file}")
    
    with open(gene_list_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                print(f"Warning: Line {line_num} format error: {line}")
                continue
                
            gene_name = parts[0]
            gene_seq = parts[1].upper()
            
            # Validate sequence length (should be 19bp)
            if len(gene_seq) != 19:
                print(f"Warning: Gene {gene_name} sequence length is not 19bp: {gene_seq}")
                continue
                
            if gene_seq in gene_dict:
                print(f"Warning: Sequence {gene_seq} duplicated: {gene_dict[gene_seq]} and {gene_name}")
            
            gene_dict[gene_seq] = gene_name
            gene_names.add(gene_name)
    
    gene_names_sorted = sorted(gene_names)
    print(f"Loaded {len(gene_dict)} gene sequences, {len(gene_names_sorted)} gene names")
    
    return gene_dict, gene_names_sorted

def match_with_mismatch(query_seq: str, gene_dict: Dict[str, str], 
                       max_mismatch: int = 0) -> Tuple[str, int]:
    """
    Gene matching with mismatch support
    Returns: (matched gene name, mismatch count) or (None, minimum mismatch count)
    """
    query_seq = query_seq.upper()
    
    # Exact match
    if query_seq in gene_dict:
        return gene_dict[query_seq], 0
    
    # If mismatch matching is not required
    if max_mismatch == 0:
        return None, 0
    
    # Fuzzy matching: find gene with minimum mismatches
    best_match = None
    min_mismatch = float('inf')
    
    for ref_seq, gene_name in gene_dict.items():
        mismatch = hamming_distance(query_seq, ref_seq)
        if mismatch < min_mismatch:
            min_mismatch = mismatch
            best_match = gene_name
            if min_mismatch == 0:  # Found perfect match (should not happen as checked above)
                break
    
    # Check if within allowed mismatch range
    if min_mismatch <= max_mismatch:
        return best_match, min_mismatch
    else:
        return None, min_mismatch

def get_gene_list_file(workdir: str, sample_name: str) -> str:
    """
    Select gene list file based on sample name
    Rule: If sample name starts with 'NC', use gene_NC.txt, otherwise use gene_treat.txt
    """
    if sample_name.startswith('NC'):
        return os.path.join(workdir, "gene_NC.txt")
    else:
        return os.path.join(workdir, "gene_treat.txt")

def process_sample(sample: str, workdir: str, max_mismatch: int = 0) -> List[Tuple]:
    """
    Process a single sample
    Returns: [(sample, gene_name, count, total, ratio, mismatch, gene_file_used), ...]
    """
    # Determine which gene list file to use based on sample name
    gene_list_file = get_gene_list_file(workdir, sample)
    
    if not os.path.exists(gene_list_file):
        print(f"  Error: Gene list file does not exist: {gene_list_file}")
        return []
    
    # Load gene dictionary
    try:
        gene_dict, gene_names = load_gene_dict(gene_list_file)
    except Exception as e:
        print(f"  Failed to load gene dictionary: {e}")
        return []
    
    print(f"  Sample {sample}: Using gene list file {os.path.basename(gene_list_file)}")
    
    # File paths
    total_file = os.path.join(workdir, "02.bwa", f"{sample}.total_mapped.txt")
    seq_file = os.path.join(workdir, "02.bwa", f"{sample}.extracted_seqs.txt")
    output_dir = os.path.join(workdir, "02.bwa")
    
    # Read total mapped reads count
    try:
        with open(total_file, 'r') as f:
            total = int(f.read().strip())
    except (FileNotFoundError, ValueError):
        print(f"  Warning: Cannot read total mapped reads file {total_file}")
        total = 0
    
    # Initialize statistics
    gene_counts = Counter()
    mismatch_info = {}  # Record maximum mismatch for each gene
    unknown_counts = Counter()
    
    # Read and count extracted sequences
    if os.path.exists(seq_file) and total > 0:
        seq_count = 0
        with open(seq_file, 'r') as f:
            for line in f:
                seq = line.strip().upper()
                if not seq or len(seq) != 19:
                    continue
                    
                seq_count += 1
                gene_name, mismatch = match_with_mismatch(seq, gene_dict, max_mismatch)
                
                if gene_name:
                    gene_counts[gene_name] += 1
                    # Record maximum mismatch
                    if gene_name not in mismatch_info or mismatch > mismatch_info[gene_name]:
                        mismatch_info[gene_name] = mismatch
                else:
                    unknown_counts[seq] += 1
        
        print(f"  Total mapped reads={total:,}, Extracted sequences={seq_count:,}")
        
        # Save unknown sequences to file
        if unknown_counts:
            unknown_file = os.path.join(output_dir, f"{sample}.unknown.tsv")
            with open(unknown_file, 'w') as f:
                f.write("sequence\tcount\tratio\n")
                for seq, count in unknown_counts.most_common():
                    ratio = count / total
                    f.write(f"{seq}\t{count}\t{ratio:.6f}\n")
            print(f"  Unknown sequences saved to: {unknown_file} ({len(unknown_counts)} unique sequences)")
    else:
        if not os.path.exists(seq_file):
            print(f"  Warning: Sequence file does not exist {seq_file}")
        elif total == 0:
            print(f"  Warning: Total mapped reads count is 0")
    
    # Generate results list
    results = []
    for gene_name in gene_names:
        count = gene_counts.get(gene_name, 0)
        max_mismatch_for_gene = mismatch_info.get(gene_name, 0)
        
        if total > 0:
            ratio = count / total
        else:
            ratio = 0.0
        
        results.append((
            sample,
            gene_name,
            count,
            total,
            f"{ratio:.6f}",
            max_mismatch_for_gene,
            os.path.basename(gene_list_file)  # Record which gene list file was used
        ))
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Gene sequence statistics script (supports mismatch matching)')
    parser.add_argument('--workdir', required=True, help='Working directory path')
    parser.add_argument('--mismatch', type=int, default=0, 
                       help='Maximum allowed mismatches (default: 0, exact match)')
    parser.add_argument('--list-file', default='list', 
                       help='Sample list file (default: list)')
    parser.add_argument('--output', default='all_stats.tsv',
                       help='Output file name (default: all_stats.tsv)')
    parser.add_argument('--include-mismatch', action='store_true',
                       help='Include mismatch column in output')
    parser.add_argument('--include-gene-file', action='store_true',
                       help='Include gene list file column in output')
    
    args = parser.parse_args()
    
    # Check files
    workdir = args.workdir
    list_file = os.path.join(workdir, args.list_file)
    output_dir = os.path.join(workdir, "02.bwa")
    
    if not os.path.exists(list_file):
        print(f"Error: Sample list file does not exist: {list_file}")
        sys.exit(1)
    
    if not os.path.exists(output_dir):
        print(f"Error: Output directory does not exist: {output_dir}")
        print("Please run the sequence extraction script first")
        sys.exit(1)
    
    # Check if gene list files exist
    gene_nc_file = os.path.join(workdir, "gene_NC.txt")
    gene_treat_file = os.path.join(workdir, "gene_treat.txt")
    
    if not os.path.exists(gene_nc_file):
        print(f"Warning: NC group gene list file does not exist: {gene_nc_file}")
    if not os.path.exists(gene_treat_file):
        print(f"Warning: Treatment group gene list file does not exist: {gene_treat_file}")
    
    if not os.path.exists(gene_nc_file) and not os.path.exists(gene_treat_file):
        print(f"Error: No gene list files found")
        print(f"Please ensure at least one file exists: {gene_nc_file} or {gene_treat_file}")
        sys.exit(1)
    
    # Load sample list
    with open(list_file, 'r') as f:
        samples = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(samples)} samples")
    print(f"Maximum allowed mismatches: {args.mismatch}")
    print(f"NC group gene list: {os.path.basename(gene_nc_file) if os.path.exists(gene_nc_file) else 'Not exists'}")
    print(f"Treatment group gene list: {os.path.basename(gene_treat_file) if os.path.exists(gene_treat_file) else 'Not exists'}")
    
    # Count sample types
    nc_samples = [s for s in samples if s.startswith('NC')]
    treat_samples = [s for s in samples if not s.startswith('NC')]
    print(f"NC group samples: {len(nc_samples)}")
    print(f"Treatment group samples: {len(treat_samples)}")
    
    # Process all samples
    all_results = []
    unknown_total = 0
    gene_files_used = set()
    
    for sample in samples:
        print(f"\nProcessing sample: {sample}")
        
        results = process_sample(sample, workdir, args.mismatch)
        
        if results:
            all_results.extend(results)
            gene_files_used.add(results[0][6])  # Record which gene file was used
            
            # Count total unknown sequences
            unknown_file = os.path.join(output_dir, f"{sample}.unknown.tsv")
            if os.path.exists(unknown_file):
                with open(unknown_file, 'r') as f:
                    lines = f.readlines()
                    unknown_total += len(lines) - 1  # Subtract header line
    
    # Write output file
    output_file = os.path.join(output_dir, args.output)
    
    # Build header
    header_parts = ["sample", "gene_name", "gene_count", "total_mapped", "ratio"]
    if args.include_mismatch:
        header_parts.append("max_mismatch")
    if args.include_gene_file:
        header_parts.append("gene_file")
    header = "\t".join(header_parts) + "\n"
    
    with open(output_file, 'w') as f:
        f.write(header)
        
        for result in all_results:
            row_parts = [
                result[0],  # sample
                result[1],  # gene_name
                str(result[2]),  # gene_count
                str(result[3]),  # total_mapped
                result[4]   # ratio
            ]
            
            if args.include_mismatch:
                row_parts.append(str(result[5]))  # max_mismatch
            if args.include_gene_file:
                row_parts.append(result[6])  # gene_file
            
            f.write("\t".join(row_parts) + "\n")
    
    print(f"\nProcessing completed!")
    print(f"Total samples: {len(samples)}")
    print(f"Gene list files used: {', '.join(sorted(gene_files_used))}")
    print(f"Total unknown sequence types: {unknown_total}")
    print(f"Results saved to: {output_file}")
    
    # Generate summary statistics
    summary_file = os.path.join(output_dir, "summary.txt")
    with open(summary_file, 'w') as f:
        f.write("Gene Statistics Summary Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Working directory: {workdir}\n")
        f.write(f"Number of samples: {len(samples)}\n")
        f.write(f"NC group samples: {len(nc_samples)}\n")
        f.write(f"Treatment group samples: {len(treat_samples)}\n")
        f.write(f"Gene list files used: {', '.join(sorted(gene_files_used))}\n")
        f.write(f"Maximum allowed mismatches: {args.mismatch}\n")
        f.write(f"Unknown sequence types: {unknown_total}\n")
        
        # Calculate detected gene ratios
        all_genes = set()
        nc_genes = set()
        treat_genes = set()
        
        for result in all_results:
            gene_name = result[1]
            count = result[2]
            if count > 0:
                all_genes.add(gene_name)
                if result[0].startswith('NC'):
                    nc_genes.add(gene_name)
                else:
                    treat_genes.add(gene_name)
        
        # Load all genes to get total count
        all_gene_files = []
        if os.path.exists(gene_nc_file):
            all_gene_files.append(gene_nc_file)
        if os.path.exists(gene_treat_file):
            all_gene_files.append(gene_treat_file)
        
        total_gene_count = 0
        for gene_file in all_gene_files:
            try:
                with open(gene_file, 'r') as gf:
                    lines = gf.readlines()
                    total_gene_count += len([l for l in lines if l.strip() and not l.startswith('#')])
            except:
                pass
        
        if total_gene_count > 0:
            f.write(f"Total genes detected: {len(all_genes)} ({len(all_genes)/total_gene_count*100:.1f}%)\n")
            f.write(f"NC group detected genes: {len(nc_genes)}\n")
            f.write(f"Treatment group detected genes: {len(treat_genes)}\n")
    
    print(f"Summary report: {summary_file}")

if __name__ == "__main__":
    main()