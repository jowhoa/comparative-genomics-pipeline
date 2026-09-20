#!/usr/bin/env python3
"""
calculate_theoretical_coverage.py
Computes yield, total reads, N50 read length, and theoretical sequencing coverage.
"""

import sys
import gzip
import json
import argparse
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate theoretical sequencing coverage from FASTQ.")
    parser.add_argument("--fastq", required=True, help="Path to raw/filtered input FASTQ file (can be gzipped)")
    parser.add_argument("--expected-size", type=int, required=True, help="Expected haploid genome size in bp")
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--output-tsv", required=True, help="Path to output TSV summary")
    parser.add_argument("--output-json", required=True, help="Path to output JSON summary")
    return parser.parse_args()

def calculate_n50(lengths, total_bases):
    lengths.sort(reverse=True)
    half_total = total_bases / 2.0
    cumulative = 0
    for length in lengths:
        cumulative += length
        if cumulative >= half_total:
            return length
    return 0

def main():
    args = parse_args()
    
    fastq_path = Path(args.fastq)
    open_func = gzip.open if fastq_path.suffix == ".gz" else open
    
    total_reads = 0
    total_bases = 0
    lengths = []
    
    print(f"Profiling {fastq_path}...", file=sys.stderr)
    
    with open_func(fastq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()
            
            if not qual:
                break
                
            read_len = len(seq)
            lengths.append(read_len)
            total_bases += read_len
            total_reads += 1

    if total_reads == 0:
        raise ValueError(f"Input FASTQ {args.fastq} contains zero reads.")

    mean_length = total_bases / total_reads
    n50_length = calculate_n50(lengths, total_bases)
    theoretical_coverage = total_bases / float(args.expected_size)

    metrics = {
        "sample_id": args.sample_id,
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mean_read_length": round(mean_length, 2),
        "n50_read_length": n50_length,
        "expected_genome_size_bp": args.expected_size,
        "theoretical_coverage": round(theoretical_coverage, 2)
    }

    with open(args.output_tsv, "w") as tsv_out:
        tsv_out.write("sample_id\ttotal_reads\ttotal_bases\tmean_length\tn50_length\texpected_genome_size\ttheoretical_coverage\n")
        tsv_out.write(f"{metrics['sample_id']}\t{metrics['total_reads']}\t{metrics['total_bases']}\t{metrics['mean_read_length']}\t{metrics['n50_read_length']}\t{metrics['expected_genome_size_bp']}\t{metrics['theoretical_coverage']}\n")

    with open(args.output_json, "w") as json_out:
        json.dump(metrics, json_out, indent=4)

    print(f"Analysis Complete: Total Bases: {total_bases:,} bp | Coverage: {theoretical_coverage:.2f}X", file=sys.stderr)

if __name__ == "__main__":
    main()
