#!/usr/bin/env python3
import pandas as pd
import argparse

def get_best_hits(blast_file):
    cols = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df = pd.read_csv(blast_file, sep='\t', names=cols)
    return df.sort_values(by=['evalue', 'bitscore'], ascending=[True, False]).drop_duplicates(subset=['qseqid'])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query_to_ref', required=True)
    parser.add_argument('-r', '--ref_to_query', required=True)
    parser.add_argument('-o', '--output', required=True)
    args = parser.parse_args()

    forward_hits = get_best_hits(args.query_to_ref)
    reverse_hits = get_best_hits(args.ref_to_query)

    # Swap query and subject columns in the reverse hits to match the forward format
    rev_pairs = reverse_hits[['qseqid', 'sseqid']].rename(columns={'qseqid': 'sseqid', 'sseqid': 'qseqid'})
    
    # Merge to filter the forward hits down to ONLY the reciprocal pairs
    rbh = pd.merge(forward_hits, rev_pairs, on=['qseqid', 'sseqid'])
    
    # Save without headers to maintain the exact BLAST format needed for plotting
    rbh.to_csv(args.output, sep='\t', index=False, header=False)

if __name__ == '__main__':
    main()
