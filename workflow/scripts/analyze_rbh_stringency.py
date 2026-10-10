import pandas as pd
import sys
from Bio import SeqIO

blast1_file = sys.argv[1]
blast2_file = sys.argv[2]
fasta_file = sys.argv[3]
out_file = sys.argv[4]

columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen']
df1 = pd.read_csv(blast1_file, sep='\t', names=columns)
df2 = pd.read_csv(blast2_file, sep='\t', names=columns)

total_predicted = sum(1 for _ in SeqIO.parse(fasta_file, "fasta"))

def calculate_rbh(d1, d2, min_id, min_cov):
    f1 = d1[(d1['pident'] >= min_id) & ((d1['length'] / d1['qlen']) >= min_cov) & ((d1['length'] / d1['slen']) >= min_cov)]
    f2 = d2[(d2['pident'] >= min_id) & ((d2['length'] / d2['qlen']) >= min_cov) & ((d2['length'] / d2['slen']) >= min_cov)]
    top1 = f1.sort_values(['qseqid', 'evalue', 'bitscore'], ascending=[True, True, False]).drop_duplicates('qseqid')
    top2 = f2.sort_values(['qseqid', 'evalue', 'bitscore'], ascending=[True, True, False]).drop_duplicates('qseqid')
    dict1 = pd.Series(top1.sseqid.values, index=top1.qseqid).to_dict()
    dict2 = pd.Series(top2.sseqid.values, index=top2.qseqid).to_dict()
    return sum(1 for q, s in dict1.items() if s in dict2 and dict2[s] == q)

thresholds = [
    (0.0, 0.00, "Base (E-value < 1e-5 only)"),
    (50.0, 0.50, "Moderate (50% Id, 50% Cov)"),
    (70.0, 0.70, "Strict (70% Id, 70% Cov)"),
    (90.0, 0.90, "Highly Stringent (90% Id, 90% Cov)")
]

with open(out_file, 'w') as f:
    f.write(f"Total BRAKER Predicted Proteins (Longest Isoforms): {total_predicted}\n\n")
    f.write(f"{'Stringency Criteria':<35} | {'RBH Count':<10} | {'No RBH':<10} | {'% with RBH':<10}\n")
    f.write("-" * 75 + "\n")
    for min_id, min_cov, label in thresholds:
        rbh_c = calculate_rbh(df1, df2, min_id, min_cov)
        no_rbh = total_predicted - rbh_c
        pct = (rbh_c / total_predicted) * 100
        f.write(f"{label:<35} | {rbh_c:<10} | {no_rbh:<10} | {pct:.2f}%\n")
