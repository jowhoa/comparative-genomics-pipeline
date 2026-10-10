import sys
from Bio import SeqIO
import re

input_fasta = sys.argv[1]
output_fasta = sys.argv[2]
source_type = sys.argv[3]

gene_to_longest = {}
for record in SeqIO.parse(input_fasta, "fasta"):
    if source_type == 'braker':
        # Strips only the terminal '.t1', '.t2', etc. transcript suffix
        gene_id = re.sub(r'\.t\d+$', '', record.id)
    else:
        match = re.search(r'\[gene=(.*?)\]', record.description)
        gene_id = match.group(1) if match else re.sub(r'\.\d+$', '', record.id)

    if gene_id not in gene_to_longest:
        gene_to_longest[gene_id] = record
    else:
        if len(record.seq) > len(gene_to_longest[gene_id].seq):
            gene_to_longest[gene_id] = record

SeqIO.write(gene_to_longest.values(), output_fasta, "fasta")
