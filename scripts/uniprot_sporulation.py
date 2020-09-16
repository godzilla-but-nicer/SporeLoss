import pandas as pd
import re
from Bio import SeqIO

# load the subtiwiki gene list
genes = pd.read_csv(snakemake.input.genes)
gene_names = set(genes['locus tag'].values)

# regular expressions and other strings for matching gene and species names
gene_re = r'\WGN=([a-zA-Z0-9\-\_]+)\W'
species = 'Bacillus subtilis (strain 168)'

# load sequence records
for r, rec in enumerate(SeqIO.parse(snakemake.input.faa, 'fasta')):
    desc = rec.description
    matches = re.findall(gene_re, desc)
    if len(matches) > 0:
        match_name = matches[0]
    else:
        match_name = ''
    if match_name in gene_names and species in desc:
        with open(snakemake.output[0] + match_name + '.faa', 'w') as fout:
            SeqIO.write(rec, fout, 'fasta')
