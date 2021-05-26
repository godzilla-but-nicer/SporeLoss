from Bio import AlignIO
from Bio.Alphabet import Gapped, generic_dna
from tqdm import tqdm
import pandas as pd

# load the start stop info
df = pd.read_csv(snakemake.input.txt)

# get the gene name from snakemake
gene_name = snakemake.wildcards.gene
gene_row = df[df['gene'] == gene_name]
start = int(gene_row[' start'] - 1)
stop = int(gene_row[' stop'])

print('trimming {0:s} to: [{1:d}, {2:d})'.format(gene_name, start, stop))

# load the sequence records set up alphabet
seqs = AlignIO.read(snakemake.input.afa, 'fasta')
alphabet = Gapped(generic_dna)

# iterate over the sequences and trim
for record in tqdm(seqs):
    record.seq = record.seq[start:stop]

AlignIO.write(seqs, snakemake.output[0], 'fasta')
