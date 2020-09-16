import pandas as pd
from Bio import SeqIO

# faa extension
fext = '.faa'

# this will become a dataframe for file writing
dict_list = []

with open(snakemake.input.ref, 'r') as ref:
    gene_labels = ref.read().split('\n')[:-1]

# for each gene we'll grab first entry in the label because thats what 
# blast output uses
for gene in gene_labels:
    # new row for the eventual dataframe
    new_row = {}
    new_row['name'] = gene

    # get the file and pull out the label
    faa = snakemake.input.faas + gene + fext
    rec = SeqIO.read(faa, 'fasta')
    new_row['label'] = rec.name

    dict_list.append(new_row)

names_df = pd.DataFrame(dict_list)

names_df.to_csv(snakemake.output[0])
