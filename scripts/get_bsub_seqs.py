import pandas as pd
import numpy as np
import gffutils
from Bio import SeqIO

# load subtiwiki gene list
df = pd.read_csv(snakemake.input.genes)
lts = df['gene'].values
names = df['locus tag'].values

# make arrays that will hold the starts and stops of all of the genes
starts = np.zeros(lts.shape[0], dtype=np.int)
stops = np.zeros(lts.shape[0], dtype=np.int)

# open gff db
db = gffutils.create_db(snakemake.input.gff, dbfn=snakemake.output.db, force=True, merge_strategy='merge')
db = gffutils.FeatureDB(snakemake.output.db)

# iterate over the database and find the genes tagged with sporulation
for gene in db.all_features(featuretype='gene'):
    locus_tag = gene.attributes['locus_tag'][0]
    if locus_tag in lts:
        gene_i = np.where(lts == locus_tag)
        starts[gene_i] = gene.start
        stops[gene_i] = gene.end

# for now we'll remove any genes that did not have matching entries in the
# GFF
keep_idx = np.where(stops != 0)
lts = lts[keep_idx]
starts = starts[keep_idx]
stops = stops[keep_idx]
names = names[keep_idx]

# load the genome
genome = SeqIO.to_dict(SeqIO.parse(snakemake.input.fna, 'fasta'))
for key in genome.keys():
    for t, tag in enumerate(lts):
        gene_seq = genome[key][starts[t]:stops[t]]
        gene_seq.id = tag
        gene_seq.name = names[t]
        gene_seq.description = names[t] + ' NCBI NC_000964.3'
        with open(snakemake.output.fna, 'a') as fout:
            SeqIO.write(gene_seq, fout, 'fasta')
