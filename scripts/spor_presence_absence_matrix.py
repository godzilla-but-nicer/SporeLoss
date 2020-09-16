import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


# first get our headers
with open(snakemake.input.header, 'r') as head_fin:
    header = head_fin.read().split(' ')

# next load the dataframe with the relabelling
name_file = snakemake.input.blast_dir + 'blastout_genes.csv'
relabels = pd.read_csv(name_file, index_col=0)

# this list of dictionaries will become our presence absence matrix
dict_list = []

# iterate through all of the blast output files
blast_files = glob(snakemake.input.blast_dir + '*.blastout')
for fi, blastf in enumerate(blast_files):
    new_row = {}
    # load the file as df
    blast_df = pd.read_csv(blastf, sep='\t', names=header, index_col=False)
    # we're going to make a new dataframe with rows as genomes with vectors of
    # gene presence-absence bits
    genome = blastf.split('/')[-1].split('.')[0]
    new_row['genome'] = genome
    for _, name, label in relabels.itertuples():
        if label in blast_df['qseqid'].values:
            new_row[name] = 1
        else:
            new_row[name] = 0
    dict_list.append(new_row)

# save this dictionary
pa_mat = pd.DataFrame(dict_list)
pa_mat.to_csv(snakemake.output.mat)

# let's do the PCA plot
pca = PCA(n_components=2)
X = pa_mat.values[:, 1:]
X_pca = pca.fit_transform(X)

fig, ax = plt.subplots()
ax.scatter(X_pca[:, 0], X_pca[:, 1])
ax.set_xlabel('PC1 ({:.2f}%)'.format(pca.explained_variance_[0]))
ax.set_ylabel('PC2 ({:.2f}%)'.format(pca.explained_variance_[1]))
fig.savefig(snakemake.output.pca)

# gene commonness
num_genomes = pa_mat.shape[0]
gene_freqs = np.sum(X, axis=0) / num_genomes

# get labels
gene_labels = pa_mat.columns[1:]

# get sorted indices by gene_frequencies
sort_i = np.argsort(gene_freqs)

fig, ax = plt.subplots()
ax.plot(np.arange(gene_freqs.shape[0]), gene_freqs[sort_i[::-1]])
ax.set_xticklabels(gene_labels[sort_i[::-1]], rotation=90)
ax.set_ylabel('Gene Frequency')
fig.savefig(snakemake.output.gene_dist)
