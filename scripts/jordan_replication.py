import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# load and combine dataframes
pam = pd.read_csv(
    'data/Jordan_MAG_list/spore_prediction/presence_absence.csv', index_col=0)

labs = pd.read_csv(
    'data/Jordan_MAG_list/spore_prediction/jor_predicted_classes_2.csv',
    index_col=0).drop(['umap_0', 'umap_1'], axis=1)

df = pam.merge(labs, on='genome')

# pull out presence absence matrix and labels
mat = pam.drop('genome', axis=1).values
labs = df['predicted_2'].values

# get columns and row sums
col_sum = np.sum(mat, axis=0)
num_genes = np.sum(mat, axis=1)

# we'll use these bins a lot of times
bin_width = 12
bin_edges = np.arange(0, num_genes.max() + bin_width, bin_width)

# plot overall histogram of gene counts
hist, edges = np.histogram(num_genes, bins=bin_edges)
plt.bar(edges[1:], hist, width=20)
plt.xlabel('Number of Sporulation Genes')
plt.ylabel('Count')
plt.savefig('overall_gene_count.png')

# then do histogram of gene counts for the sporulators and nonsporulators
spor_mat = mat[labs == 0, :]
spor_num = np.sum(spor_mat, axis=1)
non_mat = mat[labs == 1, :]
non_num = np.sum(non_mat, axis=1)

spo_hist, spo_edges = np.histogram(spor_num, bins=bin_edges)
non_hist, non_edges = np.histogram(non_num, bins=bin_edges)
plt.figure()
plt.bar(spo_edges[1:], spo_hist, label='Sporulators', width=bin_width, alpha=0.6)
plt.bar(non_edges[1:], non_hist, label='Nonsporulators', width=bin_width, alpha=0.6)
plt.xlabel('Number of Sporulation Genes')
plt.ylabel('Count')
plt.legend()
plt.savefig('type_gene_count.png')
