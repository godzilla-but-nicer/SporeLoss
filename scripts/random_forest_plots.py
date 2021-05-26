import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA

# load data
labels_df = pd.read_csv(snakemake.input.classes, index_col=0).drop(['umap_0', 'umap_1'], axis=1)
matrix_df = pd.read_csv(snakemake.input.mat, index_col=0)
genes = pd.read_csv(snakemake.input.genes, index_col=0)

# combine data and pull out relevant subsets for random forest
df = matrix_df.merge(labels_df, on='genome')
X = df.drop(['genome', 'predicted_2'], axis=1).values
y = df['predicted_2'].values
names = df.drop(['genome', 'predicted_2'], axis=1).columns

# train random forest, get ranked gene names by importance
rfc = RandomForestClassifier(n_estimators=1000, n_jobs=4, random_state=420)
rfc.fit(X, y)

# because we can lets show a decision plot after PCA (not just UMAP)
pca = PCA(2)
X_pca = pca.fit_transform(X)
exp_var = pca.explained_variance_ratio_ * 100

fig, ax = plt.subplots()
colors = ['C0', 'C1']
shapes = ['o', 's']
labels = ['Sporulator', 'Nonsporulator']
for coords, lab in zip(X_pca, y):
    ax.scatter(coords[0], coords[1], c=colors[lab], marker=shapes[lab], label=labels[lab])

ax.set_xlabel('PC1 ({:.3f}%)'.format(exp_var[0]))
ax.set_ylabel('PC2 ({:.3f}%)'.format(exp_var[1]))
ax.legend()
plt.savefig('plots/decision/pca_space_umap_labels.png')


important = rfc.feature_importances_
ranked = np.argsort(important)[::-1]
important_names = names[ranked]

rfc_df = pd.DataFrame({'geneName': names, 'rfc_import': important})
# ok one thing that is important is knowing whether genes are correlated
# lets make a big dumb correlation matrix
corr_array = df.drop(['genome', 'predicted_2'], axis=1).corr()
fig, ax = plt.subplots()
sns.heatmap(corr_array, center=0, ax=ax)
plt.savefig('plots/random_forest/corr_mat.png')

# variable to set how many point to plot in the scatter plots
n_pts = 609

# lets do outdegree with ranked gene importances
small_genes = genes[['Outdegree', 'betweenness', 'geneName', 'infomap', 'pagerank']].merge(rfc_df, on='geneName').set_index(['geneName'])
small_genes['rfc_nonzero'] = small_genes['rfc_import'].values > 0
ranked_genes = small_genes.reindex(important_names).dropna()

# Separate zero and nonzero importances for plotting
zeros = ranked_genes[ranked_genes['rfc_nonzero'] == False]
non = ranked_genes[ranked_genes['rfc_nonzero'] == True]

fig, ax = plt.subplots()
ax.scatter(zeros['rfc_import'], zeros['Outdegree'], c='grey')
ax.scatter(non['rfc_import'], non['Outdegree'])
# for i in range(0, 5):
#     ax.text(non.iloc[i]['rfc_import'], non.iloc[i]['Outdegree'], non.index[i], rotation=45)
ax.set_xlabel('Random Forest Importance')
ax.set_ylabel('Outdegree')
plt.savefig('plots/random_forest/reg_outdegree.png')

# betweenness too. These are the interpretable ones!
fig, ax = plt.subplots()
ax.scatter(zeros['rfc_import'], zeros['betweenness'], c='grey')
ax.scatter(non['rfc_import'], non['betweenness'])
ax.set_xlabel('Random Forest Importance')
ax.set_ylabel('Importance')
plt.savefig('plots/random_forest/reg_betw.png')

# infomap for funsies
fig, ax = plt.subplots(figsize=(10.5, 4))
ax.bar(ranked_genes.index[:50], ranked_genes['rfc_import'][:50], color='grey', width=0.1)
ax.scatter(ranked_genes.index[:50], ranked_genes['rfc_import'][:50])
ax.set_xticklabels(ranked_genes.index[:50], rotation=60)
ax.set_ylabel('Random Forest Importance')
plt.tight_layout()
plt.savefig('plots/random_forest/importance.png')

# ok one thing that is important is knowing whether genes are correlated
# lets make a big dumb correlation matrix for just the top 50 genes
# I think pearsons correlation is fine but there seems to be some debate
corr_array = df[ranked_genes.index[:50]].corr()
fig, ax = plt.subplots(figsize=(8,7))
sns.heatmap(corr_array, center=0, ax=ax)
plt.tight_layout()
plt.savefig('plots/random_forest/corr_mat.png')

# pagerank might be more informative than betweenness
fig, ax = plt.subplots()
ax.scatter(zeros['rfc_import'], zeros['pagerank'], c='grey')
ax.scatter(non['rfc_import'], non['pagerank'])
ax.set_xlabel('Random Forest Importance')
ax.set_ylabel('PageRank')
plt.tight_layout()
plt.savefig('plots/random_forest/reg_pagerank.png')
