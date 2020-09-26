import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler

# load the dataframes and attach labels to the weller and wu set
ww = pd.read_csv(snakemake.input.ww, index_col=0)
jor = pd.read_csv(snakemake.input.jor, index_col=0)
labs = pd.read_csv(snakemake.input.labs, index_col=0)
genes = pd.read_csv(snakemake.input.genes, index_col=0)


labels = labs[['assembly_id', 'spore_forming']]
labels = labels.rename({'assembly_id' : 'genome'}, axis=1)

ww_labels = pd.merge(ww, labels, on='genome')

# copy each of the dataframes for each of our network features
ww_betw = ww_labels.copy()
ww_pg = ww_labels.copy()
ww_deg = ww_labels.copy()

jor_betw = jor.copy()
jor_pg = jor.copy()
jor_deg = jor.copy()

# add a little bit of noise to the network properties
genes['betweenness'] += (genes['betweenness'] + 
                         np.random.normal(scale=1e-3, size=genes.shape[0]))
genes['Outdegree'] += (genes['Outdegree'] + 
                         np.random.normal(scale=1e-3, size=genes.shape[0]))
genes['pagerank'] += (genes['pagerank'] + 
                         np.random.normal(scale=1e-3, size=genes.shape[0]))



# scale the dataframes by the network properties
# there are about 50 genes in the genome that do not appear in our
# network properties list. Thats actually fine we just need to be aware
for col in jor.drop('genome', axis=1).columns:
    if genes['geneName'].str.contains(col).any():
        # betweenness
        ww_betw[col] *= genes.loc[genes['geneName'] == col, 'betweenness'].values
        jor_betw[col] *= genes.loc[genes['geneName'] == col, 'betweenness'].values

        # pagerank
        ww_pg[col] *= genes.loc[genes['geneName'] == col, 'pagerank'].values
        jor_pg[col] *= genes.loc[genes['geneName'] == col, 'pagerank'].values

        # degree
        ww_deg[col] *= genes.loc[genes['geneName'] == col, 'Outdegree'].values
        jor_deg[col] *= genes.loc[genes['geneName'] == col, 'Outdegree'].values
    else:
        # drop genes we can't scale
        ww_betw.drop(col, axis=1, inplace=True)
        jor_betw.drop(col, axis=1, inplace=True)
        ww_pg.drop(col, axis=1, inplace=True)
        jor_pg.drop(col, axis=1, inplace=True)
        ww_deg.drop(col, axis=1, inplace=True)
        jor_deg.drop(col, axis=1, inplace=True)

print(ww_betw.head())


# going to define a function that does all of the stuff
def df_to_umap(labelled_df, unlabelled_df, net_prop):
    # Get all of this into array form
    y = labelled_df['spore_forming'].values
    X_ww = labelled_df.drop(['spore_forming', 'genome'], axis=1).values
    ww_names = labelled_df['genome'].values
    
    X_jor = unlabelled_df.drop('genome', axis=1).values
    jor_names = unlabelled_df['genome'].values

    # standardize
    ss = StandardScaler()
    jor_std = ss.fit_transform(X_jor)
    ww_std = ss.transform(X_ww)

    # fit UMAP
    um = umap.UMAP()
    jor_umap = um.fit_transform(jor_std)
    ww_umap = um.transform(ww_std)

    # plot the reduced samples
    fig, ax = plt.subplots()
    ax.scatter(jor_umap[:, 0], jor_umap[:, 1])

    legend_lab = 'Spore Forming: '
    for yn in np.unique(y):
        these_samples = ww_umap[y == yn, :]
        ax.scatter(these_samples[:, 0], these_samples[:, 1], label=legend_lab + yn)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.legend()

    return ax


# betweenness scaled plot
fig, ax = plt.subplots()
ax = df_to_umap(ww_betw, jor_betw, 'betweenness')
plt.savefig(snakemake.output.umap_betw)

# pagerank
fig, ax = plt.subplots()
ax = df_to_umap(ww_pg, jor_pg, 'pagerank')
plt.savefig(snakemake.output.umap_pg)

fig, ax = plt.subplots()
ax = df_to_umap(ww_deg, jor_deg, 'degree')
plt.savefig(snakemake.output.umap_deg)