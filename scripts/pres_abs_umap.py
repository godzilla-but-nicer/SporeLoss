import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler

# load the dataframes and attach labels to the weller and wu set
ww = pd.read_csv(snakemake.input.ww, index_col=0)
jor = pd.read_csv(snakemake.input.jor, index_col=0)
labs = pd.read_csv(snakemake.input.labs, index_col=0)


labels = labs[['assembly_id', 'spore_forming']]
labels = labels.rename({'assembly_id' : 'genome'}, axis=1)

ww_labels = pd.merge(ww, labels, on='genome')

# Get all of this into array form
y = ww_labels['spore_forming'].values
X_ww = ww_labels.drop(['spore_forming', 'genome'], axis=1).values
ww_names = ww_labels['genome'].values

X_jor = jor.drop('genome', axis=1).values
jor_names = jor['genome'].values

# standardize
ss = StandardScaler()
jor_std = ss.fit_transform(X_jor)
ww_std = ss.transform(X_ww)

# fit UMAP
um = umap.UMAP()
jor_umap = um.fit_transform(jor_std)
ww_umap = um.transform(ww_std)

# save these reduced versions (because umap take long?)
ww_dict_umap = {'genome': ww_names, 'umap_0': ww_umap[:, 0],
               'umap_1': ww_umap[:, 1], 'spore_forming': y}
jor_dict_umap = {'genome': jor_names, 'umap_0': jor_umap[:, 0],
                'umap_1': jor_umap[:, 1]}

wwdf = pd.DataFrame(ww_dict_umap)
wwdf.to_csv(snakemake.output.ww_umap)
jordf = pd.DataFrame(jor_dict_umap)
jordf.to_csv(snakemake.output.jor_umap)

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

plt.savefig(snakemake.output.umap_plot)
