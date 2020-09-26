import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# load umap data
umap_ww = pd.read_csv(snakemake.input.ww, index_col=0)
umap_jor = pd.read_csv(snakemake.input.jor, index_col=0)

# extract the jor data as a numpy array and cluster
jor = umap_jor.drop('genome', axis=1).values

km = KMeans(n_clusters=2)
km.fit(jor)

# predictions for the jor genomes
jor_predictions = km.predict(jor)
umap_jor['predicted_2'] = jor_predictions
umap_jor.to_csv(snakemake.output.classes)

# we want to figure out where our decision regions are
# so we need to make a grid
step = 0.005
x_min, x_max = jor[:, 0].min() - 10*step, jor[:, 0].max() + 10*step
y_min, y_max = jor[:, 1].min() - 10*step, jor[:, 1].max() + 10*step
x_grid = np.arange(x_min, x_max, step)
y_grid = np.arange(y_min, y_max, step)
xx, yy = np.meshgrid(x_grid, y_grid)

# get our grid into the right shape to make predictions
# then reshape it for plotting
grid_predict = np.zeros((xx.ravel().shape[0], 2))
grid_predict[:, 0] = xx.ravel()
grid_predict[:, 1] = yy.ravel()

Z = km.predict(grid_predict)
Z = Z.reshape(xx.shape)

print(km.cluster_centers_)

# plot this bad boy
fig, ax = plt.subplots()
ax.contourf(xx, yy, Z, alpha = 0.3, cmap='viridis')
ax.scatter(jor[:, 0], jor[:, 1], marker='o', facecolors='none',
           edgecolor='grey', alpha = 0.8)
ax.scatter(km.cluster_centers_[:, 0], km.cluster_centers_[:, 1], marker='+',
           s=300, linewidth=0.8, color='black')
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')

fig.savefig(snakemake.output.descision)
