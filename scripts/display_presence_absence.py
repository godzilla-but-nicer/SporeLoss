import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# load data
ww = pd.read_csv(snakemake.input.ww, index_col=0)
labs = pd.read_csv(snakemake.input.labs, index_col=0)

# fix column label on labs
labs = labs.rename({'assembly_id': 'genome'}, axis=1)
ww = ww.merge(labs[['taxa_name', 'genome', 'spore_forming']], how='inner')

print(ww['spore_forming'])

ww_sorted = ww.sort_values('spore_forming', axis=0)
ww_sorted = ww_sorted.set_index('taxa_name')
ww_plot = ww_sorted.drop(['genome', 'spore_forming'], axis=1)
sns.heatmap(ww_plot, cbar=False)
plt.tight_layout()
plt.savefig(snakemake.output.mat)