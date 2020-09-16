import pandas as pd

df = pd.read_csv(snakemake.input.table)
hg = df[df['metagenome_source'] == 'gut metagenome']
hg_sras = hg['Run'].values

with open(snakemake.output[0], 'w') as fout:
    for sra in hg_sras:
        fout.write('{}\n'.format(sra))