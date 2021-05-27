import os
import numpy as np
import pandas as pd
from glob import glob

# get file name bases
files = glob(snakemake.input.fna + '*.fna')
print(len(np.unique(files)))
base_names = []
for file in files:
    before_ext = file.split('/')[-1].split('.')[:-1]
    base_name = '.'.join(before_ext)
    base_names.append(base_name)

# I'll use pandas to save my name lookup table
df_rows = []
for bn in base_names:
    out_bn = bn.replace('_genomic', '')

    # pull the first line in and split it into words. We want genus_species
    # when possible but sometimes we'll settle for genus_sp_code
    with open(snakemake.input.fna + bn + '.fna', 'r') as fin:
        first_line = fin.readline()
    first_words = first_line.split(' ')
    if first_words[2] == 'sp.':
        clean_name = '_'.join([first_words[1], 'sp', first_words[3]])
    elif first_words[1][0] == '[':
        clean_name = '_'.join(first_words[2:5]).replace('_strain', '')
    else:
        clean_name = '_'.join(first_words[1:4]).replace('_strain', '')

    row = {'acc_name': out_bn, 
           'clean_name': clean_name.replace(',', '').lower()}

    df_rows.append(row)

    # now we'll do the actually renamining with just shell commands I think
    os.system('mv {0} {1}'.format(snakemake.input.fna + bn + '.fna', 
                                  snakemake.input.fna + clean_name + '.fna'))
    os.system('mv {0} {1}'.format(snakemake.input.faa + out_bn + '_protein.faa',
                                  snakemake.input.faa + clean_name + '.faa'))

out_df = pd.DataFrame(df_rows)
out_df.to_csv('data/genome_list.csv')
print(out_df['clean_name'].unique().shape[0])