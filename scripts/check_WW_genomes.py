import pandas as pd
from glob import glob

# load the data and get the path to the genomes
gen_path = 'data/weller_wu_labelled/genomes/'
wwdf = pd.read_csv('data/WW_ascension_taxa.csv')

# figure out which genomes are already downloaded
found_genomes = [path.split('/')[-1] for path in glob(gen_path + '*.fna')]
found_ids = [gen.split('.fna')[0] for gen in found_genomes]
print(found_ids[:5])

# pull out taxon ids
ascs = wwdf['assembly_id']

# remove the genomes already downloaded from ascs
ascs = [asc for asc in ascs if asc not in found_ids]

asc_dict = {'assembly_id' : ascs}
out_df = pd.DataFrame(asc_dict)
out_df.to_csv('data/missing_ww.csv')