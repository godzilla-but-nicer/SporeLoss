import pandas as pd
import re

# load weller and wu's list of sporulating genomes
WW_ref = pd.read_excel(snakemake.input.xlsx)

with open(snakemake.input.afa, 'r') as afa_fin:
    lines = afa_fin.readlines()


def get_asc_and_name(line):
    # get the taxa name
    inside_braces = re.findall(r'\{(.+)\}', line)[0]
    substrings = inside_braces.split(' ')
    # with unnamed species we have 'sp.' in the second entry
    # but these are used in the WW list
    if substrings[1] == 'sp.':
        name = ' '.join(substrings[:3])
    else:
        name = ' '.join(substrings[:2])

    # get ascension number
    after_dash = line.split('-')[1]
    asc = after_dash.replace('/', '.')

    return (asc, name)

WW_ref['ascension_number'] = '-1'
for line in lines:
    if line[0] == '>':
        asc, name = get_asc_and_name(line)
        WW_ref.loc[WW_ref['Taxon Name'].str.strip() == name.strip(), 'ascension_number'] = asc
    else:
        continue

# make a smaller cleaner dataframe with names that I like
rename = {'Taxon Name': 'taxa_name', 'NCBI Taxonomy ID': 'taxonomy_id',
          'Spore-Forming?': 'spore_forming'}
WW_out = WW_ref.drop(['Sporulation Gene Count',
                      'Synonymous Substution Rate (dS Sum)',
                      'Amino Acid Substitution Rate (Branch Length Sum)',
                      'Codon Bias Index'], axis=1)
WW_out = WW_out.rename(rename, axis=1)
print(WW_out[WW_out['ascension_number'] == '-1'])
# remove missing ascension numbers
WW_out = WW_out[WW_out['ascension_number'] != '-1']
WW_out.to_csv(snakemake.output.taxa)

# write out a simpler file with the query for NCBI
asc_nums = WW_out['ascension_number']
with_sep = asc_nums.apply(lambda x: x + ' OR\n')
with open(snakemake.output.asc, 'w') as asc_fout:
    for line in with_sep:
        asc_fout.write(line)