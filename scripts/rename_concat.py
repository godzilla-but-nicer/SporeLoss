import pandas as pd
from tqdm import tqdm
from Bio import AlignIO, SeqIO

# assign input files
genome_dir = snakemake.input.gen_dir
concat = open(snakemake.input.con, 'r')

with open(snakemake.input.gen_list) as fin:
    genomes = fin.readlines()
    genomes = [genome.strip() for genome in genomes]

# list of dictionaries to become df
tree_key = []

for g, genome in tqdm(enumerate(genomes)):
    new_name = genome.replace('_genomic', '_protein')

    # we need a bunch of entries corresponding to each assembly name
    # because we have many contigs (scaffolds?) in each genome
    for rec in SeqIO.parse(genome_dir + genome + '.fna', 'fasta'):
        contig_name = rec.name
        list_desc = rec.description.split(' ')
        species = ' '.join(list_desc[1:3])
        entry = {'downstream_name': new_name, 'contig_name': contig_name, 
                 'species': species}
        tree_key.append(entry)
    
# write out the dataframe in case I need to refer to these ids later
rename_df = pd.DataFrame(tree_key)
rename_df.to_csv(snakemake.output.csv)

# now we have to actually rename the labels on the concatenated ribosomal genes
out_align = open(snakemake.output.renamed, 'w')
aligns = AlignIO.read(concat, 'fasta')

# counters
count = 0
name_change = 0

for rec in tqdm(aligns[:-1]):
    count += 1
    relabel = rename_df[rename_df['contig_name'].str.contains(rec.name)]['downstream_name']
    if len(relabel) > 0:
        name_change += 1
        rec.description = relabel.values[0]
        rec.id = relabel.values[0]

print('Changed {}% of names'.format(name_change / count * 100))
    

AlignIO.write(aligns, out_align, 'fasta')
out_align.close()
concat.close()
