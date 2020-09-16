import os
import pandas as pd
import numpy as np
from glob import glob

# file extensions we will pass into our commands
zext = '.zip'
pext = '.faa'
aext = '.gff'
gext = '.fna'

# load the WW dataframe
wwdf = pd.read_csv(snakemake.input.csv)

# figure out which genomes are already downloaded
found_genomes = [path.split('/')[-1] for path in glob(snakemake.output.gen_dir + '*.fna')]
found_ids = [gen.split[0] for gen in found_genomes]

# pull out taxon ids
ascs = wwdf['assembly_id']

# remove the genomes already downloaded from ascs
ascs = [asc for asc in ascs if asc not in found_ids]

# this list will track taxa that have more than one associated genome
multiple_genomes = []

# iterate through these and download them
for ti, asc in enumerate(ascs):
    # download and unzip the data
    zip_name = str(asc) + zext
    os.system(
        './datasets download assembly {0} -pgf {1}'.format(asc, zip_name))
    os.system('unzip {0} -x {1}'.format(zip_name,
                                        'README.md ncbi_dataset/data/'
                                        'dataset_catalog.json ncbi_dataset/'
                                        'fetch.txt'))
    # get some paths and names to organize things later
    genomes = glob('ncbi_dataset/data/GCF*')
    # this is a problem because I'm literally just grabbing which ever genome
    # happens to be first. They're all from refseq and should be good but could
    # be genomes for strains with different functional behavior
    if len(genomes) < 1:
        print('skipping ascension {}! No genomes found'.format(asc))
        genome_name = np.nan
    else:
        if len(genomes) > 1:
            print('Warning! Multiple genomes for ascension {}!'.format(asc))
        genome_path = genomes[0]
        genome_name = genome_path.split('/')[-1]
        print(genome_name)
        # move things into a structure that matches the MAG directory
        os.system('mv {0} {1}'.format(genome_path + '/protein.faa',
                                      (snakemake.output.prot_dir + genome_name
                                       + pext)))
        os.system('mv {0} {1}'.format(genome_path + '/genomic.gff',
                                      (snakemake.output.annot_dir + genome_name
                                       + aext)))
        os.system('mv {0} {1}'.format(genome_path + '/*.fna',
                                      (snakemake.output.gen_dir + genome_name
                                       + gext)))
    # cleanup
    os.system('rm -r ncbi_dataset')
    os.system('rm {}.zip'.format(asc))

print('ok')
