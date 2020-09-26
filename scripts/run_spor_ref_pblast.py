import os
import glob

# set up some strings we'll use down the line
fext = '.faa'
dbext = '.pbdb'
boutext = '.txt'

# open the file with the list of ref genes, assemble list of paths
with open(snakemake.input.ref, 'r') as ref_file:
    ref_genes = [ref.strip() for ref in ref_file.readlines()]
ref_paths = [snakemake.input.ref_dir + r + fext for r in ref_genes]

# open the file with the list of genomes, assemble a list of paths
with open(snakemake.input.prot, 'r') as prot_file:
    prot_gen = [prot.strip() for prot in prot_file.readlines()]
prot_paths = [snakemake.input.prot_dir + p + fext for p in prot_gen]


# first we make blast dbs out of the sequences
dbs = []
for p, ppath in enumerate(prot_paths):
    pdb_path = snakemake.output.db + prot_gen[p] + dbext
    dbs.append(pdb_path)
    os.system('makeblastdb -in {0} -out {1} -dbtype prot'.format(ppath, pdb_path))


for g, gpath in enumerate(dbs):
    genome_dir = snakemake.output.out + prot_gen[g] + '/'
    os.mkdir(genome_dir)
    for r, rpath in enumerate(ref_paths):
        bout_path = genome_dir + ref_genes[r] + boutext
        os.system('blastp -query {0} -db {1} -out {2} -num_threads 4 -outfmt 6 -evalue 0.001'.format(rpath, gpath, bout_path))

    blast_outs = glob.glob(genome_dir + '*')
    for blastf in sorted(blast_outs):
        os.system('cat {0} >> {1}'.format(blastf, snakemake.output.out + prot_gen[g] + '.blastout'))
    os.system('rm -r {0}'.format(genome_dir))
