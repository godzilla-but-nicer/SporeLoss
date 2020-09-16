from Bio import AlignIO


input_files = snakemake.input
names = [f.split('/')[-1].split('.')[0] for f in input_files]
print(names)

# open alignment files
afa = {}
for gene, fin in zip(names, input_files):
    afa[gene] = AlignIO.read(fin, 'fasta')

# loop over each entry for the first gene
for i, record in enumerate(afa[names[0]]):
    # and in that loop over all of the other genes
    for gene in names[1:]:
        record.seq += afa[gene][i].seq

# save it
AlignIO.write(afa[names[0]], snakemake.output[0], 'fasta')
