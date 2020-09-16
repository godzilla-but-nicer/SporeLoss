from glob import glob

anno_path = snakemake.input.gff
geno_path = snakemake.input.fna

gffs = glob(anno_path + '*.gff')
fastas = glob(geno_path + '*.fna')

# use to compare files
clean_gffs = []
clean_fastas = []

# we want to see if the file names without extension match for GFFs and fastas
for gff in gffs:
    gff_name_list = gff.split('/')[-1].split('.')[:-1]
    name = '.'.join(gff_name_list)
    clean_gffs.append(name)

for fasta in fastas:
    fasta_name_list = fasta.split('/')[-1].split('.')[:-1]
    name = '.'.join(fasta_name_list)
    clean_fastas.append(name)

print(sorted(clean_gffs)[:5])
print(sorted(clean_fastas)[:5])

if sorted(clean_gffs) == sorted(clean_fastas):
    nice_format = '\n'.join(clean_fastas)
    with open(snakemake.output[0], 'w') as fout:
        fout.write(nice_format)
else:
    print('File Mismatch@!!')
