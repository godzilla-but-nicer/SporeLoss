from Bio import SeqIO
import re
import gffutils

# these would come from snakemake
full_gff_path = "data/Jordan_MAG_list/GCF_000005845.2_ASM584v2_genomic.ecoli_outgroup.gff"
full_fasta_path = "data/Jordan_MAG_list/GCF_000005845.2_ASM584v2_genomic.ecoli_outgroup.fna"

# protein name prefix (different for this genome, sick standard)
prefix = 'ribosomal subunit protein '

# would be piped in from config
gene_labels = ['L2', 'L3', 'L4', 'L5', 'L14', 'L16', 'L18', 'L22', 'L22', 'L24',
               'S3', 'S8', 'S10', 'S17', 'S18']

# these dicts will keep things for later
start_stop = {}

# Open GFF sqlite3 database, name doesn't matter, it will be deleted always
with open(full_gff_path) as fin:
    db = gffutils.create_db(full_gff_path, dbfn='data/Jordan_MAG_list/test.db',
                            force=True, merge_strategy='merge')
db = gffutils.FeatureDB('data/Jordan_MAG_list/test.db')

# pull out contig id, start, and stop coords for each gene
for gene in db.all_features(featuretype='CDS'):
    for target in gene_labels:
        # we need to
        regex = prefix + target + r'\D*$'
        if re.search(regex, gene.attributes['product'][0]):
            start_stop[target] = (gene.seqid, gene.start, gene.end)

print(start_stop)

# if we don't find all of the ribosomal genes than our genome is probably
# not good enough
if set(gene_labels) == start_stop.keys():
    print('we in here')

    # load up the sequences and extract the specific genes
    sequences = SeqIO.to_dict(SeqIO.parse(full_fasta_path, 'fasta'))
    for gene in start_stop:
        coords = start_stop[gene]
        gene_rec = sequences[coords[0]][coords[1]:coords[2]]
        with open('data/Jordan_MAG_list/ribosomal_genes/ecoli_'+gene+'.fna',
                  'a') as fo:
            SeqIO.write(gene_rec, fo, 'fasta')