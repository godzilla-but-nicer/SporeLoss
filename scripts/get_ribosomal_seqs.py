from Bio import SeqIO
from tqdm import tqdm
import re
import gffutils

# really we would read a file with a list of these
with open(snakemake.input.genomes, 'r') as fin:
    genomes = [line.strip() for line in fin]

# these would come from snakemake
gff_path = snakemake.input.gff
fasta_path = snakemake.input.fna

gene_labels = ['L2', 'L3', 'L4', 'L5', 'L14', 'L16', 'L18', 'L22', 'L22', 'L24',
               'S3', 'S8', 'S10', 'S17', 'S18']
               
for i in tqdm(range(len(genomes))):
    start_stop = {}  # reset the gene coords
    # Open GFF sqlite3 database, name doesn't matter, it will be deleted always
    full_gff_path = gff_path + genomes[i] + '.gff'
    with open(full_gff_path) as fin:
        db = gffutils.create_db(full_gff_path, dbfn=snakemake.output.db,
                                force=True, merge_strategy='merge')
    db = gffutils.FeatureDB(snakemake.output.db)

    # pull out contig id, start, and stop coords for each gene
    for gene in db.all_features(featuretype='CDS'):
        for target in gene_labels:
            regex = 'ribosomal protein ' + target + r'\D*$'
            if re.search(regex, gene.attributes['product'][0]):
                start_stop[target] = (gene.seqid, gene.start, gene.end)

    # if we don't find all of the ribosomal genes than our genome is probably
    # not good enough
    if set(gene_labels) == start_stop.keys():

        with open('data/Jordan_MAG_list/successful.txt', 'a') as flog:
            flog.write(genomes[i] + '\n')

        # load up the sequences and extract the specific genes
        sequences = SeqIO.to_dict(SeqIO.parse(
            fasta_path + genomes[i] + '.fna', 'fasta'))
        for gene in start_stop:
            coords = start_stop[gene]
            gene_rec = sequences[coords[0]][coords[1]:coords[2]]
            with open('data/Jordan_MAG_list/ribosomal_genes/'+gene+'.fna',
                      'a') as fo:
                SeqIO.write(gene_rec, fo, 'fasta')
    
