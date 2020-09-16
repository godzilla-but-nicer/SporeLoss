from Bio import SeqIO

records = SeqIO.parse(snakemake.input.fna, 'fasta')
for rec in records:
    # translate with bacterial code and give it some labels
    prot = rec.translate(table=11, to_stop=True)
    if len(prot) < 10:
        prot = rec.reverse_complement().translate(table=11, to_stop=True)
    prot.id = rec.id
    prot.name = rec.name
    prot.description = rec.description
    
    # write it to our output file
    with open(snakemake.output.faa, 'a') as fout:
        SeqIO.write(prot, fout, 'fasta')
