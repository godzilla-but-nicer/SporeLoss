import numpy as np
import pandas as pd
from Bio import AlignIO, Seq

# parameter to determine the maximum missing proportion that we keep
missing_thresh = 0.4

# load the alignments and turn them into a numpy array
alignments = AlignIO.read(snakemake.input[0], 'fasta')
align_arr = np.array([list(rec) for rec in alignments])

# get a list of missing values per base
missing_bases = []

# iterate over the whole alignment counting missing bases
for base in range(align_arr.shape[1]):
    missing = 0
    for seq in range(align_arr.shape[0]):
        if alignments[seq, base] not in ['A', 'T', 'G', 'C']:
            missing += 1
    missing_bases.append(missing)

# calculate the proportion of missing bases for each column
missing_prop = np.array([m / align_arr.shape[0] for m in missing_bases])
align_arr = align_arr[:, missing_prop < missing_thresh]

for r, rec in enumerate(alignments):
    joined_seq = ''.join(align_arr[r])
    print(joined_seq[:10])
    rec.seq = Seq.Seq(joined_seq)

with open(snakemake.output[0], 'w') as fout:
    AlignIO.write(alignments, fout, 'fasta')
