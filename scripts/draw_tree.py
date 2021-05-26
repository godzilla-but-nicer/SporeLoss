from ete3 import Tree, TreeStyle
import pandas as pd

# load jordan predictions
jor = pd.read_csv(snakemake.input.jor, index_col=0)
sample = jor.sample(100)
sample_genomes = sample['genome']
root = 'NC_000913.3'

tree = Tree(snakemake.input.tree)
tree.set_outgroup(root)

# keep_list = []
# for n, node in enumerate(tree.get_leaf_names()):
#     if node in sample_genomes:
#         keep_list.append(node)
# 
# tree.prune(keep_list)

cs = TreeStyle()
cs.mode = 'c' # draw tree in circular mode
cs.arc_start = -180
cs.arc_span = 180
cs.show_leaf_name = True
tree.render("mytree.pdf", w=183, units="mm", tree_style=cs)
