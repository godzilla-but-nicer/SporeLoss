from os.path import join as j

configfile: "workflow/config.yaml"

DATA_DIR = config["data_dir"]

PAPER_DIR = config["paper_dir"]
PAPER_SRC, SUPP_SRC = [j(PAPER_DIR, f) for f in ("main.tex", "supp.tex")]
PAPER, SUPP = [j(PAPER_DIR, f) for f in ("main.pdf", "supp.pdf")]

RIBOSOMAL = config["ribosomal_genes"]

# HPC_DIR = config["hpc_dir"]
# FASTAs = config["fastas"]

rule all:
    input:
        PAPER, SUPP

rule paper:
    input:
        PAPER_SRC, SUPP_SRC
    params:
        paper_dir = PAPER_DIR
    output:
        PAPER, SUPP
    shell:
        "cd {params.paper_dir}; make"


# rule some_data_processing:
    # input:
        # "data/some_data.csv"
    # output:
        # "data/derived/some_derived_data.csv"
    # script:
        # "workflow/scripts/process_some_data.py"

rule get_file_names:
    input:
        gff="data/Jordan_MAG_list/annotations/",
        fna="data/Jordan_MAG_list/genomes/"
    output:
        "data/Jordan_MAG_list/genomes_list.txt"
    script:
        "scripts/get_file_names.py"

# Get ribosomal genes from B. subtilis ref genome
# this step is BIG and really should probably be broken up
rule get_ribosomal_fastas:
    input:
        gff="data/Jordan_MAG_list/annotations/",
        fna="data/Jordan_MAG_list/genomes/",
        genomes="data/Jordan_MAG_list/genomes_list.txt"
    output:
        expand("data/Jordan_MAG_list/ribosomal_genes/{gene}.fna", 
               gene=config["ribosomal_genes"]),
        db=temp("data/Jordan_MAG_list/test.db")
    script:
        "scripts/get_ribosomal_seqs.py"

# also for ribosomal gene alignment we're going to need an outgroup
rule add_ecoli_outgroup:
    input:
        gff="data/Jordan_MAG_list/GCF_000005845.2_ASM584v2_genomic.ecoli_outgroup.gff",
        fna="data/Jordan_MAG_list/GCF_000005845.2_ASM584v2_genomic.ecoli_outgroup.fna"
    output:
        "data/dummy.txt",
        db=temp("data/Jordan_MAG_list/test.db")
    script:
        "scripts/get_ribosomal_seqs_ecoli.py"


# alignment of genes so that we can build a phylogeny
rule multiple_align_ribosomal:
    input:
        "data/Jordan_MAG_list/ribosomal_genes/{gene}.fna"
    output:
        "data/Jordan_MAG_list/ribosomal_aligned/{gene}.afa"
    threads: 4
    shell:
        "clustalo --threads {threads} -i {input} -o {output}"

# requires manually edited trimming file
# I need to figure out a way to automate this
rule trim_ribosomal:
    input:
        afa="data/Jordan_MAG_list/ribosomal_aligned/{gene}.afa",
        txt="data/Jordan_MAG_list/ribosomal_aligned/manual_start_stop.txt"
    output:
        "data/Jordan_MAG_list/ribosomal_trimmed/{gene}.afa"
    script:
        "scripts/trim_aligned_genes.py"

rule concatenate_trimmed:
    input:
        expand("data/Jordan_MAG_list/ribosomal_trimmed/{gene}.afa",
               gene=RIBOSOMAL)
    output:
        "data/Jordan_MAG_list/tree_files/concat_ribosomal.afa"
    script:
        "scripts/concatenate_trimmed.py"

# the IDs in the sequences are not the ones we want to use for our tree
rule rename_seqs:
    input:
        gen_dir="data/Jordan_MAG_list/genomes/",
        gen_list="data/Jordan_MAG_list/genomes_list.txt",
        con="data/Jordan_MAG_list/tree_files/concat_ribosomal.afa"
    output:
        renamed="data/Jordan_MAG_list/tree_files/renamed_concat.afa",
        csv="data/Jordan_MAG_list/tree_files/tree_key.csv"
    script:
        "scripts/rename_concat.py"

# drop any site where 90% or more have missing base
rule clean_renamed_concat:
    input:
        "data/Jordan_MAG_list/tree_files/renamed_concat.afa"
    output:
        seq="data/Jordan_MAG_list/tree_files/renamed_concat_clean.afa",
        html="data/Jordan_MAG_list/tree_files/renamed_concat_clean.html"
    shell:
        "trimal -in {input} -out {output.seq} -htmlout {output.html} -strict"


rule fasttree:
    input:
        "data/Jordan_MAG_list/tree_files/renamed_concat_clean.afa"
    output:
        "data/Jordan_MAG_list/tree_files/concat_ribosomal.tree"
    shell:
        "fasttree -fastest -nt -gamma -gtr {input} > {output}"

# getting the sequences of sporulation genes from a reference B.subtilis genome
rule bsub_sporulation:
    input:
        fna='data/Jordan_MAG_list/GCF_000009045.1_ASM904v1_genomic.bsubtilisref.fna',
        gff='data/Jordan_MAG_list/GCF_000009045.1_ASM904v1_genomic.bsubtilisref.gff',
        genes='data/subtiwiki-sporulation-genes.csv'
    output:
        fna="data/ref_spor_genes/bsub_ref.fna",
        db=temp("data/Jordan_MAG_list/test.db")
    script:
        "scripts/get_bsub_seqs.py"

# an alternate plan to get these protein sequences. Search uniprot!
rule sporulation_proteins_uniprot:
    input:
        faa='data/uniprot_sprot.fasta',
        genes='data/subtiwiki-sporulation-genes.csv'
    output:
        directory('data/bsub_ref_spor_genes/')
    script:
        'scripts/uniprot_sporulation.py'

rule spor_protein_blast:
    input:
        ref='data/bsub_ref_spor_genes/ref_gene_list.txt',
        ref_dir='data/bsub_ref_spor_genes/',
        prot="data/Jordan_MAG_list/proteins_list.txt",
        prot_dir="data/Jordan_MAG_list/proteins/"
    output:
        db=directory('data/Jordan_MAG_list/pblastdbs/'),
        out=directory('data/Jordan_MAG_list/ref_spor_pblast_out/')
    script:
        'scripts/run_spor_ref_pblast.py'
    
# the blast output uses different labels for the gene names se we need to get
# them
rule blastout_gene_names:
    input:
        ref='data/bsub_ref_spor_genes/ref_gene_list.txt',
        faas='data/bsub_ref_spor_genes/'
    output:
        'data/Jordan_MAG_list/ref_spor_pblast_out/blastout_genes.csv'
    script:
        'scripts/get_blast_output_gene_labels.py'
    
# make a giant matrix for presence or absence of these genes
rule presence_absence_matrix:
    input:
        blast_dir='data/Jordan_MAG_list/ref_spor_pblast_out/',
        header='data/Jordan_MAG_list/blastout_header_fmt_6.txt'
    output:
        mat='data/Jordan_MAG_list/spore_prediction/presence_absence.csv',
        pca='plots/spor_presence_absence/pca.pdf',
        gene_dist='plots/spor_presence_absence/gene_dist.pdf'
    script:
        'scripts/spor_presence_absence_matrix.py'

# we need to get the genomes weller and wu used which will be our test set
# we can pull these from one of their the alignment fasta's
# This might all be irrelevant actually, I can download using taxonomy id
# through the ncbi datasets commandline tool
rule get_WW_ascension_nums:
    input:
        afa='data/WW_prot_alignment.aln',
        xlsx='data/weller_and_wu_sporulators.xlsx'
    output:
        asc='data/WW_ascensions.txt',
        taxa='data/WW_ascension_taxa.csv'
    script:
        'scripts/pull_WW_ascension_numbers.py'

rule download_WW_genomes:
    input:
        # csv='data/missing_ww.csv'
        csv='data/WW_ascension_taxa.csv'
    output:
        prot_dir=directory('data/weller_wu_labelled/proteins/'),
        gen_dir=directory('data/weller_wu_labelled/genomes/'),
        annot_dir=directory('data/weller_wu_labelled/annotations/')
    script:
        'scripts/download_WW_genomes.py'
    
# Now that I have the WW genomes I can get a precence absence matrix for them
rule ww_protein_blast:
    input:
        ref='data/bsub_ref_spor_genes/ref_gene_list.txt',
        ref_dir='data/bsub_ref_spor_genes/',
        prot="data/weller_wu_labelled/ww_genomes.txt",
        prot_dir="data/weller_wu_labelled/proteins/"
    output:
        db=directory('data/weller_wu_labelled/pblastdbs/'),
        out=directory('data/weller_wu_labelled/ref_spor_pblast_out/')
    script:
        'scripts/run_spor_ref_pblast.py'

rule ww_presence_absence_matrix:
    input:
        blast_dir='data/weller_wu_labelled/ref_spor_pblast_out/',
        header='data/Jordan_MAG_list/blastout_header_fmt_6.txt'
    output:
        mat='data/weller_wu_labelled/spore_prediction/presence_absence.csv',
        pca='plots/spor_presence_absence/ww_pca.pdf',
        gene_dist='plots/spor_presence_absence/ww_gene_dist.pdf'
    script:
        'scripts/spor_presence_absence_matrix.py'

# lets actually explore this stuff now
rule pres_abs_umap:
    input:
        ww='data/weller_wu_labelled/spore_prediction/presence_absence.csv',
        jor='data/Jordan_MAG_list/spore_prediction/presence_absence.csv',
        labs='data/WW_ascension_taxa.csv'
    output:
        umap_plot='plots/umap/umap_labelled.pdf',
        ww_umap='data/weller_wu_labelled/spore_prediction/ww_pres_abs_umap.csv',
        jor_umap='data/Jordan_MAG_list/spore_prediction/jor_pres_abs_umap.csv'
    threads: 4
    script:
        'scripts/pres_abs_umap.py'

rule weighted_umap:
    input:
        ww='data/weller_wu_labelled/spore_prediction/presence_absence.csv',
        jor='data/Jordan_MAG_list/spore_prediction/presence_absence.csv',
        labs='data/WW_ascension_taxa.csv',
        genes='data/delta6-network-genes.csv'
    output:
        umap_betw='plots/umap/betweenness_umap.pdf',
        umap_pg='plots/umap/pagerank_umap.pdf',
        umap_deg='plots/umap/degree_umap.pdf'
    script:
        'scripts/weighted_umap.py'

rule umap_kmeans_2:
    input:
        ww='data/weller_wu_labelled/spore_prediction/ww_pres_abs_umap.csv',
        jor='data/Jordan_MAG_list/spore_prediction/jor_pres_abs_umap.csv'
    output:
        descision='plots/decision/kmeans.pdf',
        classes='data/Jordan_MAG_list/spore_prediction/jor_predicted_classes_2.csv'
    script:
        'scripts/umap_k_means_2.py'

rule umap_kmeans_4:
    input:
        ww='data/weller_wu_labelled/spore_prediction/ww_pres_abs_umap.csv',
        jor='data/Jordan_MAG_list/spore_prediction/jor_pres_abs_umap.csv'
    output:
        descision='plots/decision/kmeans4.pdf',
        classes='data/Jordan_MAG_list/spore_prediction/jor_predicted_classes_4.csv'
    script:
        'scripts/umap_k_means_4.py'

rule display_pres_abs_ww:
    input:
        ww='data/weller_wu_labelled/spore_prediction/presence_absence.csv',
        labs='data/WW_ascension_taxa.csv'
    output:
        mat='plots/presence_absence/full.pdf'
    script:
        'scripts/display_presence_absence.py'

rule random_forest_vars:
    input:
        classes = 'data/Jordan_MAG_list/spore_prediction/jor_predicted_classes_2.csv',
        mat = 'data/Jordan_MAG_list/spore_prediction/presence_absence.csv',
        genes = 'data/delta6-network-genes.csv'
    output:
        corr_mat='plots/random_forest/corr_mat.png',
        reg_degree='plots/random_forest/reg_outdegree.png',
        reg_betw='plots/random_forest/reg_betw.png',
        reg_pagerank='plots/random_forest/reg_pagerank.png',
        importance='plots/random_forest/importance.png',
        dec_pca='plots/decision/pca_space_umap_labels.png'
    threads: 4
    script:
        'scripts/random_forest_plots.py'

# I'll try to reproduce some of jordans results
rule copy_jordan_figs:
    input:
        mat = 'data/Jordan_MAG_list/spore_prediction/presence_absence.csv'
    output:
        'plots/count_only/number_spor_genes_dist.png'
    script:
        'scripts/jordan_replication.py'

# lets try and draw a tree with python
rule draw_a_tree:
    input:
        tree = 'data/Jordan_MAG_list/tree_files/concat_ribosomal.tree',
        jor='data/Jordan_MAG_list/spore_prediction/jor_predicted_classes_2.csv'
    output:
        tree_viz='plots/trees/first_tree.pdf'
    script:
        'scripts/draw_tree.py'

