import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
# do not show figure in window Linux mode
import matplotlib as mpl
mpl.use('Agg')

celltype_markers = {}

h5ad_file = sys.argv[1]
celltype_file = sys.argv[2]
outfigure = sys.argv[3]
with open(celltype_file,'r') as indata:
    for line in indata:
        celltype,markers = line.strip().split('\t')
        if celltype not in celltype_markers:
            celltype_markers[celltype] = markers.split(',')

if not os.path.exists(outfigure):
    os.makedirs(outfigure)
h5ad = sc.read(h5ad_file)
sc.settings.set_figure_params(dpi=100, 
        frameon=True, 
        figsize=(5, 5), 
        facecolor='white',
        )
# sc.settings.autoshow = False
sc.settings.figdir = outfigure

for celltype in celltype_markers:
    show_gene_list = []
    missed_gene_list = []
    for gene in celltype_markers[celltype]:
        if gene in h5ad.var_names:
            show_gene_list.append(gene)
        else:
            missed_gene_list.append(gene)
    if missed_gene_list:
        print('{celltype} have missed gene: {missed_gene_list}'.format(
            celltype=celltype,missed_gene_list = missed_gene_list
        ))
    sc.pl.umap(h5ad,color=["louvain"] + show_gene_list,
        color_map=mpl.cm.coolwarm,save='_{celltype}.pdf'.format(celltype=celltype.replace(" ","_")))
    sc.pl.stacked_violin(h5ad,show_gene_list,groupby='louvain',cmap='Reds',
        save='_{celltype}.violin.pdf'.format(celltype=celltype.replace(" ","_")))
