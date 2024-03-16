import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
# do not show figure in window Linux mode
import matplotlib
matplotlib.use('Agg')

celltype_markers = {
    "0-CD8":["GZMK","CD8B","CD8A"],
    "1-B cell":["CD79A","HLA-DRA","BCL11A","CD79B"],
    "2 CD8+ T-EM Cell":["DUSP2","CD8A","GZMK"],
    "3-CD4+ T-Naïve":["LEF1","IL7R","CD4","SELL","TCF7"],
    "4-CD4+ T-EM Cell":["IL7R","KLRB1","CD40LG","KLRB1","CCR6","RORC"],
    "5 Lymphoid_Epithelial?":["TFF3","CXCL2","LYZ","KRT18","TAGLN","ANKAR","CXCL3","KRT19",],
    "6-CD8":["CD8A","CD8B"],
    "7-CD4-Treg":["IL1R2","TNFRSF4","IL2RA","TNFRSF18","CD4","IL7R"],
    "8-Plasma B":["JCHAIN"],
    "9-NK":["GNLY","NKG7","TYROBP","KLRD1"],
    "10-CD4+ T-H":["IL7R","KLRB1","CAPG","KLRB1","CCR6","RORC"],
    "11-NK":["TYROBP","XCL1","NCAM1","CD56","KLRD1"], #NCAM1 (CD56)
    "12-B cell":["HLA-DRA","CD79B","CD79A"],
    "13 Proliferating cells":["TYMS","MKI67","KIAA0101","TOP2A","ZWINT"] }

h5ad_file = sys.argv[1]
outfigure = sys.argv[2]
if not os.path.exists(outfigure):
    os.makedirs(outfigure)
h5ad_file = sc.read(h5ad_file)
sc.settings.set_figure_params(dpi=100, 
        frameon=True, 
        figsize=(5, 5), 
        facecolor='white',
        )
sc.settings.autoshow = False
sc.settings.figdir = outfigure
for celltype in celltype_markers:
    show_gene_list = []
    missed_gene_list = []
    for gene in celltype_markers[celltype]:
        if gene in h5ad_file.var_names:
            show_gene_list.append(gene)
        else:
            missed_gene_list.append(gene)
    if missed_gene_list:
        print('{celltype} have missed gene: {missed_gene_list}'.format(
            celltype=celltype,missed_gene_list = missed_gene_list
        ))
    sc.pl.umap(h5ad_file,color=show_gene_list,save='{celltype}.pdf'.format(celltype=celltype))
