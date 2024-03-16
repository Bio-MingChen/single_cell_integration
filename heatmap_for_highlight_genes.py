import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
# do not show figure in window Linux mode
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
    "11-NK":["TYROBP","XCL1","NCAM1","KLRD1"], #NCAM1 (CD56)
    "12-B cell":["HLA-DRA","CD79B","CD79A"],
    "13 Proliferating cells":["TYMS","MKI67","KIAA0101","TOP2A","ZWINT"] }

cluster_names = ["CD8_C0_CD99","B_C1_BANK1","CD8_C2_FKBP4","CD4+ T-Naïve","CD4+ T-EM Cell",
    "Lym_Epi","CD8_C6_ETV1","CD4-Treg","Plasma B","NK_C9_FCGR3A","CD4+ T-H","NK_C11_CD38",
    "B_C12_IL3RA","Proliferating cells"]
marker_genes = set([j for k,v in celltype_markers.items() for j in v])
marker_genes = list(marker_genes)
print(marker_genes)
h5ad_file = sys.argv[1]
outfigure = sys.argv[2]
h5ad = sc.read(h5ad_file)
# sample
print(h5ad.obs.groupby('louvain').count())
sample_index = h5ad.obs.groupby('louvain').sample(6000).index
# add color group
color_number = len(np.unique(h5ad.obs['louvain']))
colors = sc.pl.palettes.vega_20[:color_number]
colors_dict = dict(zip(np.unique(h5ad.obs['louvain']),colors))
row_colors = h5ad.obs.loc[sample_index,:]['louvain'].map(colors_dict).to_list()
#prep data
heatmap_data = h5ad[sample_index,marker_genes].to_df()
#plot heatmap
plt.figure()
sns.clustermap(heatmap_data,row_cluster=True,col_cluster=False,cmap="vlag",standard_scale=1,
	yticklabels=False,row_colors=row_colors)
plt.savefig(outfigure)
