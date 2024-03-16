import sys
import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad

if len(sys.argv) != 3:
    print("""
    Usage: 
        python genrate_3D_umap_data.py h5ad ofile 
    """)
h5ad = sys.argv[1]
ofile = sys.argv[2]

h5ad_data = sc.read(h5ad)

#Re-run umap with 3-dimensions
sc.tl.umap(h5ad_data,n_components=3)
umap_df = pd.DataFrame(h5ad_data.obsm['X_umap'],index=h5ad_data.obs_names,
    columns=['UMAP_1','UMAP_2','UMAP_3'])

if 'louvain' in h5ad_data.obs:
    umap_df['cluster'] = h5ad_data.obs['louvain']
elif 'leiden' in h5ad_data.obs:
    umap_df['cluster'] = h5ad_data.obs['leiden']

umap_df.to_csv(ofile,sep='\t',index=True,index_label='label')