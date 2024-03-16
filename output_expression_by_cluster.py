import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad


cancers = sc.read('./data/pbmc68k.h5ad')
cancers_df = cancers.to_df()
cancers_df['cluster'] = cancers.obs['louvain']
clusters = np.unique(cancers.obs['louvain'])
for cluster in clusters:
    cluster_df = cancers_df[cancers_df['cluster'] == cluster]
    cluster_df.drop(labels=['cluster'],axis=1,inplace=True)
    outfile = 'cancer_cluster{}.tsv'.format(cluster)
    cluster_df.T.to_csv(outfile,sep='\t',index=True,index_label='gene')