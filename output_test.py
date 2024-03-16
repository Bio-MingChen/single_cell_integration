import os
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc

from memory_profiler import profile
# ecc_n_7 = sc.read('./data/ECC_N_7_QC.h5ad')
# ecc_n_8 = sc.read('./data/ECC_N_8_QC.h5ad')

# print(ecc_n_7)
# print(ecc_n_7.var_names)
# print(len(ecc_n_7.var_names))
# print(len(set(ecc_n_7.var_names)))

def output_data(concat_h5ad,prefix,cluster,outdir):
    """
    output umap / cluster / diff data
    outdir/data : umap /cluster
    outdir/Diff : diff data
    """
    data_path = os.path.join(outdir,'data')
    diff_path = os.path.join(outdir,'DIFF')
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    if not os.path.exists(diff_path):
        os.makedirs(diff_path)
    umap_path = os.path.join(data_path,prefix + '_UMAP.csv')
    cluster_path = os.path.join(data_path,prefix + '_cluster.csv')
    # data
    df = pd.DataFrame(concat_h5ad.obsm['X_umap'],columns=["UMAP_1","UMAP_2"],index=concat_h5ad.obs_names)
    df.to_csv(umap_path,index=True,index_label='Barcode')
    cluster_df = pd.DataFrame({'Cluster':concat_h5ad.obs.loc[:,cluster]})
    cluster_df.to_csv(cluster_path,index=True,index_label='Barcode')
    # diff
    rank_gene = concat_h5ad.uns['rank_genes_groups']
    cluster_number = len(rank_gene['names'][0])
    rank_gene_dict = defaultdict(dict)
    for idx in range(cluster_number):
        for j in ['names','logfoldchanges','pvals','pvals_adj']:
            rank_gene_dict[str(idx)][j] = [i[idx] for i in rank_gene[j]]
    for key in rank_gene_dict:
        gene_order = rank_gene_dict[key]['names']
        rank_gene_dict[key]['Gene'] = concat_h5ad.var.loc[gene_order,'gene_ids']
        rank_gene_dict[key]['pct.1'] = rank_gene['pts'].loc[gene_order,key].tolist()
        rank_gene_dict[key]['pct.2'] = rank_gene['pts_rest'].loc[gene_order,key].tolist()

        diff_df = pd.DataFrame({
            'Gene':rank_gene_dict[key]['Gene'],
            'Gene_name':rank_gene_dict[key]['names'],
            'cluster': key,
            'avg_logFC':rank_gene_dict[key]['logfoldchanges'],
            'p_val':rank_gene_dict[key]['pvals'],
            'p_val_adj':rank_gene_dict[key]['pvals_adj'],
            'pct.1':rank_gene_dict[key]['pct.1'],
            'pct.2':rank_gene_dict[key]['pct.2'],
        })
        up_df = diff_df[diff_df['p_val_adj']<0.05 & diff_df['avg_logFC'] > 0]
        down_df = diff_df[diff_df['p_val_adj']<0.05 & diff_df['avg_logFC'] < 0]
        up_df.sort_values(by='avg_logFC',ascending=False)
        top50_df = up_df.iloc[:50,:]
        
        up_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_up.xls')
        down_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_down.xls')
        top50_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_top50.xls')
        up_df.to_csv(up_outpath,sep='\t',index=False)
        down_df.to_csv(down_outpath,sep='\t',index=False)
        top50_df.to_csv(top50_outpath,sep='\t',index=False)

ecc = sc.read('./data/ECC.h5ad')
print(ecc)
print(ecc.obs)
print(ecc.var)
print(ecc.obsm['X_umap'])
df = pd.DataFrame(ecc.obsm['X_umap'],columns=["UMAP_1","UMAP_2"],index=ecc.obs_names)
df.to_csv('ECC_UMAP.csv',index=True,index_label='Barcode')
cluster_df = pd.DataFrame({'Cluster':ecc.obs.loc[:,'louvain']})
print(cluster_df)
cluster_df.to_csv('ECC_cluster.csv',index=True,index_label='Barcode')
output_data(ecc,'ECC','louvain','.')
# df['Barcode'] = ecc.obs_names
# reorder_df = df.iloc[:,[2,0,1]]
# reorder_df.to_csv('ECC_UMAP.csv')

