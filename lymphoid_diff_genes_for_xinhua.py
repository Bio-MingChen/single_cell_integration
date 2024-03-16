import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
import math
import seaborn as sns

# do not show figure in window Linux mode
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

h5ad_file = sys.argv[1]
sample_type_file = sys.argv[2]
outdir = sys.argv[3]

def volcano_plot(diff_df,output_path):
    """
    """
    diff_df['-logQ'] = [-math.log(i,10) if i else np.nan for i in diff_df['p_val_adj']]
    # fill na
    max_value = max(diff_df['-logQ'][~diff_df['-logQ'].isna()])
    diff_df['-logQ'] = diff_df['-logQ'].fillna(max_value)
    diff_df['significant'] = 'normal'
    diff_df['significant'][(diff_df['avg_logFC']>0) & (diff_df['-logQ']>=-math.log(0.05,10))] = 'up'
    diff_df['significant'][(diff_df['avg_logFC']<0) & (diff_df['-logQ']>=-math.log(0.05,10))] = 'down'
    plt.figure() 
    ax = sns.scatterplot(x="avg_logFC", y="-logQ",
                        hue='significant',
                        hue_order = ('up','normal','down'),
                        palette=("#E41A1C","grey","#377EB8",),
                        data=diff_df)
    ax.set_ylabel('-log(pvalue)',fontweight='bold')
    ax.set_xlabel('FoldChange',fontweight='bold')
    #fig = ax.get_figure()
    plt.savefig(output_path)

def save_diff(concat_h5ad,diff_path):
    if not os.path.exists(diff_path):
        os.makedirs(diff_path)
    rank_gene = concat_h5ad.uns['rank_genes_groups']
    print(rank_gene)
    cluster_number = len(rank_gene['names'][0])
    rank_gene_dict = defaultdict(dict)
    for idx in range(cluster_number):
        for j in ['names','logfoldchanges','pvals','pvals_adj']:
            rank_gene_dict[str(idx)][j] = [i[idx] for i in rank_gene[j]]
    all_diff_df = None
    all_up_diff_df = None
    all_top100_diff_df = None
    for key in rank_gene_dict:
        gene_order = rank_gene_dict[key]['names']
        rank_gene_dict[key]['Gene'] = concat_h5ad.var.loc[gene_order,'gene_ids']
        rank_gene_dict[key]['pct.1'] = rank_gene['pts'].loc[gene_order,'T'].tolist()
        rank_gene_dict[key]['pct.2'] = rank_gene['pts'].loc[gene_order,'N'].tolist()
        diff_df = pd.DataFrame({
            'Gene':rank_gene_dict[key]['Gene'],
            'Gene_name':rank_gene_dict[key]['names'],
            'cluster': key,
            'avg_logFC':np.array(rank_gene_dict[key]['logfoldchanges'],dtype='float64'),
            'p_val':np.array(rank_gene_dict[key]['pvals'],dtype='float64'),
            'p_val_adj':np.array(rank_gene_dict[key]['pvals_adj'],dtype='float64'),
            'pct.1':rank_gene_dict[key]['pct.1'],
            'pct.2':rank_gene_dict[key]['pct.2'],
        })
        if all_diff_df is None:
            all_diff_df = diff_df
        else:
            all_diff_df = all_diff_df.append(diff_df) # need to rename
        up_df = diff_df[(diff_df['p_val_adj']<0.05) & (diff_df['avg_logFC'] > 0) &(diff_df['pct.1'] > 0.1)]
        up_df.sort_values(by='avg_logFC',ascending=False,inplace=True)
        top100_df = up_df.iloc[:100,:]
        if all_up_diff_df is None:
            all_up_diff_df = up_df
        else:
            all_up_diff_df = all_up_diff_df.append(up_df) # need to rename
        if all_top100_diff_df is None:
            all_top100_diff_df = top100_df
        else:
            all_top100_diff_df = all_top100_diff_df.append(top100_df) # need to rename
        up_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_up.xls')
        top100_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_top100.xls')
        #up_df.to_csv(up_outpath,sep='\t',index=False)
        #top100_df.to_csv(top100_outpath,sep='\t',index=False)
    all_diff_outpath = os.path.join(diff_path,'Cluster_diff.xls')
    all_up_diff_outpath = os.path.join(diff_path,'Cluster_up_diff.xls')
    all_top100_diff_outpath = os.path.join(diff_path,'Cluster_top100_diff.xls')
    all_up_diff_df.to_csv(all_up_diff_outpath,sep='\t',index=False)
    all_top100_diff_df.to_csv(all_top100_diff_outpath,sep='\t',index=False)
    all_diff_df.to_csv(all_diff_outpath,sep='\t',index=False)
    return all_diff_df

h5ad = sc.read(h5ad_file)
sample_type_data = pd.read_csv(sample_type_file)
h5ad.obs['Barcode'] = h5ad.obs_names.to_list()
h5ad.obs = pd.merge(h5ad.obs,sample_type_data,on=['Barcode'],left_index=True)
h5ad.obs.index = h5ad.obs['Barcode']
print(h5ad.obs)
volcano_dir = os.path.join(outdir,'volcano_plot')
if not os.path.exists(volcano_dir):
    os.makedirs(volcano_dir)
diff_dir = os.path.join(outdir,'DIFF')
if not os.path.exists(diff_dir):
    os.makedirs(diff_dir)
clusters = h5ad.obs['louvain'].cat.categories
for cluster in clusters:
    cluster_index = h5ad.obs[h5ad.obs['louvain'] == cluster]['Barcode']
    print(cluster_index)
    subset_h5ad = h5ad[cluster_index]
    sc.tl.rank_genes_groups(subset_h5ad,groupby='sample_type',groups=['T'],
        reference='N',method='wilcoxon',pts=True)
    cluster_outdir = os.path.join(diff_dir,cluster)
    all_diff_df = save_diff(subset_h5ad,cluster_outdir)
    volcano_output_path = os.path.join(volcano_dir,'Cluster{}_volcano_plot.pdf'.format(cluster))
    volcano_plot(all_diff_df,volcano_output_path)
