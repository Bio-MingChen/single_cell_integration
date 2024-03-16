#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import argparse
from subprocess import call
from collections import defaultdict

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from memory_profiler import profile

import matplotlib
# do not show figure in window Linux mode
matplotlib.use('Agg')

def get_palette(h5ad,column):
    color_number = len(set(h5ad.obs[column]))
    if color_number <= 20:
        return sc.pl.palettes.vega_20_scanpy
    # elif (color_number > 20) and (color_number <= 28):
    #     return sc.pl.palettes.zeileis_28
    else:
        return sc.pl.palettes.godsnot_102

def output_data(concat_h5ad,prefix,cluster,outdir,r):
    """
    output umap / cluster / diff data
    outdir/data : umap /cluster
    outdir/Diff : diff data
    """
    diff_path = os.path.join(outdir,'DIFF_{}'.format(r))
    data_path = os.path.join(outdir,'data_{}'.format(r))
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    if not os.path.exists(diff_path):
        os.makedirs(diff_path)
    umap_path = os.path.join(data_path,prefix + '_{}_UMAP.csv'.format(r))
    cluster_path = os.path.join(data_path,prefix + '_{}_cluster.csv'.format(r))
    # data
    df = pd.DataFrame(concat_h5ad.obsm['X_umap'],columns=["UMAP_1","UMAP_2"],index=concat_h5ad.obs_names)
    df.to_csv(umap_path,index=True,index_label='Barcode')
    cluster_df = pd.DataFrame({'Cluster':concat_h5ad.obs.loc[:,cluster]})
    cluster_df.to_csv(cluster_path,index=True,index_label='Barcode')
    # highly variable genes
    highly_variable_genes = concat_h5ad.var_names[concat_h5ad.var['highly_variable']].to_list()
    highly_variable_ofile = os.path.join(data_path,prefix + '_highly_variable.txt')
    with open(highly_variable_ofile,'w') as odata:
        for gene in highly_variable_genes:
            odata.write(gene + '\n')
    # all genes and their gene ids
    gene_ids_ofile = os.path.join(data_path,prefix + '_gene_ids.txt')
    print(concat_h5ad.var)
    concat_h5ad.var['gene_ids'].to_csv(gene_ids_ofile,index=True,index_label='Gene')
    # tsne
    print(concat_h5ad)
    if 'X_tsne' in concat_h5ad.obsm:
        df = pd.DataFrame(concat_h5ad.obsm['X_tsne'],columns=["tSNE_1","tSNE_2"],index=concat_h5ad.obs_names)
        tsne_path = os.path.join(data_path,prefix + '_TSNE.csv')
        df.to_csv(tsne_path,index=True,index_label='Barcode')
    # diff
    rank_gene = concat_h5ad.uns['rank_genes_groups']
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
        rank_gene_dict[key]['pct.1'] = rank_gene['pts'].loc[gene_order,key].tolist()
        rank_gene_dict[key]['pct.2'] = rank_gene['pts_rest'].loc[gene_order,key].tolist()

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
    all_diff_outpath = os.path.join(diff_path,'Cluster_diff.xls')
    all_up_diff_outpath = os.path.join(diff_path,'Cluster_up_diff.xls')
    all_top100_diff_outpath = os.path.join(diff_path,'Cluster_top100_diff.xls')
    all_up_diff_df.to_csv(all_up_diff_outpath,sep='\t',index=False)
    all_top100_diff_df.to_csv(all_top100_diff_outpath,sep='\t',index=False)
    all_diff_df.to_csv(all_diff_outpath,sep='\t',index=False)

def transfer_h5ad(h5ad,ofile):
    """
    transfer h5ad to rds 
    """
    cmd = """
    source /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/bin/activate /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/sceasy_env 
    Rscript /TJPROJ6/SC/personal_dir/chenming/research/scanpy_pipeline/scripts/sc_transfer.R \
        -i {h5ad} \
        -o {ofile} \
        -m scale.data
    """.format(h5ad=h5ad,ofile=ofile)

    assert not call(cmd,shell=True)
    print('Transfer {h5ad} to {ofile} successfully!'.format(h5ad=h5ad,ofile=ofile))

@profile
def run_cluster_and_diff(in_h5ad,r,method,odir,diff_method,prefix,save_data,use_raw,output_rds,batch_name):
    """
    Re run cluster and find different expression genes 
    """
    print('run cluster with {resolution} for method {method}'.format(
        resolution = r,method = method
    ))
    sc.settings.figdir = os.path.join(odir,'figures_with_resolution_{}'.format(r))
    if method == 'louvain':
        sc.tl.louvain(in_h5ad,resolution=r)
    elif method == 'leiden':    
        sc.tl.leiden(in_h5ad,resolution=r)
    print('cluster with {resolution}: {clusters}'.format(
        clusters = in_h5ad.obs[method].cat.categories,
        resolution = r,
    ))
    if batch_name in in_h5ad.obs:
        batch_palette = get_palette(in_h5ad,batch_name)
        sc.pl.umap(in_h5ad,color='batch',save='_{}_batch.png'.format(r),palette=batch_palette)
    
    cluster_palette = get_palette(in_h5ad,method)
    sc.pl.umap(in_h5ad,color=method,save='_{r}_{cluster}.png'.format(r=r,cluster=method),palette=cluster_palette)
    # diff
    print(f"use_raw is {use_raw}")
    sc.tl.rank_genes_groups(in_h5ad,method,method=diff_method,use_raw=use_raw,pts=True)
    # save diff
    output_data(in_h5ad,prefix,method,odir,r)
    # save data
    if save_data:
        h5ad_dir = os.path.join(odir,'h5ads')
        if not os.path.exists(h5ad_dir):
            os.makedirs(h5ad_dir)
        ofile = os.path.join(odir,'h5ads',prefix + '_resolution{r}.h5ad'.format(r=r))
        in_h5ad.write(ofile, compression = "gzip")
        if output_rds:
            rds_ofile = ofile.replace('.h5ad','.rds')
            transfer_h5ad(ofile,rds_ofile)

def main(args):
    infile = args.get('infile')
    in_h5ad = sc.read(infile)
    method = args.get('method')
    odir = args.get('odir')
    if not os.path.exists(odir):
        os.makedirs(odir)
    diff_method = args.get('diff_method')
    prefix = args.get('prefix')
    save_data = args.get('save_data')
    use_raw = args.get('use_raw')
    output_rds =args.get('output_rds')
    resolutions = [float(i) for i in args.get('resolution').split(',')]
    batch_name = args.get('batch_name')
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=100, 
        frameon=True, 
        figsize=(5, 5), 
        facecolor='white',
        )
    for r in resolutions:
        run_cluster_and_diff(in_h5ad,r,method,odir,diff_method,prefix,save_data,use_raw,output_rds,batch_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''
        This program is used to find the optimized clusters by indicating
        different resolution for louvain or leiden method
        the Seurat suggest range is 0.4~1.2
        '''
    )
    parser.add_argument('-i','--infile',help='h5ad file')
    parser.add_argument('-r','--resolution',help='one or more resolutions to use for example: 0.5 or 1,1.2')
    parser.add_argument('-m','--method',default='louvain',choices=['louvain','leiden'],
        help='method to cluster choose from %(choices)s')
    parser.add_argument('-o','--odir',default=os.environ['PWD'],help='output directory,default is current directory')
    parser.add_argument('--diff_method','-d',default="wilcoxon",
        choices=["logreg","t-test","wilcoxon","t-test_overestim_var"],
        help='method to rank genes choose from %(choices)s [default:%(default)s]')
    parser.add_argument('--prefix','-p',help='Integration prefix to use')
    parser.add_argument('--batch_name','-bn',default='batch',help='name of batch column, default:%(default)s')
    parser.add_argument('--save_data','-s',action='store_true',help='whether output h5ad for each resolution [default:%(default)s]')
    parser.add_argument('--output_rds','-or',action='store_true',help='whether transfer h5ad to rds if save_data is True [default:%(default)s]')
    parser.add_argument('--use_raw','-u',action='store_true',default=None,
        help='If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale, you can also get away without using .raw at all.')
    args = vars(parser.parse_args())
    main(args)