#!/usr/bin/env python
# -*- coding=utf-8 -*-

#This script aim to integrate big number cells with bbknn/ingest/harmony integration algorithm
import os
from uuid import uuid1
import argparse
from collections import OrderedDict,defaultdict
from subprocess import call

import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from memory_profiler import profile
# do not show figure in window Linux mode
import matplotlib
matplotlib.use('Agg')

class TitleParser():
    """
    Reading in a file title and than getting indicated column element.
    Title names are all CaseIgnore
    """
    
    def __init__(self,title):
        self.title_list = [i.lower() for i in title.strip().split('\t')]

    def get_field(self,line_list,colname,check=True):
        """
        Reading in a list and returning the element with the 
        index which colname in title's list
        """
        if check:
            if len(self.title_list) != len(line_list):
                raise Exception("Title length differs with line!")
        try:
            idx = self.title_list.index(str(colname).lower())
        except:
            print('{colname} not in title！'.format(colname=colname))
            return None

        return line_list[idx]

    def have_title(self,colname):
        """
        Judging whether a colname is in title
        """
        if str(colname).lower() in self.title_list:
            return True
        return False
    
    def get_idx(self,colname):
        """
        return index of columns by their name in title list
        """
        idx = self.title_list.index(str(colname).lower())

        return idx

def parse_config(config,outdir):
    """
    config format
    sample\th5ad_path\trds_path
    SampleA\t/path/to/SampleA.h5ad
    SampleB\t/path/to/SampleB.h5ad
    if rds file,transfer it to .h5ad by sceasy
    """
    h5ad_dir = os.path.join(outdir,'h5ads')
    if not os.path.exists(h5ad_dir):
        os.makedirs(h5ad_dir)
    h5ad_dict = OrderedDict()
    with open(config,'r') as indata:
        title = indata.readline()
        tp = TitleParser(title)
        title_list = title.strip().split('\t')
        if 'h5ad_path' in title_list:
            for line in indata:
                line_list = line.strip().split('\t')
                sample = tp.get_field(line_list,'sample')
                h5ad = tp.get_field(line_list,'h5ad_path')
                h5ad_dict[sample] = sc.read(h5ad)
        elif 'rds_path' in title_list:
            config_transfer = config + '.' + str(uuid1())
            raw_h5ad_dir = os.path.join(h5ad_dir,'raw_h5ads')
            if not os.path.exists(raw_h5ad_dir):
                os.makedirs(raw_h5ad_dir)
            with open(config_transfer,'w') as odata:
                title = '\t'.join(title_list + ['h5ad_path']) + '\n'
                odata.write(title)
                for line in indata:
                    line_list = line.strip().split('\t')
                    sample = tp.get_field(line_list,'sample')
                    rds_path = tp.get_field(line_list,'rds_path')
                    h5ad_file = os.path.splitext(os.path.basename(rds_path))[0] + '.h5ad'
                    h5ad_path = os.path.join(raw_h5ad_dir,h5ad_file)
                    odata.write("\t".join(line_list + [h5ad_path]) + '\n')
            h5ad_dict = transfer_rds(config_transfer,outdir)
    return h5ad_dict

def transfer_rds(config_transfer,outdir):
    """
    run rds to h5ad transfer
    """
    cmd = """
    source /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/bin/activate /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/sceasy_env 
    Rscript /TJPROJ6/SC/personal_dir/chenming/research/scanpy_pipeline/scripts/sc_transfer.R \
        -c {config_transfer}
    """.format(config_transfer=config_transfer)
    assert not call(cmd,shell=True)
    h5ad_dict = parse_config(config_transfer,outdir)
    return h5ad_dict

def get_palette(h5ad,column):
    color_number = len(set(h5ad.obs[column]))
    if color_number <= 20:
        return sc.pl.palettes.vega_20_scanpy
    # elif (color_number > 20) and (color_number <= 28):
    #     return sc.pl.palettes.zeileis_28
    else:
        return sc.pl.palettes.godsnot_102

def remove_batch_effect(concat_h5ad,method,cluster,run_tsne,resolution):
    """
    method : bbknn / ingest / harmony
    cluster : louvain / leiden
    """
    if method == 'bbknn':
        sc.pp.highly_variable_genes(concat_h5ad,n_top_genes=2000)
        sc.tl.pca(concat_h5ad)
        sc.external.pp.bbknn(concat_h5ad,batch_key='batch')
        sc.tl.umap(concat_h5ad)
        if cluster == 'louvain':
            sc.tl.louvain(concat_h5ad,resolution=resolution)
        elif cluster == 'leiden':
            sc.tl.leiden(concat_h5ad,resolution=resolution)
        batch_palette = get_palette(concat_h5ad,'batch')
        cluster_palette = get_palette(concat_h5ad,cluster)
        if run_tsne:
            sc.tl.tsne(concat_h5ad)
            sc.pl.tsne(concat_h5ad,color='batch',save='_batch.png',palette=batch_palette)
            sc.pl.tsne(concat_h5ad,color=cluster,save='_{cluster}.png'.format(cluster=cluster),palette=cluster_palette)
        sc.pl.umap(concat_h5ad,color='batch',save='_batch.png',palette=batch_palette)
        sc.pl.umap(concat_h5ad,color=cluster,save='_{cluster}.png'.format(cluster=cluster),palette=cluster_palette)

    elif method == 'ingest':
        batches = concat_h5ad.obs.batch.cat.categories
        print('batches is {batches}'.format(batches=batches))
        concat_ref = concat_h5ad[concat_h5ad.obs.batch == batches[0]]
        sc.pp.highly_variable_genes(concat_ref,n_top_genes=2000)
        sc.pp.pca(concat_ref)
        sc.pp.neighbors(concat_ref)
        if cluster == 'louvain':
            sc.tl.louvain(concat_ref,resolution=resolution)
        elif cluster == 'leiden':
            sc.tl.leiden(concat_ref,resolution=resolution)
        sc.tl.umap(concat_ref)
        sc.pl.umap(concat_ref, color=cluster,save='_ingest_ref.png')
        adatas = [concat_h5ad[concat_h5ad.obs.batch == i].copy() for i in batches[1:]]
        sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
        for idx,adata in enumerate(adatas):
            print('...integrating batch {idx}'.format(idx=idx+1))
            sc.tl.ingest(adata,concat_ref,obs=cluster)
            adata.uns['louvain_colors'] = concat_ref.uns['louvain_colors']
            # sc.pl.umap(adata,color=cluster,save='_batch_{idx}.png'.format(idx=idx+1))
        concat_h5ad = concat_ref.concatenate(adatas,batch_categories=batches)
        concat_h5ad.obs[cluster] = concat_h5ad.obs[cluster].astype('category')
        concat_h5ad.obs[cluster].cat.reorder_categories(concat_ref.obs[cluster].cat.categories,inplace=True)
        concat_h5ad.uns[cluster +'_colors'] = concat_ref.uns[cluster + '_colors']
        batch_palette = get_palette(concat_h5ad,'batch')
        cluster_palette = get_palette(concat_h5ad,cluster)
        sc.pl.umap(concat_h5ad,color='batch',save='_batch.png',palette=batch_palette)
        sc.pl.umap(concat_h5ad,color=cluster,save='_{cluster}.png'.format(cluster=cluster),palette=cluster_palette)
    elif method == 'harmony':
        sc.pp.highly_variable_genes(concat_h5ad,n_top_genes=2000)
        sc.pp.pca(concat_h5ad)
        sc.external.pp.harmony_integrate(concat_h5ad,'batch')
        sc.pp.neighbors(concat_h5ad,use_rep='X_pca_harmony')
        sc.tl.umap(concat_h5ad)
        if cluster == 'louvain':
            sc.tl.louvain(concat_h5ad,resolution=resolution)
        elif cluster == 'leiden':
            sc.tl.leiden(concat_h5ad,resolution=resolution)
        batch_palette = get_palette(concat_h5ad,'batch')
        cluster_palette = get_palette(concat_h5ad,cluster)
        sc.pl.umap(concat_h5ad,color='batch',save='_batch.png',palette=batch_palette)
        sc.pl.umap(concat_h5ad,color=cluster,save='_{cluster}.png'.format(cluster=cluster),palette=cluster_palette)
        if run_tsne:
            sc.tl.tsne(concat_h5ad)
            sc.pl.tsne(concat_h5ad,color='batch',save='_batch.png',palette=batch_palette)
            sc.pl.tsne(concat_h5ad,color=cluster,save='_{cluster}.png'.format(cluster=cluster),palette=cluster_palette)
    return concat_h5ad

def output_data(concat_h5ad,prefix,cluster,outdir):
    """
    output highly variable genes
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
    # highly variable genes
    highly_variable_genes = concat_h5ad.var_names[concat_h5ad.var['highly_variable']].to_list()
    highly_variable_ofile = os.path.join(data_path,prefix + '_highly_variable.txt')
    with open(highly_variable_ofile,'w') as odata:
        for gene in highly_variable_genes:
            odata.write(gene + '\n')
    # tsne
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
    all_top100_dff_df = None
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
        up_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_up.xls')
        top100_outpath = os.path.join(diff_path,'Cluster_' + key + '_diff_top100.xls')
        up_df.to_csv(up_outpath,sep='\t',index=False)
        top100_df.to_csv(top100_outpath,sep='\t',index=False)
    all_diff_outpath = os.path.join(diff_path,'Cluster_diff.xls')
    all_up_diff_outpath = os.path.join(diff_path,'Cluster_up_diff.xls')
    all_top100_diff_outpath = os.path.join(diff_path,'Cluster_top100_diff.xls')
    all_up_diff_df.to_csv(all_up_diff_outpath,sep='\t',index=False)
    all_top100_diff_df.to_csv(all_top100_diff_outpath,sep='\t',index=False)
    all_diff_df.to_csv(all_diff_outpath,sep='\t',index=False)

@profile  
def integrate_pipeline(h5ad_dict,method,cluster,outdir,prefix,diff_method,gene_ids_df,run_tsne,resolution,barcodes):
    """
    """
    # make sure all the data have the same variables(features/genes/transcripts)
    var_names = ''
    for _,h5ad in h5ad_dict.items():
        if not isinstance(var_names,pd.core.indexes.base.Index):
            var_names = h5ad.var_names
        else:
            var_names = var_names & h5ad.var_names
    var_names = var_names & gene_ids_df.index
    #save var_names
    for k,h5ad in h5ad_dict.items():
        h5ad_dict[k] = h5ad[:,var_names]
    #rename the observations by index
    for idx,(_,h5ad) in enumerate(h5ad_dict.items()):
        h5ad.obs_names = np.array([i.replace('-1','-{idx}'.format(idx=idx+1)) for i in h5ad.obs_names])
    #concatenation
    concat_h5ad = ad.concat(h5ad_dict.values(),label='batch',keys=h5ad_dict.keys())
    # subset barcodes
    if barcodes:
        concat_h5ad = concat_h5ad[barcodes,:]
    # remove batch effect with bbknn/ingest/harmony
    concat_h5ad = remove_batch_effect(concat_h5ad,method,cluster,run_tsne,resolution)
    # add gene_ids
    concat_h5ad.var['gene_ids'] = gene_ids_df.loc[concat_h5ad.var_names,'gene_ids']
    # diff
    sc.tl.rank_genes_groups(concat_h5ad,cluster,method=diff_method,pts=True)
    # output data
    output_data(concat_h5ad,prefix,cluster,outdir)
    # save h5ad
    ofile = os.path.join(outdir,'h5ads',prefix + '.h5ad')
    concat_h5ad.write(ofile)

def read_gene_ids(gene_ids_file):
    gene_ids_dict = {'gene_ids':[],'genenames':[]}
    uniq_gene_set = set()
    with open(gene_ids_file,'r') as indata:
        for line in indata:
            gene_id,genename = line.strip().split('\t')
            # remove repeat genename
            if genename not in uniq_gene_set:
                uniq_gene_set.add(genename)
                gene_ids_dict['gene_ids'].append(gene_id)
                gene_ids_dict['genenames'].append(genename)
    gene_ids_df = pd.DataFrame(gene_ids_dict,index=gene_ids_dict['genenames'])
    return gene_ids_df

def main(args):

    outdir = args.get('outdir')
    outfigure = os.path.join(outdir,'figures')
    # parse config
    config = args.get('config')
    h5ad_dict = parse_config(config,outdir)
    print('loading h5ad completion!')
    # exit(0)
    # global settings
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=100, 
        frameon=True, 
        figsize=(5, 5), 
        facecolor='white',
        )
    sc.settings.autoshow = False
    # sc.settings.autosave = True
    #figdir
    #Directory for saving figures (default './figures/').
    #cachedir
    #Directory for cache files (default './cache/').
    #datasetdir
    #Directory for example datasets (default './data/').
    sc.settings.figdir = outfigure
    # integration
    method = args.get('method')
    cluster = args.get('cluster')
    prefix = args.get('prefix')
    outdir = args.get('outdir')
    diff_method = args.get('diff_method')
    gene_ids = args.get('gene_ids')
    gene_ids_df = read_gene_ids(gene_ids)
    run_tsne = args.get('use_tsne')
    resolution = args.get('resolution')
    barcodes_file = args.get('barcodes')
    barcodes = []
    if barcodes_file:
        with open(barcodes_file,'r') as indata:
            for line in indata:
                if line.strip() == 'Barcode':
                    continue
                barcodes.append(line.strip())
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    integrate_pipeline(h5ad_dict,method,cluster,outdir,prefix,diff_method,gene_ids_df,run_tsne,resolution,barcodes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--config','-c',help='config file with samplenames and their h5ad/rds files')
    parser.add_argument('--prefix','-p',help='Integration prefix to use')
    parser.add_argument('--gene_ids','-g',default='/TJPROJ5/SC/pipeline/scRNA/pipeline1.0/lib/05.Seurat/GRCh38_1.2.0_genes.tsv',help='file to save the Ensembl ID of genes [default:%(default)s]')
    parser.add_argument('--outdir','-o',help='Output directory')
    parser.add_argument('--method','-m',default='bbknn',choices=['bbknn','ingest','harmony'],
        help='method for batch effect remove choose from %(choices)s [defalut:%(default)s]')
    parser.add_argument('--cluster',default='louvain',choices=['louvain','leiden'],
        help='cluster method to use choose from %(choices)s [defalut:%(default)s]')
    parser.add_argument('--diff_method','-d',default="wilcoxon",
        choices=["logreg","t-test","wilcoxon","t-test_overestim_var"],
        help='method to rank genes choose from %(choices)s [default:%(default)s]')
    parser.add_argument('--use_tsne','-u',action='store_true',help='run tsne at the same time,only useful to bbknn and harmony')
    parser.add_argument('--resolution','-r',type=float,default=1,help='resolution for cluster method louvain and leiden [default:%(default)s]')
    parser.add_argument('--barcodes','-b',help='a file contains barcodes to analyse')
    args = vars(parser.parse_args())
    main(args)