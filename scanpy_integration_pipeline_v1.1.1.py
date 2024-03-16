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
import bbknn

# do not show figure in window Linux mode
import matplotlib
matplotlib.use('Agg')

class TitleParser():
    """
    Reading in a file title and than getting indicated column element.
    Title names are all CaseIgnore
    """
    
    def __init__(self,title,sep='\t'):
        self.title_list = [i.lower() for i in title.strip().split(sep)]

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
    source /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/bin/activate /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/transfer_env
    Rscript /TJPROJ6/SC/personal_dir/chenming/research/scanpy_pipeline/scripts/sc_transfer.R \
        -c {config_transfer}
    """.format(config_transfer=config_transfer)
    assert not call(cmd,shell=True)
    h5ad_dict = parse_config(config_transfer,outdir)
    return h5ad_dict

def transfer_h5ad(h5ad,ofile):
    """
    transfer h5ad to rds 
    """
    cmd = """
    source /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/bin/activate /TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/transfer_env
    Rscript /TJPROJ6/SC/personal_dir/chenming/research/scanpy_pipeline/scripts/sc_transfer.R \
        -i {h5ad} \
        -o {ofile} \
        -m scale.data
    """.format(h5ad=h5ad,ofile=ofile)

    assert not call(cmd,shell=True)
    print('Transfer {h5ad} to {ofile} successfully!'.format(h5ad=h5ad,ofile=ofile))

def get_palette(h5ad,column):
    color_number = len(set(h5ad.obs[column]))
    if color_number <= 20:
        return sc.pl.palettes.vega_20_scanpy
    # elif (color_number > 20) and (color_number <= 28):
    #     return sc.pl.palettes.zeileis_28
    else:
        return sc.pl.palettes.godsnot_102

def remove_batch_effect(concat_h5ad,method,cluster,run_tsne,resolution,n_comps,fontsize,batch_name):
    """
    method : bbknn / ingest / harmony
    cluster : louvain / leiden
    """
    if method == 'bbknn':
        sc.tl.pca(concat_h5ad,use_highly_variable=True,n_comps=n_comps)
        #sc.external.pp.bbknn(concat_h5ad,batch_key=batch_name)
        bbknn.bbknn(concat_h5ad,batch_key=batch_name)
        sc.tl.umap(concat_h5ad)
        if cluster == 'louvain':
            sc.tl.louvain(concat_h5ad,resolution=resolution)
        elif cluster == 'leiden':
            sc.tl.leiden(concat_h5ad,resolution=resolution)
        batch_palette = get_palette(concat_h5ad,batch_name)
        cluster_palette = get_palette(concat_h5ad,cluster)
        if run_tsne:
            sc.tl.tsne(concat_h5ad,n_pcs=n_comps,use_rep="X_pca")
            sc.pl.tsne(concat_h5ad,color=batch_name,save='_batch.png',palette=batch_palette,size=fontsize)
            sc.pl.tsne(concat_h5ad,
                color=cluster,
                save='_{cluster}.png'.format(cluster=cluster),
                palette=cluster_palette,
                size=fontsize)
        sc.pl.umap(concat_h5ad,color=batch_name,save='_batch.png',palette=batch_palette,size=fontsize)
        sc.pl.umap(concat_h5ad,
            color=cluster,
            save='_{cluster}.png'.format(cluster=cluster),
            palette=cluster_palette,
            size=fontsize)

    elif method == 'harmony':
        sc.pp.pca(concat_h5ad,use_highly_variable=True,n_comps=n_comps)
        sc.external.pp.harmony_integrate(concat_h5ad,batch_name)
        sc.pp.neighbors(concat_h5ad,use_rep='X_pca_harmony')
        sc.tl.umap(concat_h5ad)
        if cluster == 'louvain':
            sc.tl.louvain(concat_h5ad,resolution=resolution)
        elif cluster == 'leiden':
            sc.tl.leiden(concat_h5ad,resolution=resolution)
        batch_palette = get_palette(concat_h5ad,batch_name)
        cluster_palette = get_palette(concat_h5ad,cluster)
        sc.pl.umap(concat_h5ad,color=batch_name,save='_batch.png',palette=batch_palette,size=fontsize)
        sc.pl.umap(concat_h5ad,
            color=cluster,
            save='_{cluster}.png'.format(cluster=cluster),
            palette=cluster_palette,
            size=fontsize)
        if run_tsne:
            print("argument run_tsne is {run_tsne},running tsne now!".format(run_tsne=run_tsne))
            sc.tl.tsne(concat_h5ad,n_pcs=n_comps,use_rep="X_pca_harmony")
            sc.pl.tsne(concat_h5ad,color=batch_name,save='_batch.png',palette=batch_palette,size=fontsize)
            sc.pl.tsne(concat_h5ad,
                color=cluster,
                save='_{cluster}.png'.format(cluster=cluster),
                palette=cluster_palette,
                size=fontsize)
    return concat_h5ad

def output_data(concat_h5ad,prefix,cluster,outdir,run_tsne):
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
    # all genes and their gene ids
    gene_ids_ofile = os.path.join(data_path,prefix + '_gene_ids.txt')
    print(concat_h5ad.var)
    concat_h5ad.var['gene_ids'].to_csv(gene_ids_ofile,index=True,index_label='Gene')
    # tsne
    print(concat_h5ad)
    if run_tsne and ('X_tsne' in concat_h5ad.obsm):
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

@profile  
def integrate_pipeline(h5ad_dict,method,cluster,outdir,prefix,diff_method,gene_ids_df,
    run_tsne,resolution,barcodes,subset_column,subset_groups,high_var_genes_to_use,n_comps,
    fontsize,not_vst,output_rds,data_layer,batch_name):
    """
    pipeline to integrate scRNA data
    ingest must have the same variables,so we just do not supply this function for now
    harmony is the proposed method to use
    key args : 
        number of components to use for pca
        resolution
        highly variable genes number
    """
    if isinstance(h5ad_dict,dict): 
        #rename the observations by index
        for idx,(_,h5ad) in enumerate(h5ad_dict.items()):
            h5ad.obs_names = np.array([i.replace('-1','-{idx}'.format(idx=idx+1)) for i in h5ad.obs_names])
        #concatenation
        concat_h5ad = ad.concat(h5ad_dict.values(),label=batch_name,keys=h5ad_dict.keys(),join='outer')
    elif isinstance(h5ad_dict,ad.AnnData):
        concat_h5ad = h5ad_dict
        if data_layer == "raw":
            # As the following analysis is based on the data matrix
            # so if input h5ad .X stores scale.data and .raw stores data
            # .raw should be transmitted to .X
            print("data_layer is raw,so transfer adata.raw to adata.X...")
            # Using .to_adata() method can not retain attributes only in previous concat_h5ad
            # like layers and varm, obsp etc. So here we directly specify .X to .raw.X
            if concat_h5ad.raw.X.shape == concat_h5ad.X.shape:
                if isinstance(concat_h5ad.raw.X,np.ndarray):
                    concat_h5ad.X = concat_h5ad.raw.X.copy()
                else:
                    concat_h5ad.X = concat_h5ad.raw.X.toarray()
            else:
                print('Different shapes of .raw and .X,so using to_adata to transfer data matrix')
                concat_h5ad = concat_h5ad.raw.to_adata()
            print(concat_h5ad)

    # subset barcodes
    if barcodes:
        concat_h5ad = concat_h5ad[barcodes,:]
    if subset_column and subset_groups:
        print('subset concat anndata by {subset_column},subset groups are {subset_groups}'.format(
            subset_column = subset_column, subset_groups = subset_groups,
        ))
        concat_h5ad = concat_h5ad[concat_h5ad.obs[subset_column].isin(subset_groups)]
    # highly variable genes run flavor=seurat_v3 if layer['counts'] exists
    if hasattr(concat_h5ad,"layers") and ('counts' in concat_h5ad.layers) and (not not_vst):
        print('concat_ad.layers["counts"] exists,so run flavor=seurat_v3')
        sc.pp.highly_variable_genes(concat_h5ad,
            layer="counts",
            flavor="seurat_v3",
            n_top_genes=high_var_genes_to_use)
    else:
        sc.pp.highly_variable_genes(concat_h5ad,n_top_genes=high_var_genes_to_use)
    highly_var_number = len(concat_h5ad.var_names[concat_h5ad.var['highly_variable']])
    print("Number of highly variable genes is {}".format(highly_var_number))
    if high_var_genes_to_use > highly_var_number: # solve the number of highly variable problem
        print("[Warning]input highly variable genes not equal with return of highly_variable_genes,add 1 to high_var_genes_to_use")
        high_var_genes_to_use += 1
        sc.pp.highly_variable_genes(concat_h5ad,n_top_genes=high_var_genes_to_use)
        highly_var_number = len(concat_h5ad.var_names[concat_h5ad.var['highly_variable']])
        print("Number of highly variable genes is {} right now".format(highly_var_number))
        
    # add raw and scale data 
    concat_h5ad.raw = concat_h5ad
    sc.pp.scale(concat_h5ad)
    # remove batch effect with bbknn/ingest/harmony
    concat_h5ad = remove_batch_effect(concat_h5ad,method,cluster,run_tsne,resolution,
        n_comps,fontsize,batch_name)
    # add gene_ids
    # concat_h5ad.var['gene_ids'] = gene_ids_df.loc[concat_h5ad.var_names,'gene_ids']
    ofile = os.path.join(outdir,'h5ads',prefix + "_mid1" +'.h5ad')
    concat_h5ad.write(ofile)
    diff_genes = set(concat_h5ad.var_names) - set(gene_ids_df['genenames'])
    print('Difference genes between concat h5ad and config is {}'.format(diff_genes))
    #BUG: rename function will be get anndata can not be saved as .h5ad normally!
    # if "gene_ids" in concat_h5ad.var: 
        # print("gene_ids is in adata.var,so rename it to gene_ids_old")
        # concat_h5ad.var = concat_h5ad.var.rename(columns={"gene_ids":"gene_ids_old"})
    if "gene_ids" not in concat_h5ad.var:
        concat_h5ad.var = concat_h5ad.var.merge(gene_ids_df,how="left",right_index=True,left_index=True)
    ofile = os.path.join(outdir,'h5ads',prefix + "_mid2" +'.h5ad')
    concat_h5ad.write(ofile)
    # diff  use_raw=None means use_raw = True
    sc.tl.rank_genes_groups(concat_h5ad,cluster,use_raw=None,method=diff_method,pts=True)
    # output data
    output_data(concat_h5ad,prefix,cluster,outdir,run_tsne)
    # save h5ad
    ofile = os.path.join(outdir,'h5ads',prefix + '.h5ad')
    concat_h5ad.write(ofile,compression="gzip")
    # transfer h5ad to rds
    if output_rds:
        rds_ofile = ofile.replace('.h5ad','.rds')
        transfer_h5ad(ofile,rds_ofile)

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
    concat_h5ad = args.get('concat_h5ad')
    if config:
        h5ad_obj = parse_config(config,outdir)
    elif concat_h5ad:
        h5ad_dir = os.path.join(outdir,'h5ads')
        if not os.path.exists(h5ad_dir):
            os.makedirs(h5ad_dir)
        h5ad_obj = sc.read(concat_h5ad)
    if (not config) and (not concat_h5ad):
        raise Exception("You have to indicate config or concat_h5ad")
    print('loading h5ad completion!')
    # exit(0)
    # global settings
    figsize = tuple(int(i) for i in args.get('figsize').split(','))
    print('figsize is {}'.format(figsize))
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=args.get('dpi'), 
        frameon=True, 
        figsize=figsize, 
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
    #outdir = args.get('outdir')
    diff_method = args.get('diff_method')
    gene_ids = args.get('gene_ids')
    gene_ids_df = read_gene_ids(gene_ids)
    run_tsne = args.get('use_tsne')
    print("run_tsne is {run_tsne}".format(run_tsne=run_tsne))
    resolution = args.get('resolution')
    barcodes_file = args.get('barcodes')
    high_var_genes_to_use = args.get('highly_variable_genes')
    n_comps = args.get('n_comps')
    fontsize=args.get("fontsize")
    print('fontsize is %s' % fontsize)
    barcodes = []
    if barcodes_file:
        with open(barcodes_file,'r') as indata:
            if barcodes_file.endswith('csv'):
                sep=','
            else:
                sep="\t"
            title = indata.readline()
            tp = TitleParser(title,sep)
            for line in indata:
                line_list = line.strip().split(sep)
                barcode = tp.get_field(line_list,'Barcode')
                barcodes.append(barcode)
        print(barcodes[:5])
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    not_vst = args.get('not_vst')
    output_rds = args.get('output_rds')
    subset_column = args.get('subset_column')
    subset_groups = args.get('subset_groups')
    data_layer = args.get('data_layer')
    batch_name = args.get('batch_name')
    if subset_groups:
        subset_groups = [i.strip() for i in args.get('subset_groups').split(',')]
    print("Run pipeline")
    integrate_pipeline(h5ad_obj,method,cluster,outdir,prefix,diff_method,
        gene_ids_df,run_tsne,resolution,barcodes,subset_column,subset_groups,
        high_var_genes_to_use,n_comps,fontsize,not_vst,output_rds,data_layer,
        batch_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ## IO arguments
    parser.add_argument('--config','-c',help='config file with samplenames and their h5ad/rds files')
    # if input is h5ad
    parser.add_argument('--concat_h5ad','-ch',help='if you already have a h5ad file with all samples')
    parser.add_argument('--data_layer','-dl',choices=["X","raw"],default="raw",help='indicating where the data matrix is saved choosing from %(choices)s')
    
    parser.add_argument('--prefix','-p',help='Integration prefix to use')
    parser.add_argument('--gene_ids','-g',
        default='/TJPROJ5/SC/pipeline/scRNA/pipeline1.0/lib/05.Seurat/GRCh38_1.2.0_genes.tsv', 
        help='file to save the Ensembl ID of genes [default:%(default)s]')
    parser.add_argument('--outdir','-o',help='Output directory')
    parser.add_argument('--output_rds','-t',action='store_true',help='transfer result h5ad to rds')
    
    ## adjustment arguments
    parser.add_argument('--batch_name','-bn',default='batch',help='name of batch column, default:%(default)s')
    parser.add_argument('--method','-m',default='harmony',choices=['bbknn','harmony'],
        help='method for batch effect remove choose from %(choices)s [defalut:%(default)s]')
    parser.add_argument('--cluster',default='louvain',choices=['louvain','leiden'],
        help='cluster method to use choose from %(choices)s [defalut:%(default)s]')
    parser.add_argument('--diff_method','-d',default="wilcoxon",
        choices=["logreg","t-test","wilcoxon","t-test_overestim_var"],
        help='method to rank genes choose from %(choices)s [default:%(default)s]')
    parser.add_argument('--not_vst',default=None,action="store_true",
        help='do NOT use vst method of seurat_v3 to extract highly variable genes method to use,default is automatically selecting by\
    whether counts data exist in adata.layers["counts"]')
    parser.add_argument('--use_tsne','-u',action='store_true',help='run tsne at the same time,only useful to bbknn and harmony')
    parser.add_argument('--resolution','-r',type=float,default=1,help='resolution for cluster method louvain and leiden [default:%(default)s]')
    parser.add_argument('--highly_variable_genes','-v',type=int,default=2000,help='number of highly variable genes to use')
    parser.add_argument('--n_comps','-n',type=int,default=50,help='Number of principal components to compute. Defaults to 50')
    ## subset arguments
    parser.add_argument('--barcodes','-b',help='a file contains barcodes to analyse')
    parser.add_argument('--subset_column','-sc',help='you can subset anndata with indicated column and groups eg:louvain')
    parser.add_argument('--subset_groups','-sg',help='which groups you want to analysis in subset columns eg:1,3,5')
    ## plot arguments
    parser.add_argument('--dpi',type=int,default=100,help='dpi used to plot figures')
    parser.add_argument('--figsize',default="5,5",help='width and height of figure size,default:%(default)s')
    parser.add_argument('--fontsize',type=float,default=None,help='figure point size. If None, is automatically computed as 120000 / n_cells.')
    args = vars(parser.parse_args())
    main(args)
