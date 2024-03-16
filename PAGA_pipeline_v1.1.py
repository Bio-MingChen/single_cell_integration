#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
from memory_profiler import profile
# do not show figure in window Linux mode
import matplotlib
matplotlib.use('Agg')

@profile
def main(args):
    #basic setting
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=100, 
        frameon=True, 
        figsize=(5, 5), 
        facecolor='white',
        )
    sc.settings.autoshow = False
    odir = args.get('odir')
    outfigure = os.path.join(odir,'figures')
    sc.settings.figdir = outfigure
    cluster_tag = args.get('cluster_tag')
    # input data
    h5ad_file = args.get('infile')
    h5ad = sc.read(h5ad_file)
    print(h5ad)
    # subset data
    clusters = args.get('clusters')
    if clusters:
        clusters = [c.strip() for c in clusters.split(',')]
        h5ad = h5ad[h5ad.obs[cluster_tag].isin(clusters),:]
        
    # subset by barcodes
    barcodes_file = args.get('barcodes')
    if barcodes_file:
        barcodes = pd.read_csv(barcodes_file)['Barcode'].tolist() # .tolist or .to_numpy
        h5ad = h5ad[barcodes,:]
    t = args.get("prefix")
    # run paga
    n_pcs = args.get('n_pcs')
    n_neighbors = args.get('n_neighbors')
    if args.get('re_analysis'):
        sc.tl.pca(h5ad)
        sc.pp.neighbors(h5ad, n_pcs=n_pcs, n_neighbors=n_neighbors)
        sc.tl.umap(h5ad)
        sc.tl.louvain(h5ad)
        sc.tl.diffmap(h5ad,n_comps=n_pcs)
        print("re-run pca and umap,cluster with louvain method,so designating louvain as cluster_tag!")
        cluster_tag = "louvain"
        print(h5ad.obs.louvain)
    else:
        # I am afraid diffmap have to run before subsetting the cells
        # which is necessary for the dpt
        sc.pp.neighbors(h5ad, n_pcs=n_pcs, n_neighbors=n_neighbors)
        sc.tl.diffmap(h5ad,n_comps=n_pcs)

    sc.pl.umap(h5ad,color=cluster_tag,palette=sc.pl.palettes.godsnot_102,size=10,save='_{}.png'.format(t))
    sc.pl.umap(h5ad,color=cluster_tag,palette=sc.pl.palettes.godsnot_102,size=10,save='_{}.pdf'.format(t))
    sc.tl.paga(h5ad, groups=cluster_tag)
    sc.pl.paga(h5ad, threshold=0.03, 
    # layout : {‘fa’, ‘fr’, ‘rt’, ‘rt_circular’, ‘drl’, ‘eq_tree’, Ellipsis}, None (default: None)
    # Plotting layout that computes positions. 
    # 'fa' stands for “ForceAtlas2”, 
    # 'fr' stands for “Fruchterman-Reingold”, 
    # 'rt' stands for “Reingold-Tilford”, 
    # 'eq_tree' stands for “eqally spaced tree”. 
    # All but 'fa' and 'eq_tree' are igraph layouts. 
    # All other igraph layouts are also permitted.
    #  See also parameter pos and draw_graph().
        layout='fa', 
        node_size_scale=1,
        node_size_power=0.5,
        max_edge_width=1,
        show=False,save='_{}_tag.png'.format(t))
    sc.pl.paga(h5ad, threshold=0.03, 
        layout='fa', 
        node_size_scale=1,
        node_size_power=0.5,
        max_edge_width=1,
        show=False,save='_{}_tag.pdf'.format(t))

    sc.pl.paga_compare(
        h5ad, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
        legend_fontsize=8, fontsize=5, edges=True, save='_{}_compare.png'.format(t))
    sc.pl.paga_compare(
        h5ad, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
        legend_fontsize=8, fontsize=5, edges=True, save='_{}_compare.pdf'.format(t))

    sc.tl.draw_graph(h5ad, init_pos='paga', layout='fa', maxiter=50)
    sc.pl.draw_graph(h5ad, color=cluster_tag,  layout='fa', legend_loc='on data',size=10,
        legend_fontsize=5, show=False,save='_{}_paga.png'.format(t))
    sc.pl.draw_graph(h5ad, color=cluster_tag,  layout='fa', legend_loc='on data',size=10,
        legend_fontsize=5, show=False,save='_{}_paga.pdf'.format(t))

    sc.tl.umap(h5ad, init_pos='paga', maxiter=100, min_dist=1)
    sc.pl.umap(h5ad, color=cluster_tag, legend_loc='on data', legend_fontsize=8,palette=sc.pl.palettes.godsnot_102,size=10,
        save='_{}_paga.png'.format(t))
    sc.pl.umap(h5ad, color=cluster_tag, legend_loc='on data', legend_fontsize=8,palette=sc.pl.palettes.godsnot_102,size=10,
        save='_{}_paga.pdf'.format(t))

    # default select the first category as root
    start_cluster = args.get('start_cluster')
    start_barcode = args.get('start_barcode')
    if (not start_cluster) and (start_barcode):
        print('set start barcode to {}'.format(start_barcode))
        h5ad.uns['iroot'] = np.flatnonzero(h5ad.obs.index == start_barcode)[0]
    else:
        if (not start_cluster):
            start_cluster = h5ad.obs[cluster_tag].cat.categories[0]
        print('{start_cluster} is used as start point to run paga'.format(start_cluster=start_cluster))
        h5ad.uns['iroot'] = np.flatnonzero(h5ad.obs[cluster_tag].isin(start_cluster.split(',')))[0]
    print(h5ad.obs[cluster_tag].cat.categories)
    print('iroot row number is {},barcode information: {}'.format(h5ad.uns['iroot'],h5ad.obs.iloc[h5ad.uns['iroot'],:]))
    
    # sc.tl.diffmap(h5ad,n_comps=10)
    sc.tl.dpt(h5ad)
    sc.pl.draw_graph(h5ad, color=[cluster_tag,'dpt_pseudotime'], layout="fa", legend_loc='on data',legend_fontsize=8,
        palette=sc.pl.palettes.godsnot_102,size=10, save='_{}_dpt.png'.format(t))
    sc.pl.draw_graph(h5ad, color=[cluster_tag,'dpt_pseudotime'], layout="fa", legend_loc='on data',legend_fontsize=8,
        palette=sc.pl.palettes.godsnot_102,size=10, save='_{}_dpt.pdf'.format(t))
    h5ad_outdir = os.path.join(odir,'h5ads')
    if not os.path.exists(h5ad_outdir):
        os.makedirs(h5ad_outdir)
    h5ad_ofile = os.path.join(h5ad_outdir,'{}.h5ad'.format(t))
    h5ad.write(h5ad_ofile)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',help='h5ad')
    parser.add_argument('-p','--prefix',help='prefix to use as output')
    parser.add_argument('-c','--cluster_tag',default="louvain",help='which column in adata.obs to use as cluster,default:%(default)s')
    parser.add_argument('-s','--start_cluster',default="0",help='which cluster to use as start %(default)s')
    parser.add_argument('-t','--start_barcode',default=None,help='indicating a barcode to start if --start_cluster is None')
    parser.add_argument('-o','--odir',help='output directory')
    parser.add_argument('-r','--re_analysis',action='store_true',
        help='re-run pca and cluster analysis for data')
    parser.add_argument('-l','--clusters',help='indicating which clusters to use,eg:1,3,5')
    parser.add_argument('-e','--n_neighbors', default=4,type=int,help='number of neighbors to consider in KNN default is 4 which is same with paga tutorial')
    parser.add_argument('-n','--n_pcs',type=int,default=20,help='number of pcs to use default:%(default)s ')
    parser.add_argument('-b','--barcodes',help='indicating which cells to use,a file contains barcodes with "Barcode" as title')
    args = vars(parser.parse_args())
    main(args)
