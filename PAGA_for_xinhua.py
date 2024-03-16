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
    # input data
    h5ad_file = args.get('infile')
    normal_barcode_file = args.get('normal_barcode')
    cancer_barcode_file = args.get('cancer_barcode')
    all_h5ad = sc.read(h5ad_file)
    # I am afraid diffmap have to run before subsetting the normal and cancer cells
    # which is necessary for the dpt
    sc.tl.diffmap(all_h5ad,n_comps=10)
    normal_barcodes = []
    cancer_barcodes = []
    with open(normal_barcode_file,'r') as indata:
        title = indata.readline()
        for line in indata:
            normal_barcodes.append(line.strip())
    with open(cancer_barcode_file,'r') as indata:
        title = indata.readline()
        for line in indata:
            cancer_barcodes.append(line.strip())
    normal_h5ad = all_h5ad[normal_barcodes]
    print(normal_h5ad)
    cancer_h5ad = all_h5ad[cancer_barcodes]
    print(cancer_h5ad)
    # run paga
    for t,h5ad in zip(['normal','cancer'],[normal_h5ad,cancer_h5ad]):
        if args.get('re_analysis'):
            sc.tl.pca(h5ad)
            sc.pp.neighbors(h5ad, n_pcs=10)
            sc.tl.umap(h5ad)
            sc.tl.louvain(h5ad)
            sc.tl.diffmap(h5ad,n_comps=10)
        sc.pl.umap(h5ad,color='louvain',save='_{}.png'.format(t))
        sc.pl.umap(h5ad,color='louvain',save='_{}.pdf'.format(t))
        sc.tl.paga(h5ad, groups='louvain')
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
            show=False,save='_{}_default.png'.format(t))
        sc.pl.paga(h5ad, threshold=0.03, 
            layout='fa', 
            node_size_scale=1,
            node_size_power=0.5,
            max_edge_width=1,
            show=False,save='_{}_default.pdf'.format(t))
        sc.pl.paga(h5ad, threshold=0.03, 
            layout='fa',
            fontsize=5,
            node_size_scale=0.5,
            node_size_power=0.5,
            max_edge_width=0.7,
            show=False,save='_{}_low.png'.format(t))
        sc.pl.paga(h5ad, threshold=0.03, 
            layout='fa',
            fontsize=5,
            node_size_scale=0.5,
            node_size_power=0.5,
            max_edge_width=0.7,
            show=False,save='_{}_low.pdf'.format(t))
        sc.tl.umap(h5ad, init_pos='paga', maxiter=100, min_dist=1)
        sc.pl.umap(h5ad, color='louvain', legend_loc='on data', legend_fontsize=5, 
            save='_{}_paga.png'.format(t))
        sc.pl.umap(h5ad, color='louvain', legend_loc='on data', legend_fontsize=5, 
            save='_{}_paga.pdf'.format(t))
        sc.tl.draw_graph(h5ad, init_pos='paga', layout='fa', maxiter=50)
        sc.pl.draw_graph(h5ad, color='louvain',  legend_loc='on data',
            legend_fontsize=5, show=False,save='_{}_paga.png'.format(t))
        sc.pl.draw_graph(h5ad, color='louvain',  legend_loc='on data',
            legend_fontsize=5, show=False,save='_{}_paga.pdf'.format(t))
        # default select the first category as root
        first_cluster = h5ad.obs.louvain.cat.categories[0]
        print(h5ad.obs.louvain.cat.categories)
        h5ad.uns['iroot'] = np.flatnonzero(h5ad.obs['louvain'] == first_cluster)[0]
        # sc.tl.diffmap(h5ad,n_comps=10)
        sc.tl.dpt(h5ad)
        sc.pl.draw_graph(h5ad, color='dpt_pseudotime', legend_loc='on data',
            save='_{}_dpt.png'.format(t))
        sc.pl.draw_graph(h5ad, color='dpt_pseudotime', legend_loc='on data',
            save='_{}_dpt.pdf'.format(t))
        h5ad_outdir = os.path.join(odir,'h5ads')
        if not os.path.exists(h5ad_outdir):
            os.makedirs(h5ad_outdir)
        h5ad_ofile = os.path.join(h5ad_outdir,'{}.h5ad'.format(t))
        h5ad.write(h5ad_ofile)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--infile',help='h5ad of xinhua project')
    parser.add_argument('-n','--normal_barcode',help='normal barcodes to use')
    parser.add_argument('-c','--cancer_barcode',help='cancer barcodes to use')
    parser.add_argument('-o','--odir',help='output directory')
    parser.add_argument('-r','--re_analysis',action='store_true',
        help='re-run pca and cluster analysis for subset data')
    args = vars(parser.parse_args())
    main(args)