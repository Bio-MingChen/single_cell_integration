import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc

h5ad_file = sys.argv[1]
outfile = sys.argv[2]
h5ad = sc.read(h5ad_file)

sc.pp.highly_variable_genes(h5ad,n_top_genes=2000)
highly_variable_genes = h5ad.var_names[h5ad.var['highly_variable'] == True].to_list()
with open(outfile,'w') as odata:
    for gene in highly_variable_genes:
        odata.write(gene+"\n")
