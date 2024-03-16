#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import scanpy as sc
import pandas as pd

if len(sys.argv) != 3:
    print('Usage: python output_var_table.py h5ad ofile')
    exit(0)

h5ad = sys.argv[1]
ofile = sys.argv[2]

indata = sc.read(h5ad)

indata.var['gene_ids'].to_csv(ofile,index=True,index_label='gene')
print('output gene ids to {ofile}'.format(ofile=ofile))