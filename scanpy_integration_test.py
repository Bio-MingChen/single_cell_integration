#!/usr/bin/env python
# -*- coding=utf-8 -*-

import scanpy as sc
import pandas as pd
import seaborn as sns


sc.settings.verbosity = 1
sc.logging.print_versions()
sc.settings.set_figure_params(
    dpi=80, frameon=False, facecolor='white')
adata_all = sc.read_h5ad('../data/objects-pancreas/pancreas_cca.h5ad')
cca_h5 = sc.read_h5ad('../data/objects-pancreas/pancreas_cca.h5ad')
print(adata_all)
print(cca_h5)
#print(adata_all.shape)
#counts = adata_all.obs.celltype.value_counts()
