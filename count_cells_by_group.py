#!/usr/bin/env python
# -*- coding=utf-8 -*-
import sys
import numpy as np
import pandas as pd

infile = sys.argv[1]
ofile = sys.argv[2]
if len(sys.argv) != 3:
    print('''
    Usage: python count_cells_by_group.py infile ofile
    ''')
df = pd.read_csv(infile,sep='\t')
count_df = df.loc[:,['label','cluster']].groupby(['cluster']).count()
count_df.columns = ['cell_number']
print(count_df.head())

count_df.to_csv(ofile,sep='\t',header=True,index=True)
print('save stat info to {ofile}'.format(
    ofile=ofile
))