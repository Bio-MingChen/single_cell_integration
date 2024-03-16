#!/usr/bin/env python
# -*- coding=utf-8 -*-
import sys
import numpy as np
import pandas as pd

infile = sys.argv[1]
fraction = sys.argv[2]
ofile = sys.argv[3]
if len(sys.argv) != 4:
    print('''
    Usage: python random_sample_by_group.py infile fraction[0~1] ofile
    ''')
df = pd.read_csv(infile,sep='\t')
print(df.head())
sample_df = df.groupby(['cluster']).sample(frac=float(fraction))
sample_df.to_csv(ofile,sep='\t',header=True,index=False)
print('save {fraction}% percent sample dataframe to {ofile}'.format(
    fraction = float(fraction)*100,ofile=ofile
))