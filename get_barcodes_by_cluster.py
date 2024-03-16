#!/usr/bin/env python
# -*- coding:utf-8 -*-
import sys
import pandas as pd

barcode_file = sys.argv[1]
cluster_file = sys.argv[2]
ofile = sys.argv[3]

b_df = pd.read_csv(barcode_file)
print(b_df.head())
if "," in cluster_file:
    clusters = [int(i) for i in cluster_file.split(",")]
else:
    clusters=[]
    with open(cluster_file,'r') as indata:
        for line in indata:
            clusters.append(line.strip())
print(clusters)
sub_b_df = b_df[b_df["Cluster"].isin(clusters)]
sub_b_df.to_csv(ofile,index=False)
