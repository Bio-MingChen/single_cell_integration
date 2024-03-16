import click
import numpy as np
import pandas as pd

@click.command()
@click.argument('barcode_file')
@click.argument('check_file')
def check_cluster_by_barcode(barcode_file,check_file):
    """\b
    check_cluster_by_barcode barcode_file check_file
    barcode column with title Barcode
    input should a csv file 
    """
    barcode = pd.read_csv(barcode_file)["Barcode"].to_list()
    check_df = pd.read_csv(check_file)
    clusters = pd.unique(check_df[check_df["Barcode"].isin(barcode)]["Cluster"])
    print(clusters)

if __name__ == "__main__":
    check_cluster_by_barcode()