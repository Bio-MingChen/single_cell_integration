library(Seurat)

# This script copied from https://satijalab.org/seurat/articles/essential_commands.html
# I use this script to generate seurat object with analysis results of standard pipeline quickly
# only seurat object will be output 

args <- commandArgs(trailingOnly = TRUE)
# two format of data can be add cellranger or HDF5 format
# for HDF5 file suffixed with .h5 
# for cellranger output format: e.g: "~/Downloads/pbmc3k/filtered_gene_bc_matrices/hg19/"
datadir <- args[1] 
ofile <- args[2]
if (endsWith(datadir,".h5")) { # endsWith is availble for R verson more than 3.3
    pbmc.counts <- Read10X_h5(datadir)
} else {
    pbmc.counts <- Read10X(data.dir = datadir)
}
pbmc <- CreateSeuratObject(counts = pbmc.counts)
pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
pbmc <- FindNeighbors(object = pbmc)
pbmc <- FindClusters(object = pbmc)
# pbmc <- RunTSNE(object = pbmc)
pbmc <- RunUMAP(object = pbmc, dims = 1:20)
saveRDS(pbmc,file=ofile)