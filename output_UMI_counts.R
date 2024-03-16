library(Seurat)

seurat_data <- readRDS('/TJPROJ5/SC/SHOUHOU/X101SC20050308-Z01/Seurat/All_cancer/gene_bar/CRC_C_1.rds')

print(seurat_data)
print(str(seurat_data))

counts_df <- as.data.frame(seurat_data@assays$RNA@counts)
print(counts_df[1:5,1:5])
write.table(counts_df,'test.counts.tsv',sep='\t',quote=FALSE)
