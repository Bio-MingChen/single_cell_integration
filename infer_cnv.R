library(infercnv)

raw_counts_matrix=system.file("extdata", 
    "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv")
head(raw_counts_matrix)
annotations_file=system.file("extdata", 
    "oligodendroglioma_annotations_downsampled.txt", package = "infercnv")
head(annotations_file)
gene_order_file=system.file("extdata", 
    "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv")
ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")