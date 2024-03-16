library(optparse)
library(Seurat)
library(dplyr)
library(purrr)

option_list <- list(
    make_option(c('-c','--config'),help='config with sample name and rds'),
    make_option(c('-m','--max_cell'),type='integer',default=NA,help='max cells number to output'),
    make_option(c('-l','--cluster_info'),help='cluster info csv file with barcodes and cluster number'),
    make_option(c('-g','--gene_list'),help='a csv file contain genes to use with gene id'),
    make_option(c('-o','--opre'),help='prefix to output'),
    make_option(c("-u", "--use_id"), action="store_true",
                default=FALSE, help="use gene id rather than gene name to output"),
    make_option(c('-n','--cluster_numbers'),default=NA,help='clusters to use in cluster_info,all for one output,NA for every cluster')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
# input all rds and extract raw counts matrix 
config <- read.table(args$config,header=TRUE,stringsAsFactors=FALSE)
gene_df <- read.table(args$gene_list,sep=',',header=TRUE,stringsAsFactors=FALSE)
total_raw_counts_list <- list()
for (i in 1:length(config$rds_path)) {
    sample_seurat <- readRDS(config$rds_path[i])
    raw_counts <- as.data.frame(sample_seurat@assays$RNA@counts)
    replace_str <- sprintf('-%d',i)
    colnames(raw_counts) <- gsub('-1',replace_str,colnames(raw_counts))
    raw_counts <- raw_counts[rownames(raw_counts) %in% gene_df$Gene,]
    if(args$use_id) { # use gene id rather than gene as rownames
        rownames(raw_counts) <- gene_df$gene_ids
    }
    # if (!is.na(args$max_cell)) { # sample max cells
    #     barcodes <- colnames(raw_counts)
    #     if (args$max_cell < length(barcodes)) {
    #         sample_bar <- sample(barcodes,args$max_cell)
    #         raw_counts <- raw_counts[,sample_bar]
    #     }
    # }
    total_raw_counts_list[[i]] <- raw_counts
}
# cbind matrix
# print(str(total_raw_counts_list))

total_raw_counts <- total_raw_counts_list %>% reduce(cbind)
# output indicated barcodes
output_func <- function(total_raw_counts,cluster_info,ofile) {
output_columns <- intersect(colnames(total_raw_counts),cluster_info$Barcode)
total_raw_counts <- total_raw_counts[,output_columns]
output_title <- c('X',colnames(total_raw_counts))
output_df <- data.frame(X=rownames(total_raw_counts),total_raw_counts)
colnames(output_df) <- output_title
write.table(output_df,
    ofile,sep='\t',quote=FALSE,row.names=FALSE)
}

cluster_info <- read.table(args$cluster_info,sep=',',header=TRUE,stringsAsFactors=FALSE)
if(is.na(args$cluster_numbers)) {
    clusters <- unique(cluster_info$Cluster)
    for (c in clusters) {
        if (!is.na(args$max_cell)) { # sample max cells
            barcodes <- cluster_info$Barcode[cluster_info$Cluster == c]
            if (args$max_cell < length(barcodes)) {
                sample_bar <- sample(barcodes,args$max_cell)
                cluster_info_sub <- cluster_info[cluster_info$Barcode %in% sample_bar,]
                ofile <- sprintf('%s_cluster%s_max_cell%s.tsv',args$opre,c,args$max_cell)
            } else {
                cluster_info_sub <- cluster_info[cluster_info$Cluster == c,]
                ofile <- sprintf('%s_cluster%s.tsv',args$opre,c)
            }
            output_func(total_raw_counts,cluster_info_sub,ofile)
        }
    }
} else if(args$cluster_numbers == 'all') {
    if (!is.na(args$max_cell)) {
        cluster_info = cluster_info[cluster_info$Barcode %in% sample(cluster_info$Barcode,args$max_cell),]
    }
    output_func(total_raw_counts,cluster_info,ofile)
} else {
    clusters = as.numeric(unlist(strsplit(args$cluster_numbers,',')))
    for (c in clusters) {
        if (!is.na(args$max_cell)) { # sample max cells
            barcodes <- cluster_info$Barcode[cluster_info$Cluster == c]
            if (args$max_cell < length(barcodes)) {
                sample_bar <- sample(barcodes,args$max_cell)
                cluster_info_sub <- cluster_info[cluster_info$Barcode %in% sample_bar,]
                ofile <- sprintf('%s_cluster%s_max_cell%s.tsv',args$opre,c,args$max_cell)
            } else {
                cluster_info_sub <- cluster_info[cluster_info$Cluster == c,]
                ofile <- sprintf('%s_cluster%s.tsv',args$opre,c)
            }
            output_func(total_raw_counts,cluster_info_sub,ofile)
        }
    }
}

