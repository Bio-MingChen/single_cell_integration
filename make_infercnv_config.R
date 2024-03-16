library(optparse)
library(dplyr)
library(purrr)
library(Seurat)

option_list <- list(
    make_option(c('-c','--config'),help='config with sample name and rds'),
    make_option(c('-l','--cluster_info'),help='cluster info csv file with barcodes and cluster number'),
    make_option(c('-g','--gene_list'),help='a csv file contain genes to use with gene id'),
    make_option(c('-o','--opre'),help='prefix to output'),
    make_option(c('-r','--control_cluster'),help='clusters to use in cluster_info as control'),
    make_option(c('-a','--case_cluster'),help='clusters to use in cluster_info as case'),
    make_option(c("-u", "--use_id"), action="store_true",default=FALSE, help="use gene id rather than gene name to output"),
    make_option(c('--control_max_cell'),type='integer',default=50,help='max cells number to use as control for each cluster'),
    make_option(c('--case_split_cell'),type='integer',default=1000,help='max cells number to combine with contrl each time')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
# input all rds and extract raw counts matrix 
config <- read.table(args$config,header=TRUE,stringsAsFactors=FALSE)
gene_df <- read.table(args$gene_list,sep=',',header=TRUE,stringsAsFactors=FALSE)
total_raw_counts_list <- list()
for (i in 1:length(config$rds_path)) {
    sample_seurat <- readRDS(config$rds_path[i])
    raw_counts <- as.data.frame(sample_seurat@assays$RNA@data)
    replace_str <- sprintf('-%d',i)
    colnames(raw_counts) <- gsub('-1',replace_str,colnames(raw_counts))
    raw_counts <- raw_counts[rownames(raw_counts) %in% gene_df$gene,]
    if(args$use_id) { # use gene id rather than gene as rownames
        rownames(raw_counts) <- gene_df$gene_ids
    }
    total_raw_counts_list[[i]] <- raw_counts
}
# cbind matrix
# print(str(total_raw_counts_list))
total_raw_counts <- total_raw_counts_list %>% reduce(cbind)
# output indicated barcodes
cluster_info <- read.table(args$cluster_info,sep=',',header=TRUE,stringsAsFactors=FALSE)

##get the control cells
control_clusters = as.numeric(unlist(strsplit(args$control_cluster,',')))
case_clusters = as.numeric(unlist(strsplit(args$case_cluster,',')))
control_list = list()
for (i in 1:length(control_clusters)) {
    c = control_clusters[i]
    print(c)
    sub_cluster <- cluster_info$Barcode[cluster_info$Cluster == c]
    if(length(sub_cluster) > args$control_max_cell) {
        control_sample = sample(sub_cluster,args$control_max_cell)
    } else {
        control_sample = sub_cluster
    }
    control_list[[i]] = cluster_info[cluster_info$Barcode %in% control_sample,]
}
total_control <- control_list %>% reduce(rbind)
control_columns <- intersect(colnames(total_raw_counts),total_control$Barcode)
print(length(control_columns))
print(control_columns[1:5])
# control_matrix <- total_raw_counts[,control_columns]
#outout control barcodes list
ofile <- sprintf('%s_control_bar.csv',args$opre)
write.table(total_control,
    ofile,sep='\t',quote=FALSE,row.names=FALSE)

## get the case cells
get_counts <- function(barcodes) {
    output_list <- list()
    case_columns <- intersect(colnames(total_raw_counts),barcodes)
    chunk_raw_counts <- total_raw_counts[,c(control_columns,case_columns)]
    sub_cluster_info_df <- cluster_info[cluster_info$Barcode %in% c(control_columns,case_columns),]
    output_title <- c('X',colnames(chunk_raw_counts))
    output_df <- data.frame(X=rownames(chunk_raw_counts),chunk_raw_counts)
    colnames(output_df) <- output_title
    output_list[[1]] <- sub_cluster_info_df
    output_list[[2]] <- output_df
    return(output_list)
}
for (i in 1:length(case_clusters)) {
    c = case_clusters[i]
    case_df = cluster_info[cluster_info$Cluster == c,]
    if (nrow(case_df) > args$case_split_cell) {
        split_bars <- split(case_df$Barcode,ceiling(1:nrow(case_df)/args$case_split_cell))
        chunks_list <- lapply(split_bars,get_counts)
        for (j in 1:length(chunks_list)) {
            anno_ofile <- sprintf('%s_cluster%d_chunk%d_annotation.tsv',args$opre,i,j)
            matrix_ofile <- sprintf('%s_cluster%d_chunk%d_matrix.tsv',args$opre,i,j)
            output_list <- chunks_list[[j]]
            write.table(output_list[[1]],anno_ofile,sep='\t',quote=FALSE,row.names=FALSE)
            write.table(output_list[[2]],matrix_ofile,sep='\t',quote=FALSE,row.names=FALSE)
        }
    } else {
        output_list <- get_counts(case_df$Barcode)
        anno_ofile <- sprintf('%s_cluster%d_annotation.tsv',args$opre,i)
        matrix_ofile <- sprintf('%s_cluster%d_matrix.tsv',args$opre,i)
        write.table(output_list[[1]],anno_ofile,sep='\t',quote=FALSE,row.names=FALSE)
        write.table(output_list[[2]],matrix_ofile,sep='\t',quote=FALSE,row.names=FALSE)
    }
}
