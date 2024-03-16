library(optparse)
library(scHCL)
library(DESeq2)
library(pheatmap)

option_list <- list(
    make_option(c('-i','--infile'),help='infile name'),
    make_option(c('-o','--opre'),default='normalized_counts',help='prefix of ofile name'),
    make_option(c('-n','--number'),default=3,type="integer",
        help='top n related cells will use default is 3'),
    make_option(c('-c','--cutoff'),default=NULL,
        help='correlation less than specified value will be remove defalut is NULL')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
infile <- args$infile
#expression matrix
infile_vector <- unlist(strsplit(infile,".",fixed=TRUE))
if (infile_vector[length(infile_vector)] == "rds") {
    seurat_obj <- readRDS(infile)
    emtx <- seurat_obj@assays$RNA@counts
} else {
    emtx <- as.matrix(read.table(infile,header=TRUE,sep="\t",row.names=1,check.names=FALSE,
        stringsAsFactors=FALSE))
}

hcl_result <- scHCL(scdata = emtx, numbers_plot = args$number)
if (!is.null(args.cutoff)) {
    cutoff <- as.numeric(args.cutoff)
    hcl_result$scHCL_probility
}
## do not consider to output matrix of pearson currently considering its large storage
#pearson_matrix_ofile <- sprintf('%s_pearson.xls',args$opre)
#write.table(hcl_result$cors_matrix,file=pearson_matrix_ofile,sep='\t',quote=FALSE,col.names=NA)

celltype_ofile <- sprintf('%s_celltype.xls',args$opre)
heatmap_ofile <- sprintf('%s_pheatmap.png',args$opre)
write.table(hcl_result$scHCL_probility,file=celltype_ofile,sep='\t',quote=FALSE,row.names=FALSE)

#heatmap for top number(default3) celltype
gettissue <- function(x,Num=3){
  top_markers <-order(x,decreasing = T)[1:Num]
  return(top_markers)
}
cors <- hcl_result$cors_matrix
cors_index <- apply(cors,2,gettissue,args$number)
cors_index <- sort(unique(as.integer(cors_index)))
data <- cors[cors_index,]
opre <- basename(args$opre)
opre <- unlist(strsplit(opre,'_'))[1]
get_plot_dims <- function(heat_map) {
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}
celltype_heatmap <- pheatmap(data,cellheight=10,show_colnames=FALSE,main=sprintf('%s Celltype Prediction',opre))
plot_dims <- get_plot_dims(celltype_heatmap)
png(heatmap_ofile, height = plot_dims$height, width = plot_dims$width, units = "in", res = 300)
celltype_heatmap
dev.off()
