library(optparse)
library(scHCL)
library(pheatmap)
library(ggplot2)
library(dplyr)

option_list <- list(
    make_option(c('-i','--infile'),help='infile name'),
    make_option(c('-o','--opre'),default='normalized_counts',help='prefix of ofile name'),
    make_option(c('-n','--number'),default=3,type="integer",
        help='top n related cells will use default is 3'),
    make_option(c('--stat_db'), 
        default="/TJPROJ6/SC/personal_dir/chenming/research/scHCL_cell_identity/data/scHCL_ref_colname_modify2.csv",
        help="config file used to compute statistics"),
    make_option(c('-c','--cutoff'),default=NULL,
        help='correlation less than specified value will be remove defalut is NULL')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
infile <- args$infile

## get expression matrix from rds or plain text which 
## columns are barcodes and rows are genenames
infile_vector <- unlist(strsplit(infile,".",fixed=TRUE))
if (infile_vector[length(infile_vector)] == "rds") {
    seurat_obj <- readRDS(infile)
    emtx <- seurat_obj@assays$RNA@counts
} else {
    emtx <- as.matrix(read.table(infile,header=TRUE,sep="\t",row.names=1,check.names=FALSE,
        stringsAsFactors=FALSE))
}
# run scHCL
hcl_result <- scHCL(scdata = emtx, numbers_plot = args$number)
if (!is.null(args$cutoff)) {
    cutoff <- as.numeric(args$cutoff)
    hcl_result$scHCL_probility <- filter(hcl_result$scHCL_probility,Score > cutoff)
}

## do not consider to output matrix of pearson currently considering its large storage
#pearson_matrix_ofile <- sprintf('%s_pearson.xls',args$opre)
#write.table(hcl_result$cors_matrix,file=pearson_matrix_ofile,sep='\t',quote=FALSE,col.names=NA)

# compute counts by celltype or organ 
# annotation file scHCL offered is somewhat confused because of its format like 
# Monocyte_FCGR3A.high.Adult.Bone.Marrow.CD34P. so we need to extract celltype and organ. 
# stat_db default is a file where celltype and organ for each barcode separately showed
# current version is contributed by caodongyan@novogene.com

celltype_ofile <- sprintf('%s_celltype_eachcell.xls',args$opre)
heatmap_ofile <- sprintf('%s_pheatmap.png',args$opre)
stat_db <- read.table(args$stat_db,sep=",",header=TRUE)
stat_df <- merge(
    x = hcl_result$scHCL_probility,
    y = stat_db,
    by.x = "Cell type",
    by.y = "Cell.type",
    all.x = TRUE, # left join
    )

colnames(stat_df) <- c("scHCLAnno", "Cell", "Score", "Celltype", "Organ", "OrganCelltype")

stat_df <- stat_df %>%
    arrange(Cell,Score)
write.table(stat_df[,1:5],file=celltype_ofile,sep='\t',quote=FALSE,row.names=FALSE)

# extract top1 correlation for each barcode and 
# compute cells' counts by scHCLAnno,Celltype,Organ column
# finally make barplot of top10 related types
top1_df <- stat_df %>%
        group_by(Cell) %>%
        top_n(1,Score)
# Celltype
celltype_df <- top1_df %>%
        group_by(Celltype) %>%
        count() %>%
        ungroup() %>%
        arrange(desc(n)) %>%
        rename(Counts=n)

celltype_df_ofile <- sprintf('%s_celltype_stat.xls',args$opre)
write.table(celltype_df,celltype_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)

celltype_bar <- ggplot(data=top_n(celltype_df,10,Counts)) +
    geom_bar(aes(x=reorder(Celltype,-Counts),y=Counts,fill=Celltype),stat="identity") +
    theme_bw() +
    theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
    labs(x="Celltype",y="Counts",title="Top 10 related celltype")

celltype_bar_ofile <- sprintf('%s_celltype_bar.png',args$opre)
ggsave(celltype_bar_ofile,celltype_bar)

#Organ
organ_df <- top1_df %>%
        group_by(Organ) %>%
        count() %>%
        ungroup() %>%
        arrange(desc(n)) %>%
        rename(Counts=n)

organ_df_ofile <- sprintf('%s_organ_stat.xls',args$opre)
write.table(organ_df,organ_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)

organ_bar <- ggplot(data=top_n(organ_df,10,Counts)) +
    geom_bar(aes(x=reorder(Organ,-Counts),y=Counts,fill=Organ),stat="identity") +
    theme_bw() +
    theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
    labs(x="Organ",y="Counts",title="Top 10 related Organ")

organ_bar_ofile <- sprintf('%s_organ_bar.png',args$opre)
ggsave(organ_bar_ofile,organ_bar)

#scHCLAnno
scHCLAnno_df <- top1_df %>%
        group_by(scHCLAnno) %>%
        count() %>%
        ungroup() %>%
        arrange(desc(n)) %>%
        rename(Counts=n)

scHCLAnno_df_ofile <- sprintf('%s_scHCLAnno_stat.xls',args$opre)
write.table(scHCLAnno_df,scHCLAnno_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)

scHCLAnno_bar <- ggplot(data=top_n(scHCLAnno_df,10,Counts)) +
    geom_bar(aes(x=reorder(scHCLAnno,-Counts),y=Counts,fill=scHCLAnno),stat="identity") +
    theme_bw() +
    theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
    labs(x="scHCLAnno",y="Counts",title="Top 10 related scHCLAnno")

scHCLAnno_bar_ofile <- sprintf('%s_scHCLAnno_bar.png',args$opre)
ggsave(scHCLAnno_bar_ofile,scHCLAnno_bar)

## heatmap for top number(default3) celltype
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
