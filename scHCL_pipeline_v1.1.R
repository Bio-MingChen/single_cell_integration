library(optparse)
library(scHCL)
library(pheatmap)
library(ggplot2)
library(Seurat)
library(dplyr)

option_list <- list(
    make_option(c('-i','--infile'),help='infile name with .rds as suffix'),
    make_option(c('-o','--odir'),default='normalized_counts',help='output directory'),
    make_option(c('-m','--matrix'),default='data',help='which matrix to use,default is data'),
    make_option(c('-n','--number'),default=3,type="integer",
        help='top n related cells will use default is 3'),
    make_option(c('-t','--cluster_str'),default='seurat_clusters',help='column name of cluster'),
    make_option(c('-s','--subset_clusters'),help='subset clusters eg. 1,3,5'),
    make_option(c('-a','--all_clusters'),default=FALSE,action="store_true",help='treat all cells as one cluster to run scHCL'),
    make_option(c('--stat_db'), 
        default="/TJPROJ6/SC/personal_dir/chenming/research/scHCL_cell_identity/data/scHCL_ref_colname_modify2.csv",
        help="config file used to compute statistics"),
    make_option(c('-c','--cutoff'),default=NULL,
        help='correlation less than specified value will be remove defalut is NULL')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
infile <- args$infile

run_scHCL <- function(emtx,odir,opre) {
    # run scHCL
    hcl_result <- scHCL(scdata = emtx, numbers_plot = args$number)
    if (!is.null(args$cutoff)) {
        cutoff <- as.numeric(args$cutoff)
        hcl_result$scHCL_probility <- filter(hcl_result$scHCL_probility,Score > cutoff)
    }

    ## do not consider to output matrix of pearson currently considering its large storage
    #pearson_matrix_ofile <- sprintf('%s/pearson.xls',odir)
    #write.table(hcl_result$cors_matrix,file=pearson_matrix_ofile,sep='\t',quote=FALSE,col.names=NA)

    # compute counts by celltype or organ 
    # annotation file scHCL offered is somewhat confused because of its format like 
    # Monocyte_FCGR3A.high.Adult.Bone.Marrow.CD34P. so we need to extract celltype and organ. 
    # stat_db default is a file where celltype and organ for each barcode separately showed
    # current version is contributed by caodongyan@novogene.com

    celltype_ofile <- sprintf('%s/%s_celltype_eachcell.xls',odir,opre)
    heatmap_ofile <- sprintf('%s/%s_pheatmap.png',odir,opre)
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
    write.table(stat_df,file=celltype_ofile,sep='\t',quote=FALSE,row.names=FALSE)

    # extract top1 correlation for each barcode and 
    # compute cells' counts by scHCLAnno,Celltype,Organ column
    # finally plot barplot of top10 related types
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

    celltype_df_ofile <- sprintf('%s/%s_celltype_stat.xls',odir,opre)
    write.table(celltype_df,celltype_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)

    celltype_bar <- ggplot(data=top_n(celltype_df,10,Counts)) +
        geom_bar(aes(x=reorder(Celltype,-Counts),y=Counts,fill=Celltype),stat="identity") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
        labs(x="Celltype",y="Counts",title="Top 10 related celltype")

    celltype_bar_ofile <- sprintf('%s/%s_celltype_bar.png',odir,opre)
    ggsave(celltype_bar_ofile,celltype_bar,width=7,height=7)

    #Organ
    organ_df <- top1_df %>%
            group_by(Organ) %>%
            count() %>%
            ungroup() %>%
            arrange(desc(n)) %>%
            rename(Counts=n)

    organ_df_ofile <- sprintf('%s/%s_organ_stat.xls',odir,opre)
    write.table(organ_df,organ_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)

    organ_bar <- ggplot(data=top_n(organ_df,10,Counts)) +
        geom_bar(aes(x=reorder(Organ,-Counts),y=Counts,fill=Organ),stat="identity") +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
        labs(x="Organ",y="Counts",title="Top 10 related Organ")

    organ_bar_ofile <- sprintf('%s/%s_organ_bar.png',odir,opre)
    ggsave(organ_bar_ofile,organ_bar,width=7,height=7)

    #OrganCelltype
    scHCLAnno_df <- top1_df %>%
            group_by(OrganCelltype) %>%
            count() %>%
            ungroup() %>%
            arrange(desc(n)) %>%
            rename(Counts=n)

    scHCLAnno_df_ofile <- sprintf('%s/%s_OrganCelltype_stat.xls',odir,opre)
    write.table(scHCLAnno_df,scHCLAnno_df_ofile,sep="\t",quote=FALSE,row.names=FALSE)
    scHCLAnno_top10_df <- scHCLAnno_df %>% 
        top_n(10,Counts) %>%
        arrange(desc(Counts))
    
    scHCLAnno_bar <- ggplot(data = scHCLAnno_top10_df) +
        geom_bar(aes(x=reorder(OrganCelltype,-Counts),y=Counts,fill=OrganCelltype),stat="identity") +
        # coord_flip() +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 30,vjust = 0.5, hjust = 1),
            plot.title = element_text(hjust = 0.5)) + 
        labs(x="OrganCelltype",y="Counts",title="Top 10 related Organ Celltype")

    scHCLAnno_bar_ofile <- sprintf('%s/%s_OrganCelltype_bar.png',odir,opre)
    ggsave(scHCLAnno_bar_ofile,scHCLAnno_bar,width=7,height=7)

    ## heatmap for top number(default3) celltype
    gettissue <- function(x,Num=3){
    top_markers <-order(x,decreasing = T)[1:Num]
    return(top_markers)
    }
    cors <- hcl_result$cors_matrix
    cors_index <- apply(cors,2,gettissue,args$number)
    cors_index <- sort(unique(as.integer(cors_index)))
    data <- cors[cors_index,]
    get_plot_dims <- function(heat_map) {
    plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
    plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
    return(list(height = plot_height, width = plot_width))
    }
    celltype_heatmap <- pheatmap(data,cellheight=10,show_colnames=FALSE,main=sprintf('%s Celltype Prediction',opre))
    plot_dims <- get_plot_dims(celltype_heatmap)
    print(heatmap_ofile)
    png(filename = heatmap_ofile, height = plot_dims$height, width = plot_dims$width, units = "in", res = 300, type = "cairo")
    # celltype_heatmap
    pheatmap(data,cellheight=10,show_colnames=FALSE,main=sprintf('%s Celltype Prediction',opre))
    dev.off()
}


## pipeline to run scHCL
## run scHCL by indicated clusters in default otherwise if args$all_clusters,run all cells 

seurat_obj <- readRDS(infile)
if (!is.null(args$subset_clusters)) {
    clusters <- unlist(FetchData(seurat_obj,vars = args$cluster_str))
    sub_clusters <- unlist(strsplit(args$subset_clusters,",",fixed=TRUE))
    seurat_obj <- seurat_obj[,which(clusters %in% sub_clusters)]
}

dir.create(args$odir)

if (args$all_clusters) {
    emtx <- GetAssayData(seurat_obj, slot = args$matrix)
    run_scHCL(emtx,args$odir,"All")
} else {
        clusters <- FetchData(seurat_obj,vars = args$cluster_str)
        print(unique(clusters))
    for (i in unlist(unique(seurat_obj[[args$cluster_str]]))) {
        print(i)
        seurat_obj_sub <- seurat_obj[,which(clusters == i)]
        print(seurat_obj_sub)
        emtx_sub <- GetAssayData(seurat_obj_sub, slot = args$matrix)
        odir <- file.path(args$odir,sprintf("Cluster%s",i))
        dir.create(odir)
        run_scHCL(emtx_sub,odir,sprintf("Cluster%s",i))
    }
}

