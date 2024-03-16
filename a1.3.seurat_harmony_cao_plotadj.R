library(Seurat)

devtools::load_all('/TJPROJ5/SC/pipeline/10X_RNA/Development/Harmony/harmony-master')
#library(harmony)
library(dplyr)
library(stringr)
library(argparse)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(cowplot)

library(reticulate)
##use_python("/TJPROJ5/SC/software/miniconda3/bin/python")
##use_python("/NJPROJ1/RNA/software/python/Python-2.7.14/bin/python")



parser = ArgumentParser()
parser$add_argument("--path", help="/path/to/samples gene bar")
parser$add_argument("--compare", help="samplevssamplevssample,using vs as the split")
parser$add_argument("--gene_path", help="path to gene.tsv")
parser$add_argument("--species", help="GRCh38 or mm10")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--Nfeatures",help="nfeatures of FindVariableFeatures()", default='2000' )
parser$add_argument("--prefix", help="prefix of results")
parser$add_argument("--resolution", help="resolution for cluster",default='0.6')
parser$add_argument("--cca_use", help="cca_use for cluster",default='20')
parser$add_argument("--min_cells", help="genes expressed in min cells",default='3')
parser$add_argument("--x_low_cutoff", help="x_low_cutoff for search high variable genes",default='0.0125')
parser$add_argument("--x_high_cutoff", help="x_high_cutoff for search high variable genes",default='3')
parser$add_argument("--y_cutoff", help="y_cutoff for search high variable genes",default='0.5')

args <- parser$parse_args()
str(args)

path=args$path
compare=args$compare
gene_path=args$gene_path
species = args$species
outdir=args$outdir
prefix=args$prefix
Nfeatures=args$Nfeatures
resolution=args$resolution
cca_use=args$cca_use
min_cells=args$min_cells
x_low_cutoff=args$x_low_cutoff
x_high_cutoff=args$x_high_cutoff
y_cutoff=args$y_cutoff

###Nfeatures = args$Nfeatures
min_cells <- as.numeric(min_cells)
resolution <- as.numeric(resolution)
min_cells <- as.numeric(min_cells)
Nfeatures<-as.numeric(Nfeatures)
cca_use <- as.numeric(cca_use)
x_low_cutoff <- as.numeric(x_low_cutoff)
x_high_cutoff <- as.numeric(x_high_cutoff)
y_cutoff <- as.numeric(y_cutoff)

if (is.null(prefix)){
	prefix=compare
}
#creat directory
if (!dir.exists(paste0(outdir,'/Anchors'))){
  dir.create(paste0(outdir,'/Anchors'))
}


if (!dir.exists(paste0(outdir,'/DIFF'))){
  dir.create(paste0(outdir,'/DIFF'))
}

if (!dir.exists(paste0(outdir,'/Marker'))){
  dir.create(paste0(outdir,'/Marker'))
}

if (!dir.exists(paste0(outdir,'/CellsRatio'))){
  dir.create(paste0(outdir,'/CellsRatio'))
}

#if (!dir.exists(paste0(outdir,'/Conserve'))){
#  dir.create(paste0(outdir,'/Conserve'))
# }
#gene annotation
gene.names <- read.table(file = paste0(gene_path, '/',species,'_genes.tsv'),sep = '\t',header = F,stringsAsFactors = FALSE)
gene.names$V2 <- make.unique(gene.names$V2)
colnames(gene.names) <- c("Gene","Gene_name")
myoutdir = outdir
#creat seurat object
ob.list <- list()
samples<-strsplit(compare,'vs')[[1]]

numsap=1
for (each in samples){
  pbmc <- readRDS(paste0(path,'/',each,'_QC.rds'))
	colnames(pbmc@assays$RNA@counts) <- str_replace_all(colnames(pbmc@assays$RNA@counts), '-1',paste0('-',numsap))
	ob <- CreateSeuratObject(counts =pbmc@assays$RNA@counts,project =each,min.cells = min_cells)
	ob$stim <-each
	ob <- NormalizeData(ob)
#  ob <- FindVariableFeatures(ob,  selection.method = "vst",nfeatures = Nfeatures)
	numsap=numsap+1
	ob.list[[each]] <- ob
}
outdir = myoutdir

seurat.obj = merge(x=ob.list[[1]], y=ob.list[[2]])

if(length(ob.list) > 2){
  for (j in 3:length(ob.list)){
    seurat.obj = merge(x=seurat.obj, y=ob.list[[j]])
  }
}


#Nfeatures=length(rownames(x =seurat.obj)) 尽量不要用所有基因跑
seurat.obj <- FindVariableFeatures(seurat.obj,  selection.method = "vst",nfeatures = Nfeatures)

seurat.obj<-ScaleData(seurat.obj,verbose = FALSE)

seurat.obj<-RunPCA(seurat.obj,pc.genes = seurat.obj@var.genes, npcs = cca_use, verbose = FALSE)

sample = compare

# run harmany batch correction
seurat.obj <- RunHarmony(seurat.obj,"stim", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(seurat.obj, 'harmony')

harmony.data<-cbind(rownames(harmony_embeddings),harmony_embeddings)
colnames(harmony.data)[1]<-"Barcode"
write.table(harmony.data,file=paste0(outdir,'/Anchors/',prefix,'_harmony.csv'),row.names=F,col.names=T,sep=',',quote=F)

## run TSNE
seurat.obj = RunTSNE(seurat.obj, dims.use=1:cca_use, do.fast=TRUE, reduction.use="pca")
#seurat.obj = RunTSNE(seurat.obj, dims.use=1:cca_use, do.fast=TRUE, reduction.use="harmony")

# run UMAP
seurat.obj <- RunUMAP(seurat.obj, reduction = "harmony", dims = 1:cca_use,umap.method = "umap-learn",metric = "correlation")
seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony", dims = 1:cca_use)
seurat.obj <- FindClusters(seurat.obj, resolution = resolution)




#saveRDS(seurat.obj,file=paste0(outdir,'/',prefix,'_combined.rds'))

#write the common high_varience_gene
high_varience_gene <- data.frame(VariableFeatures(seurat.obj))
colnames(high_varience_gene)[1] <- 'Gene_name'
high_varience_gene <- left_join(high_varience_gene, gene.names, by="Gene_name")

high_varience_gene <- high_varience_gene[,c(2,1)]
write.table(high_varience_gene,file = paste0(outdir,'/Anchors/',prefix,"_high_varience_gene.csv"),sep=',',quote = F,row.names =F)


umap.data<-seurat.obj@reductions$umap@cell.embeddings
umap.data<-cbind(rownames(umap.data),umap.data)
colnames(umap.data)[1]<-"Barcode"
write.table(umap.data,file=paste0(outdir,'/Anchors/',prefix,'_UMAP.csv'),row.names=F,col.names=T,sep=',',quote=F)


p1 <- DimPlot(seurat.obj, reduction = "umap", group.by = "stim")
p2 <- DimPlot(seurat.obj, reduction = "umap", label = TRUE)
p3 <- DimPlot(seurat.obj, reduction = "umap",split.by = "stim")

pdf(paste0(outdir,'/Anchors/',prefix,'_UMAP.pdf'))
plot_grid(p1, p2)
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_UMAP.png'),type="cairo-png")
plot_grid(p1, p2)
dev.off()

#svg(paste0(outdir,'/Anchors/',prefix,'_UMAP.svg'))
#plot_grid(p1, p2)
#dev.off()


pdf(paste0(outdir,'/Anchors/',prefix,'_sample_UMAP.pdf'))
p1
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_sample_UMAP.png'),type="cairo-png")
p1
dev.off()

#svg(paste0(outdir,'/Anchors/',prefix,'_sample_UMAP.svg'))
#p1
#dev.off()


pdf(paste0(outdir,'/Anchors/',prefix,'_clusters_UMAP.pdf'))
p2
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_clusters_UMAP.png'),type="cairo-png")
p2
dev.off()


pdf(paste0(outdir,'/Anchors/',prefix,'_split_UMAP.pdf'))
p3
dev.off()

#svg(paste0(outdir,'/Anchors/',prefix,'clusters_UMAP.svg'))
#p2
#dev.off()

tsne.data<-seurat.obj@reductions$tsne@cell.embeddings
tsne.data<-cbind(rownames(tsne.data),tsne.data)
colnames(tsne.data)[1]<-"Barcode"
write.table(tsne.data,file=paste0(outdir,'/Anchors/',prefix,'_tSNE.csv'),row.names=F,col.names=T,sep=',',quote=F)

p1 <- DimPlot(seurat.obj, reduction = "tsne", group.by = "stim")
p2 <- DimPlot(seurat.obj, reduction = "tsne", label = TRUE)


pdf(paste0(outdir,'/Anchors/',prefix,'_TSNE.pdf'))
plot_grid(p1, p2)
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_TSNE.png'),type="cairo-png")
plot_grid(p1, p2)
dev.off()

#svg(paste0(outdir,'/Anchors/',prefix,'_TSNE.svg'))
#plot_grid(p1, p2)
#dev.off()


pdf(paste0(outdir,'/Anchors/',prefix,'_sample_TSNE.pdf'))
p1
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_sample_TSNE.png'),type="cairo-png")
p1
dev.off()

#svg(paste0(outdir,'/Anchors/',prefix,'_sample_TSNE.svg'))
#p1
#dev.off()


pdf(paste0(outdir,'/Anchors/',prefix,'_clusters_TSNE.pdf'))
p2
dev.off()

png(paste0(outdir,'/Anchors/',prefix,'_clusters_TSNE.png'),type="cairo-png")
p2
dev.off()



###cluster result
write.table(Idents(seurat.obj),file=paste0(outdir,'/Anchors/',prefix,'_cluster.csv'),row.names=T,col.names=c('Barcode,Cluster'),quote=F,sep=',')



#samples cluster percent
sample_clus <- data.frame(Idents(seurat.obj))
sample_clus <- cbind(rownames(sample_clus),sample_clus)
colnames(sample_clus) <- c('Barcode','Cluster')

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

sample_number <- length(samples)
num <- sample_number + 1

b <- mutate(sample_clus,Sample=as.numeric(str_split(sample_clus$Barcode,'-',simplify = TRUE)[,2]))

# stat counts of each sample in each cluster
c <- b %>% 
  group_by(Cluster) %>% 
  count(Sample) %>% 
  spread(key=Sample,value=n)
colnames(c)[2:num] <- samples

#calculate the relative percent of each sample in each cluster
df <- c[,-1]
rowsum <- rowSums(df,na.rm=T)
df_r <- df/rowsum

#print out the relative table with heads
df_r <- cbind(as.numeric(rownames(df_r))-1, df_r)
colnames(df_r)[1] <- "Cluster"
colnames(df_r)[2:num] <- samples

write.table(c,file=paste0(outdir,'/CellsRatio/',prefix,'_cluster_abundance.csv'),sep=',',row.names = F, quote=F)
write.table(df_r,file=paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.csv'),sep=',',row.names = F, quote=F)

#draw barplot

colour1 <- if(sample_number>40) hue_pal()(sample_number) else allcolour[1:sample_number]
td <- gather(df_r,key="Cluster Name",value="Cells Ratio",-Cluster)
td[,1] <- factor(td[,1], levels = sort(as.numeric(df_r$Cluster)))
td[,2] <- factor(td[,2], levels = samples)


plt<- ggplot(td,aes(x=td[,1],y=td[,3],fill=td[,2]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Cluster Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour1)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Sample'))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45))


pdf(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.pdf'))
plt
dev.off()

pdf(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent_adj.pdf'),height=17,width = 25)
plt
dev.off()

png(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.png'),type="cairo-png")
plt
dev.off()

#svg(paste0(outdir,'/CellsRatio/',prefix,'_cluster_persent.svg'))
#plt
#dev.off()

#clusters in sample percent
df_c <- t(t(df)/rowSums(t(df),na.rm=T))
row.names(df_c) <-c(1:dim(df_c)[1])
df_c <- as.data.frame(cbind((as.numeric(rownames(df_c))-1), df_c))
colnames(df_c)[1] <- "Cluster"
write.table(df_c,file=paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.csv'),sep=',',row.names = F, quote=F)

cluster_number <- length(row.names(df_c))
colour2 <- if(cluster_number>40) hue_pal()(cluster_number) else allcolour[1:cluster_number]
td_c <- gather(df_c,key="Sample Name",value="Cells Ratio",-Cluster)
td_c[,1] <- factor(td_c[,1], levels = sort(as.numeric(df_c$Cluster)))
td_c[,2] <- factor(td_c[,2], levels = samples)

plt_c<- ggplot(td_c,aes(x=td_c[,2],y=td_c[,3],fill=td_c[,1]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="Sample Name",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour2)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45))


pdf(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.pdf'))
plt_c
dev.off()

png(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.png'),type="cairo-png")
plt_c
dev.off()


pdf(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent_adj.pdf'),height=17,width=25)
plt_c
dev.off()

#svg(paste0(outdir,'/CellsRatio/',prefix,'_sample_persent.svg'))
#plt_c
#dev.off()




#Identify conserved cell type markers 
#library(dplyr)
#DefaultAssay(combined) <- "RNA"
#myfilt<-combined@meta.data %>%select(orig.ident,seurat_clusters)%>%group_by(orig.ident) %>%unique() %>% table()%>% as.data.frame() %>% filter(Freq==0)
#ConserveCluster<-setdiff(as.vector(levels(combined@meta.data$seurat_clusters)),as.vector(myfilt$seurat_clusters))
#for(eachC in  ConserveCluster){
#	nk.markers <- FindConservedMarkers(combined, ident.1 = eachC , grouping.var = "stim", verbose = FALSE)
#	mymarkers<-data.frame(nk.markers)
#	mymarkers<-cbind(rownames(mymarkers),mymarkers)
#	colnames(mymarkers)[1]<-"Gene_name"
#	topgene<-(mymarkers%>%top_n(4,max_pval))$Gene_name
#	topgene<-as.vector(topgene)
#
#	write.table(mymarkers,file=paste0(outdir,'/Conserve/','Cluster_',eachC,'_Conserved.xls'),quote=F,row.names=F,col.names=T,sep='\t')
#	write.table(mymarkers,file=paste0(outdir,'/Conserve/','Cluster_',eachC,'_Conserved.csv'),quote=F,row.names=F,col.names=T,sep=',')
#
#	p1<-FeaturePlot(combined, features = topgene, min.cutoff = "q9")
#	pdf(paste0(outdir,'/Conserve/',prefix,'Cluster_',eachC,'_Conserved.pdf'))
#	print(p1)
#	dev.off()
#
#	png(paste0(outdir,'/Conserve/',prefix,'Cluster_',eachC,'_Conserved.png'),type="cairo-png")
#	print(p1)
#	dev.off()
#
#	svg(paste0(outdir,'/Conserve/',prefix,'Cluster_',eachC,'_Conserved.svg'))
#	print(p1)
#	dev.off()

#}





##Find All Markers
seurat.obj.markers <- FindAllMarkers(object = seurat.obj, only.pos = F, min.pct = 0.25)

colnames(seurat.obj.markers)[7] <- 'Gene_name'
seurat.obj.markers <- left_join(seurat.obj.markers, gene.names, by="Gene_name")
seurat.obj.markers <- seurat.obj.markers[,c(8,7,6,2,1,5,3,4)]
write.table(seurat.obj.markers,file=paste0(outdir,'/DIFF/','Cluster_diff.xls'),quote=F,row.names=F,col.names=T,sep='\t')
write.table(seurat.obj.markers,file=paste0(outdir,'/DIFF/','Cluster_diff.csv'),quote=F,row.names=F,col.names=T,sep=',')

cluster <- unique(seurat.obj.markers$cluster)
for (i in cluster) {
  data <- filter(seurat.obj.markers,cluster == i)
  write.table(data,paste0(outdir,'/DIFF/','Cluster_',i,'_diff.xls'),quote=F,row.names=F,col.names=T,sep='\t')
  data1 <- filter(seurat.obj.markers,cluster == i, p_val_adj < 0.05, avg_logFC > 0)
  write.table(data1,paste0(outdir,'/DIFF/','Cluster_',i,'_diff_significant.xls'),quote=F,row.names=F,col.names=T,sep='\t')
}

###Marker Analysis
##top genes 
top_markers <- seurat.obj.markers %>% group_by(cluster) %>% top_n(4, avg_logFC)
colnames(top_markers)[2] <-'gene'
clus<-unique(top_markers$cluster)

# violin and tsne
for (each in clus){
  top_clus<-top_markers[top_markers$cluster==each,]
  genes<-top_clus$gene
  
  p3<-VlnPlot(seurat.obj,genes,ncol =2,pt.size = 0)
  p3<-p3+xlab("Cluster")+ylab("log(UMI)")
  ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_violin.pdf'), p3,height=5,width =6)
#  ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_violin.svg'), p3,height=5,width =6)
  ggsave(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_violin.png'), p3,type = 'cairo-png')
  
  pdf(paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_tsne.pdf'), height=6,width =7)
  print(FeaturePlot(object =seurat.obj,genes,cols = c("grey", "blue"), reduction = "tsne"))
  dev.off()
	mypdf = paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_tsne.pdf')
	mypng  = paste0(outdir,'/Marker/',prefix,'_Cluster_',each,'_tsne.png')
	system(paste("/usr/bin/convert  -density 600",   mypdf ,mypng  ))
  
}

##Heatmap
heatgene = intersect(top_markers$gene,VariableFeatures(seurat.obj))
DoHeatmap(object = seurat.obj,  features= heatgene)
ggsave(paste0(outdir,'/Marker/',prefix,'_Heatmap.pdf'),width = 15, height=10)
ggsave(paste0(outdir,'/Marker/',prefix,'_Heatmap.png'),type="cairo-png",width = 10)


saveRDS(seurat.obj, file=paste0(outdir,'/',prefix,'.rds'))
save.image(file=paste0(outdir,'/',prefix,'.RData'))
sessionInfo()
