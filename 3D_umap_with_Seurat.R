

# Load plot_ly
library(Seurat)
library(plotly)
library(RColorBrewer)
# Construct a dataframe using data from your pre-clustered Seurat v3.1.1 object
# Here 'seurat_clusters' is list of numeric cluster identities, you can find it here: yourseuratobject[["seurat_cluster"]], 
# or yourseuratobject$seurat_clusters, where 'yourseuratobject' is a Seurat object created with Seurat v3.1.1 (works for v3.0.0 as well)

yourseuratobject <- readRDS('./data/pbmc3k.rds')

# Re-run UMAPs that you have accurate calculations for all UMAP(s)
yourseuratobject <- RunUMAP(yourseuratobject,
                            dims = 1:10,
                            n.components = 3L)

# Extract tSNE information from Seurat Object
umap_1 <- yourseuratobject[["umap"]]@cell.embeddings[,1]
umap_2 <- yourseuratobject[["umap"]]@cell.embeddings[,2]
umap_3 <- yourseuratobject[["umap"]]@cell.embeddings[,3]

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = yourseuratobject, reduction = "umap")

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# get number of clusters
cluster_number = length(levels(plot.data$seurat_clusters))
sprintf('cluster_number is %s', cluster_number)
all_colors <- c(
    # "#000000",  # remove the black, as often, we have black colored annotation
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
    "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
    "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
    "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
    "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
    "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
    "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
    "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
    "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
    "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72"
)
manual_colors <- all_colors[0:cluster_number]
# Plot your data, in this example my Seurat object had 21 clusters (0-20)
p <- plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~seurat_clusters, 
        colors = manual_colors,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

widget_file_size <- function(p) {
  d <- getwd()
  withr::with_dir(d, htmlwidgets::saveWidget(p, "index.html"))
  f <- file.path(d, "index.html")
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
  sprintf('save plot to %s', paste(d,'index.html',sep='/'))
}
widget_file_size(partial_bundle(p))
# Say you wanto make a gene-expression 3D plot, where you can plot gene expression against a color scale
# Here using the same seurat object as above, we extract gene expression information for beta-actin 'ACTB'
# Here we concentrate on SCT normalized data, or log normalized RNA NOT raw counts.
# In addition if you want, you may look at normalised-RNA, SCT or integrated slots, to look at gene expression
# Setting your DefaultAssay() will inform R which assay to pick up expression data from.
# DefaultAssay(object = yourseuratobject)
# DefaultAssay(object = yourseuratobject) <- "RNA"
# DefaultAssay(object = yourseuratobject) <- "integrated"
# DefaultAssay(object = yourseuratobject) <- "SCT"

# create a dataframe
# plot.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "ACTB"), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'changed' column of your dataframe
# plot.data$changed <- ifelse(test = plot.data$ACTB <1, yes = plot.data$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
# plot.data$label <- paste(rownames(plot.data)," - ", plot.data$ACTB, sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
# plot_ly(data = plot.data, 
#         x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
#         color = ~changed, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
#         opacity = .5,
#         colors = c('darkgreen', 'red'), 
#         type = "scatter3d", 
#         mode = "markers",
#         marker = list(size = 5, width=2), 
#         text=~label,
#         hoverinfo="text"
# )

# On running this code the HTML output should appear in RStudio. You can save the output as a
# HTML file. Once you have saved, just open the HTML file in any web browser (double click on the html- file
# and if asked select to open with any web browser like google chrome/safari/mozilla/explorer etc).
# It should be have all of the integrated features you saw in the RStudio output file.

########## #
########## #

# Alternative method as designed by @vertesy (Thanks for the suggestions!)
# create a dataframe
# goi <- "TOP2A"
# plotting.data <- FetchData(object = yourseuratobject, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=goi), slot = 'data')

# Say you want change the scale, so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'Expression' column of your dataframe
# Cutoff <- 2
# Cutoff <- quantile(plotting.data[,goi], probs = .95)
# plotting.data$"ExprCutoff" <- ifelse(test = plotting.data[,goi] <Cutoff, yes = plotting.data[,goi], no = Cutoff)

# Add the label column, so that now the column has 'cellname-its expression value'
# plotting.data$label <- paste(rownames(plotting.data)," - ", plotting.data[,goi], sep="")

# Plot your data, in this example my Seurat object had 21 clusters (0-20), and cells express a gene called ACTB
# plot_ly(data = plotting.data,
#         # name = goi,
#         x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
#         color = ~ExprCutoff, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
#         opacity = .5,
#         colors = c('darkgrey', 'red'), 
#         type = "scatter3d", 
#         mode = "markers",
#         marker = list(size = 1), 
#         text=~label,
#         hoverinfo="text"
# ) %>%layout(title=goi)

# Thank you for reading and using this code to further your scRNAseq analysis!
# If you liked it, dont forget to acknowledge, fork and star!
# Citation information is within the Readme, please dont forget to cite!
# Have a wonderful day!!