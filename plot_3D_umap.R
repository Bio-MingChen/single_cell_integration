library(optparse)
library(plotly)

option_list <- list(
    make_option(c('-i','--infile'),action='store',default=NA,type='character',
    help='infile with umap_1 umap_2 umap_3 cluster label as title'),
    make_option(c('-o','--ofile'),action='store',default="index.html",type='character',
    help='ofile name with .html as suffix,default is index.html')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

infile <- opt$infile
ofile <- opt$ofile

plot.data <- read.table(infile,sep='\t',header=TRUE,stringsAsFactors=FALSE)

# get number of clusters
plot.data$cluster <- factor(plot.data$cluster)
cluster_number = length(levels(plot.data$cluster))

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
        color = ~cluster, 
        colors = manual_colors,
        type = "scatter3d", 
        mode = "markers", 
        marker = list(size = 5, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="cluster") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

widget_file_size <- function(p) {
  d <- getwd()
  withr::with_dir(d, htmlwidgets::saveWidget(p, ofile))
  f <- file.path(d, ofile)
  mb <- round(file.info(f)$size / 1e6, 3)
  message("File is: ", mb," MB")
  sprintf('save plot to %s', paste(d,ofile,sep='/'))
}
widget_file_size(partial_bundle(p))