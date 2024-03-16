library(sceasy)
library(reticulate)
#use_virtualenv("transfer_env")
#use_python("/TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/sceasy_env/bin/python")
use_python("/TJPROJ6/SC/personal_dir/chenming/software/Miniconda/miniconda202012/envs/transfer_env/bin/python")
library(optparse)
# I changed the part of source code of sceasy for functions of anndata2seurat and seurat2anndata
# anndata2seurat now will check the adata.layers["counts"] and put it as counts if exists
# seurat2anndata will add scale.data to adata.uns['seurat_scale.data'] if it exists
source("/TJPROJ6/SC/personal_dir/chenming/research/scanpy_pipeline/scripts/sceasy_adjustment.R")

option_list <- list(
    make_option(c('-i','--infile'),action='store',default=NA,type='character',
        help='file to convert choose from file suffixed with [.h5ad,.rds]'),
    make_option(c('-o','--ofile'),action='store',default=NA,type='character',
        help='converted output file choose from file suffixed with [.h5ad,.rds]'),
    make_option(c('-c','--config'),action='store',default=NA,type='character',
        help='a file saved one or more files to transfer and their outputfile'),  
    make_option(c('-m','--main_layer'),default="data",type='character',
        help='which matrix to convert,choose from [data,counts,scale.data],default is data,if
        you run h5ad to rds,main_layer indicates which type of data is in your adata.X; if transfering
        rds to h5ad, main_layer designates which type of data to put to adata.X'),
    make_option(c('-a','--add_scale_to_uns'),action="store_true",default=FALSE,
        help='whether or not add scale.data to anndata.uns'),
    make_option(c('-s','--assay'),default="RNA",type="character",
        help="which assay to use for seurat to anndata,default is RNA")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)
allow_suffix <- c("h5ad","rds")
data_type <- c('anndata','seurat')
main_layer <- opt$main_layer
add_scale_to_uns <- opt$add_scale_to_uns
easy_convert <- function(infile,ofile) {
    # supported data format
    # "seurat" for seurat object "anndata" for .h5ad file 
    # "sce" for .rds file and "loom" for .loom file
    if (!is.na(infile) && !is.na(ofile)) {
        infile_vector <- unlist(strsplit(infile,".",fixed=TRUE)) # fixed=True to shut down regular expression
        ofile_vector  <- unlist(strsplit(ofile,".",fixed=TRUE))
        from = data_type[infile_vector[length(infile_vector)] == allow_suffix]
        to = data_type[ofile_vector[length(ofile_vector)] == allow_suffix]
        print(from)
        print(to)
        if ((from == 'seurat') && (to == "anndata")) {
            infile <- readRDS(infile)
            layers <- c("counts","data", "scale.data")
            if (opt$assay == "RNA") {
                transfer_layers <- layers[layers %in% slotNames(infile@assays$RNA)]
            } else if (opt$assay == "integrated") {
                transfer_layers <- layers[layers %in% slotNames(infile@assays$integrated)]
            }
            print(sprintf("transfer_layers is %s",transfer_layers))
            if ( (length(from) == 1) && (length(to) == 1) && (from != to) ) {
            seurat2anndata(infile,outFile=ofile,assay=opt$assay,main_layer=main_layer,transfer_layers=transfer_layers,add_scale_to_uns=add_scale_to_uns)
            # seurat2anndata(infile,from=from,to=to,outFile=ofile,main_layer=main_layer,transfer_layers=transfer_layers)
            sprintf("Converted %s to %s,output file to %s",from,to,ofile)
            }
        } else if(from == 'anndata' && to == "seurat") {
            if ( (length(from) == 1) && (length(to) == 1) && (from != to) ) {
                anndata2seurat(infile,outFile=ofile,main_layer=main_layer)
                # anndata2seurat(infile,from=from,to=to,outFile=ofile,main_layer=main_layer)
                sprintf("Converted %s to %s,output file to %s",from,to,ofile)
            } 
        }

    }
}

infile <- opt$infile
ofile <- opt$ofile
easy_convert(infile,ofile)
print('Convertion completed!')

config <- opt$config

# config format
# sample\trds_path\th5ad_path

if (!is.na(config)) {
    config_df <- read.table(config,sep='\t',header=TRUE,stringsAsFactors=FALSE)
    for (i in 1:nrow(config_df)) {
        infile <- config_df$rds_path[i]
        ofile <- config_df$h5ad_path[i]
        print(infile)
        print(ofile)
        easy_convert(infile,ofile)
    }
    print('Convertion completed!')
}
