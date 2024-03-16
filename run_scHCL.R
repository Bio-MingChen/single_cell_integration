library(optparse)
library(scHCL)
library(DESeq2)

option_list <- list(
    make_option(c('-i','--infile'),help='infile name'),
    make_option(c('-o','--ofile'),default='normalized_counts.txt',help='ofile name'),
    make_option(c('-n','--number'),default=3,type="integer",
        help='top n related cells will use default is 3'),
    make_option(c('-s','--samples_file'),help='a file contain samples information')
)

args <- parse_args(OptionParser(option_list=option_list))
print(args)
#expression matrix
emtx <- as.matrix(read.table(args$infile,header=TRUE,sep="\t",row.names=1,check.names=FALSE,
    stringsAsFactors=FALSE))
# remove columns
emtx
barcodes <-colnames(emtx)
samples_df <- read.table(args$samples_file,sep=',',header=TRUE,stringsAsFactors=FALSE)
meta <- samples_df[samples_df$Barcode %in% barcodes,]
rownames(meta) <- meta$Barcode
meta <- meta[barcodes,]
dds <- DESeqDataSetFromMatrix(countData = emtx, colData = meta, design = ~samplenames)
dds <- estimateSizeFactors(dds,type = 'iterate')
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file=args$ofile, sep="\t", quote=F, col.names=NA)
print(summary(normalized_counts))
