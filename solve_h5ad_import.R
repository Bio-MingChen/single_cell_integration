library(Seurat)
library(SeuratDisk)

custom_as.sparse.H5Group <- function(x, ...) {
  CheckDots(...)
  for (i in c('data', 'indices', 'indptr')) {
    if (!x$exists(name = i) || !is(object = x[[i]], class2 = 'H5D')) {
      stop("Invalid H5Group specification for a sparse matrix, missing dataset ", i)
    }
  }
  if ('h5sparse_shape' %in% hdf5r::h5attr_names(x = x)) {
    return(Matrix::sparseMatrix(
      i = x[['indices']][] + 1,
      p = x[['indptr']][],
      x = x[['data']][]
      #dims = rev(x = hdf5r::h5attr(x = file[['raw.X']], which = 'h5sparse_shape'))
    ))
  }
  return(Matrix::sparseMatrix(
    i = x[['indices']][] + 1,
    p = x[['indptr']][],
    x = x[['data']][]
  ))
}
environment(custom_as.sparse.H5Group) <- asNamespace('Seurat')
assignInNamespace("as.sparse.H5Group", custom_as.sparse.H5Group, ns = "Seurat")
getAnywhere(as.sparse.H5Group)

pbmc3k <- ReadH5AD('./pbmc3k.h5ad')
print(pbmc3k)