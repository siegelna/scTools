remove_low_count_cells <- function(seurat_obj, metadata_column, threshold = 20) {

  total_counts <- table(seurat_obj@meta.data[[metadata_column]])

  low_count_cells <- names(total_counts[total_counts < threshold])

  seurat_obj <- subset(seurat_obj, idents = low_count_cells, invert = TRUE)

  return(seurat_obj)
}

# Remove cells from Seurat object if their total sum is less than a specified amount.
#
#   seurat_obj = <any seurat object that has annotated cell types>
#   metadata_column = <The metadata column that contains the cells to be removed>
#   threshold = <The min number of cells to keep for each cell type. Default 20 cells>