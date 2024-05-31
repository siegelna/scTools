read_seurat_objects_from_files <- function(file_paths, default_assay = "RNA") {
  seurat_objects <- list()
  for (file_path in file_paths) {
    # Extract file extension
    file_extension <- tolower(tools::file_ext(file_path))
    
    # Read Seurat object based on file extension
    R_Data <- c("Rda", "rda", "rdata", "RData")
    if (file_extension %in% R_Data) {  # Add more variations here if needed
      # Read Seurat object from RDA file
      seurat_obj <- get(load(file_path))
      DefaultAssay(object = seurat_obj) <- default_assay  # Set default assay
      seurat_objects[[length(seurat_objects) + 1]] <- seurat_obj
    } else if (file_extension == "h5seurat") {
      # Read Seurat object from h5Seurat file
      seurat_obj <- LoadH5Seurat(file_path, assays = default_assay)
      seurat_objects[[length(seurat_objects) + 1]] <- seurat_obj
    } else {
      # File extension not supported
      warning(paste("Unsupported file extension:", file_extension))
    }
  }
  return(seurat_objects)
}
