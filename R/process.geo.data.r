suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
})

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

library(tools)

process.GEO.data <- function(GEOID) {
  gse <- getGEO(GEOID)
  
  # Get the sample metadata
  sample_metadata <- pData(phenoData(gse[[1]]))
  
  study <- getGEOSuppFiles(GEOID)
  
  # Untar the file
  untar(rownames(study), exdir = GEOID)
  
  # Check if files are already extracted
  file_paths <- list.files(path = GEOID, full.names = TRUE, pattern = "genes|barcodes|matrix", recursive = TRUE, ignore.case = TRUE)
  
  if (length(file_paths) == 0) {
    # Extract files from tarballs if file_paths is empty
    tar_files <- list.files(path = GEOID, pattern = "\\.tar\\.gz$", full.names = TRUE)
    
    for (tarball_file in tar_files) {
      output_dir <- file_path_sans_ext(basename(tarball_file))
      output_dir <- gsub("\\.tar$", "", output_dir)
      untar(tarball_file, exdir = output_dir)
      
      files <- list.files(path = output_dir, full.names = TRUE, recursive = TRUE)
      
      for (file in files) {
        directory_name <- basename(output_dir)
        file_name <- basename(file)
        directory_prefix <- sub("^(.*?)_.*", "\\1", directory_name)
        new_file_name <- paste(directory_prefix, file_name, sep = "_")
        file.rename(file, file.path(GEOID, new_file_name))
      }

      # Remove the original .tar.gz file
      file.remove(tarball_file)
      
      unlink(output_dir, recursive = TRUE)
    }
    
    # Update file_paths after extraction
    file_paths <- list.files(path = GEOID, full.names = TRUE, pattern = "features|genes|barcodes|matrix", recursive = TRUE, ignore.case = TRUE)
  }
  
  # Proceed with the rest of the function logic if files are found
  if (length(file_paths) > 0) {
    # Proceed with the rest of the function logic
    
    # List all files in the directory with full path, excluding files with extension ".h5"
    file_paths <- file_paths[!grepl("\\.h5$", file_paths)]
    
    # Extract GSM numbers from the file paths
    gsm_numbers <- gsub(".*GSM([0-9]+)_.*", "\\1", file_paths)
    
    # Count occurrences of each GSM number
    gsm_counts <- table(gsm_numbers)
    
    # Find GSM numbers that don't have all three files
    incomplete_gsms <- names(gsm_counts)[gsm_counts != 3]
    
    # Print the removed file paths
    cat("Removed file paths:\n")
    removed_paths <- character(0)  # Initialize empty vector to store removed paths
    for (gsm in incomplete_gsms) {
      paths_to_remove <- file_paths[grep(gsm, file_paths)]
      cat(paths_to_remove, sep = "\n", "\n")
      removed_paths <- c(removed_paths, paths_to_remove)  # Append paths to removed_paths vector
    }
    
    # Remove incomplete paths from file_paths
    file_paths <- setdiff(file_paths, removed_paths)
    
    # Initialize a list to store the Seurat objects
    seurat_objects <- list()
    
    # Loop over the files
    for (i in seq(1, length(file_paths), by = 3)) {
      # Extract sample identifier
      sample_id <- sub("(barcodes|features|genes|matrix)\\.tsv\\.gz", "", basename(file_paths[i]))
      sample_id <- gsub("[^A-Za-z0-9]$", "", sample_id)
      
      # Read the files for the current sample
      result <- ReadMtx(
        mtx = file_paths[i + 2],
        cells = file_paths[i],
        features = file_paths[i + 1]
      )
      
      result@Dimnames[[1]] <- make.unique(result@Dimnames[[1]])
      result@Dimnames[[2]] <- make.unique(result@Dimnames[[2]])
      
      # Create Seurat object from the result
      seurat_obj <- CreateSeuratObject(
        counts = result,
        assay = "RNA",
        min.cells = 3,
        min.features = 300
      )
      
      # Assign the Seurat object to a variable based on the sample identifier
      assign(sample_id, seurat_obj, envir = .GlobalEnv)
      
      # Store the Seurat object in the list of Seurat objects
      seurat_objects[[sample_id]] <- seurat_obj
      
      # Remove the variable created with ReadMtx
      rm(result)
    }
    
    # Loop over the Seurat objects to add Idents based on sample metadata
    for (sample_id in names(seurat_objects)) {
      # Extract the accession number from the Seurat object name
      geo_accession <- regmatches(sample_id, regexpr("GSM[0-9]{7}", sample_id))
      
      # Check if the accession number is present in the sample metadata
      if (geo_accession %in% sample_metadata$geo_accession) {
        # Get the row index in the sample metadata dataframe
        row_index <- which(sample_metadata$geo_accession == geo_accession)
        
        # Get the column names and corresponding row values
        col_names <- colnames(sample_metadata)
        row_values <- as.character(sample_metadata[row_index, ])
        names(row_values) <- col_names
        
        # Set the Idents in the Seurat object based on the column names and row values
        for (col_name in col_names) {
          if (col_name %in% names(row_values)) {
            seurat_objects[[sample_id]]@meta.data[[col_name]] <- row_values[[col_name]]
          }
        }
      }
    }
    
    # Create a directory called "objects" if it doesn't already exist
    if (!file.exists("objects")) {
      dir.create("objects")
    }
    
    # Save the seurat_objects as a .rda file in the "objects" directory with the GEOID as the filename
    save(seurat_objects, file = file.path("objects", paste0(GEOID, ".rda")))
    
    return(seurat_objects)
  } else {
    # If file_paths is empty, return NULL or an appropriate message
    return(NULL)
  }
}