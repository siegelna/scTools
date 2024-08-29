library(Seurat)
library(openxlsx)
library(dplyr)

# Function to clean and truncate sheet names if they are too long
clean_sheetname <- function(name, max_length = 31) {
    if (nchar(name) > max_length) {
        name <- trimws(name)
        name <- gsub("[^A-Za-z0-9]", "", name)
    }
    return(name)
}

# Function to replace an existing sheet or create a new one
replace_or_add_sheet <- function(wb, sheet_name) {
  # Check if the sheet already exists
  if (sheet_name %in% names(wb)) {
    # Remove the existing sheet
    removeWorksheet(wb, sheet = sheet_name)
  }
  # Add the new sheet
  addWorksheet(wb, sheetName = sheet_name)
}

# Function to perform differential expression analysis and save results to Excel
perform_DE_analysis <- function(seurat_obj, broad_types, condition_col, split_by, output_dir) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop over each broad immune cell type
  for (broad_type in broad_types) {
    # Create a new workbook for each broad type
    wb <- createWorkbook()
    
    # Subset the object for the current broad type
    obj_subset <- subset(seurat_obj, broad.type == broad_type)
    
    # Get unique fine types for the current broad type
    fine_types <- unique(obj_subset@meta.data$fine.type)
    
    # Split the object by fine type
    SeuratObject_splitlist <- SplitObject(obj_subset, split.by = "fine.type")
    
    # Set identity to condition for each subset
    SeuratObject_splitlist <- lapply(SeuratObject_splitlist, function(x) {
      SetIdent(x, value = x@meta.data[[condition_col]])
    })
    
    # Loop over each fine type
    for (fine_type in fine_types) {
      # Prepare and clean the sheet name
      sheetname <- gsub("^\\s*\\d+\\.\\s*", "", fine_type)
      sheetname <- clean_sheetname(sheetname)
    #   sheetname <- make.names(sheetname, unique = TRUE) # Ensure unique and valid sheet names

      if (fine_type == "7. IGHG+ IGL- Plasma") sheetname <- paste("7. IGHG+ IGL- Plasma", "2")
      
      # Replace or add the sheet for each fine type within the broad type
      replace_or_add_sheet(wb, sheetname)
      
      # Run FindMarkers for the current fine type and handle errors
      tryCatch({
        cell_markers <- FindMarkers(SeuratObject_splitlist[[fine_type]], ident.1 = 'LN', ident.2 = 'Control', min.pct = 0.25, logfc.threshold = 0.25)
        
        # Filter significant results
        # cell_markers <- cell_markers %>% filter(p_val_adj <= 0.01)
        
        if (nrow(cell_markers) == 0) {
          cell_markers <- data.frame(Message = "No significant markers found")
        }
        
        # Write the data to the corresponding sheet
        writeData(wb, sheet = sheetname, x = cell_markers, rowNames = TRUE)
        
      }, error = function(e) {
        # Handle errors by writing the error message
        error_message <- paste("Error:", e$message)
        writeData(wb, sheet = sheetname, x = data.frame(Message = error_message), rowNames = FALSE)
      })
    }
    
    # Save the workbook to a file named after the broad type
    filename <- paste("AMP_II_SLE", gsub(" ", "_", broad_type), "Markers.xlsx", sep = "_")
    saveWorkbook(wb, file = file.path(output_dir, filename), overwrite = TRUE)
  }
}

# # Example usage:
# # Set parameters
# seurat_obj <- obj # Your Seurat object
# broad_types <- c('T Cell', 'NK Cell', 'Plasma Cell', 'B Cell', 'Myeloid Cell')
# condition_col <- 'Type' # Metadata column representing conditions
# split_by <- 'Type' # Metadata column to split by
# output_dir <- "Markers/Test"

# # Call the function
# perform_DE_analysis(seurat_obj, broad_types, condition_col, split_by, output_dir)
