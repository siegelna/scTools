suppressPackageStartupMessages({
    library(Seurat)
    library(stringr)
    library(data.table)
})

express.meta.analysis <- function(index, seurat_obj, treatment_col, treatments, genes, output = "percentage", cell_types = NULL) {

  processed_dfs_list <- list() 
  
  for (gene in genes) {
    input_gene <- strsplit(gene, "_")[[1]]
    gene1 <- as.character(input_gene[1])
    gene2 <- as.character(input_gene[2])
    
    cells <- unique(seurat_obj@meta.data[['monaco.fine']])
    cells <- cells[!is.na(cells)]
    
    result_df <- data.frame(matrix(ncol = length(cells), nrow = length(treatments)))
    colnames(result_df) <- paste(cells, gene, sep = "_")
    
    row_index <- 1
    
    suppressWarnings({
      for (treatment in treatments) {
        for (cell in cells) {
          tryCatch({
            if (!str_detect(gene, "_")) {
              expr_test <- tryCatch({
                subset(x = seurat_obj, monaco.fine == cell & !!sym(treatment_col) == treatment & !!sym(gene1) > 0)
              }, error = function(e) {
                data.frame()
              })
            } else {
              expr_test <- tryCatch({
                subset(x = seurat_obj, monaco.fine == cell & !!sym(treatment_col) == treatment & !!sym(gene1) > 0 & !!sym(gene2) > 0)
              }, error = function(e) {
                data.frame()
              })
            }

            all_test <-
             subset(x = seurat_obj, subset = monaco.fine == cell & !!sym(treatment_col) == treatment)
            
            if (ncol(all_test) > 0) {
              if(output == "percentage") {
                percentage <- (ncol(expr_test) / ncol(all_test)) * 100
                result <- round(percentage, 2)
              } else if(output == "ratio") {
                result <- paste0(ncol(expr_test), "/", ncol(all_test))
              } else {
                warning("Invalid output type. Defaulting to percentage.")
                percentage <- (ncol(expr_test) / ncol(all_test)) * 100
                result <- round(percentage, 2)
              }
            } else {
              result <- 0
            }
            
            result_df[row_index, c("Accession", "Cohort", "Seurat_Group", "Seurat_Group_Name", paste(cell, gene, sep = "_"))] <-
              c(accession[[index]], cohort[[index]], treatment_col, treatment, result)
            
            row_index <- row_index + 1
          }, error = function(e) {
            result_df[row_index, c("Accession", "Cohort", "Seurat_Group", "Seurat_Group_Name", paste(cell, gene, sep = "_"))] <-
              c(accession[[index]], cohort[[index]], treatment_col, treatment, 0)
            
            row_index <- row_index + 1
          })
        }
      }
    })
    
    result_df[is.na(result_df)] <- 0
    processed_dfs_list[[gene]] <- result_df
  }
  
  study_md <- select(processed_dfs_list[[1]], Accession, Cohort, Seurat_Group, Seurat_Group_Name)
  study_md[study_md == "0"] <- 0
  study_md <- study_md[rowSums(study_md == 0, na.rm = TRUE) < ncol(study_md), ]
  study_md <- unique(study_md)
  
  desired_columns <- c("Seurat_Group", "Seurat_Group_Name")
  
  removed_dfs <- processed_dfs_list[!sapply(processed_dfs_list, function(df) all(desired_columns %in% colnames(df)))]
  
  for (i in seq_along(removed_dfs)) {
    colnames(removed_dfs[[i]]) <- paste(colnames(removed_dfs[[i]]), ifelse(output == "percentage", "%", " Ratio"))
  }
  
  removed_dfs <- lapply(removed_dfs, function(df) {
    for (col in desired_columns) {
      if (!(col %in% colnames(df))) {
        df[[col]] <- NA
      }
    }
    return(df)
  })
  
  subset_dfs <- processed_dfs_list[sapply(processed_dfs_list, function(df) all(desired_columns %in% colnames(df)))]
  
  subset_dfs <- lapply(subset_dfs, function(x) {
    x %>%
      group_by(Seurat_Group, Seurat_Group_Name) %>%
      summarize(across(everything(), list(max)))
  })
  
  df <- append(subset_dfs, removed_dfs)
  
  # Merge cell dataframes
  combined_df <- do.call(cbind, df)
  
  original_columns <- grep("^Seurat_Group|^Seurat_Group_Name|^Accession|^Cohort", names(combined_df), value = TRUE)

  # Select cell-related columns
  if (!is.null(cell_types)) {
    data <- combined_df %>% select(!all_of(original_columns)) %>%
      select(grepl(paste(cell_types, collapse = "|"), colnames(.)))
  } else {
    data <- combined_df %>% select(!all_of(original_columns))
  }

  # Bind data with metadata
  final <- cbind(study_md, data)

  # Append percentage
  colnames(final) <- colnames(final) %>% 
  str_replace(ifelse(output == "percentage", "_1", "_1"), ifelse(output == "percentage", " %", " Ratio"))

  # Convert to numeric if output is percentage
  if(output == "percentage") {
    final[, -(1:4)] <- sapply(final[, -(1:4)], as.numeric)
  }
  
  return(final)
}