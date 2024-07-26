frequency_table <- function(seurat_object, meta.col, group, gene_name, filter_column, ident) {
  # Set Idents of Seurat objects
  Idents(seurat_object) <- ident
  
  # Initialize an empty list to store dataframes
  dfs <- list()
  
  # Loop through gene_name and filter_column
  for (i in seq_along(gene_name)) {
    # Extract gene-specific details
    gene <- gene_name[i]
    filter_col <- filter_column[i]
    
    # Loop through group to create dataframes
    for (grp in group) {
      # Filter data based on group
      filtered_data <- subset(seurat_object, !!sym(meta.col) == grp)
      
      # Calculate total number of cells for this group and each Ident
      total_cells <- table(Idents(filtered_data))
      
      # If filter_col is provided, apply additional filtering
      if (!is.null(filter_col)) {
        filtered_data <- subset(filtered_data, !!sym(filter_col) > 0)
      }
      
      # Create dataframe with frequencies of cells expressing the gene
      df <- as.data.frame(table(Idents(filtered_data)))
      
      # Add Group, Gene, and Total_Cells columns
      df$Group <- grp
      df$Gene <- gene
      
      # Match and add Total_Cells column
      df$Total_Cells <- total_cells[match(df$Var1, names(total_cells))]
      
      # Store dataframe in the list
      dfs[[paste(gene, grp, sep = "_")]] <- df
    }
  }
  
  # Combine all dataframes into one
  df_combined <- bind_rows(dfs) %>%
    select(Gene, Var1, Freq, Total_Cells, Group)
  
  # Arrange the dataframe by Group and Var1 (optional)
  df_combined <- df_combined %>%
    arrange(Group, Var1)
  
  # Return the combined dataframe
  return(df_combined)
}
