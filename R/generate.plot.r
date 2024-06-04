library(gridExtra)
library(ggplot2)
library(Seurat)
library(patchwork)

generate.plot <- function(seurat_object, project, genes, cell_annotation, cell_types, treatment_col, treatment = NULL) {
  output_dir <- paste0(project, "_SC_Expression")
  dir.create(output_dir, showWarnings = FALSE)  # Create directory if it doesn't exist
  
  for (gene in genes) {
    gene_output_dir <- file.path(output_dir, paste(gene))
    dir.create(gene_output_dir, showWarnings = FALSE)  # Create gene-specific directory
    
    for (treat in unique(seurat_object@meta.data[[treatment_col]]) %||% NA) {
      treatment_dir <- ifelse(is.null(treatment), gene_output_dir, file.path(gene_output_dir, treat))
      dir.create(treatment_dir, showWarnings = FALSE)  # Create treatment-specific directory
      
      for (cell_type in cell_types) {
        input_gene <- strsplit(gene, "_")[[1]]
        gene1 <- as.character(input_gene[1])
        gene2 <- as.character(input_gene[2])
        
        # Use tryCatch to handle the subset error
        tryCatch({
          # Compute the expression of the genes of interest
          g1 <- subset(seurat_object, subset = !!sym(cell_annotation) == cell_type & (!!sym(treatment_col) == treat | is.null(treatment)) & !!sym(gene1) > 0 & !!sym(gene2) == 0)
          g1 <- AddModuleScore(g1, features = gene1, name = "module_score")
          g1 <- WhichCells(g1, expression = module_score1 != 0)
          
          g2 <- subset(seurat_object, subset = !!sym(cell_annotation) == cell_type & (!!sym(treatment_col) == treat | is.null(treatment)) & !!sym(gene2) > 0 & !!sym(gene1) == 0)
          g2 <- AddModuleScore(g2, features = gene2, name = "module_score")
          g2 <- WhichCells(g2, expression = module_score1 != 0)
          
          both <- subset(seurat_object, subset = !!sym(cell_annotation) == cell_type & (!!sym(treatment_col) == treat | is.null(treatment)) & !!sym(gene1) > 0 & !!sym(gene2) > 0)
          both <- AddModuleScore(both, features = c(gene1, gene2), name = "module_score")
          both <- WhichCells(both, expression = module_score1 != 0 & module_score2 != 0)
          
          neither <- subset(seurat_object, subset = !!sym(cell_annotation) == cell_type & (!!sym(treatment_col) == treat | is.null(treatment)) & !!sym(gene2) == 0 & !!sym(gene1) == 0)
          
          # Create tabl
          lg1 <- length(g1)
          lg2 <- length(g2)
          bt1 <- length(both)
          non <- length(colnames(neither))
          
          non_p <- round(1 - (bt1 + lg1 + lg2) / non, 4)*100
          bt1_p <- round(bt1 / (lg1 + lg2 + non), 4)*100
          lg1_p <- round(lg1 / (bt1 + lg2 + non), 4)*100
          lg2_p <- round(lg2 / (bt1 + lg1 + non), 4)*100
          
          dtab <- data.frame(
            nCells = c(non, bt1, lg1, lg2),
            pctExpress =  c(non_p, bt1_p, lg1_p, lg2_p),
            row.names = c("Neither", "Both", gene1, gene2)
          )
          
          # Generate plot
          labels <- c("Neither", "Co-expression", gene1, gene2)
          values <- c("grey", "purple", "red", "blue")
          
          # Generate the table plot
          table_plot <- ggplot() +
            theme_void() +  # Remove background and axes
            annotation_custom(tableGrob(dtab, rows = rownames(dtab)), 
                              xmin = 0, xmax = 1, ymin = 0, ymax = 1)

          # Generate the main plot
          main_plot <- suppressMessages({suppressWarnings({
            DimPlot(seurat_object, cells.highlight = list(g1, g2, both), sizes.highlight = 0.2, shuffle = TRUE) + 
              scale_color_manual(labels = labels, values = values) +
              labs(color = cell_type) +
              theme(plot.width = unit(12, "in"), plot.height = unit(6, "in"), legend.position = c(1.4, 0.9))
          })})

          # Combine the main plot and table plot using Patchwork
        plot <- main_plot + table_plot  + plot_annotation(
                    title = paste("Co-expression of", gsub("_", " and ", gene)),
                    subtitle = treatment,
                    caption = project
                )

          sanitized_cell_type <- gsub("[^A-Za-z0-9_]", "_", cell_type)
          filename <- ifelse(is.null(treatment), file.path(treatment_dir, paste0(sanitized_cell_type, "_", gene, ".pdf")), 
                            file.path(treatment_dir, paste0(treat, "_", sanitized_cell_type, "_", gene, ".pdf")))
          
          ggsave(filename, plot = plot, width = 8, height = 6, units = "in", dpi = 300)
          
        }, error = function(e) {
          # If no cells are found, generate a blank PDF
          blank_plot <- ggplot() + theme_void()
          filename <- ifelse(is.null(treatment), file.path(treatment_dir, paste0("NONE_FOR_", gsub("[^A-Za-z0-9_]", "_", cell_type), "_", gene, ".pdf")), 
                            file.path(treatment_dir, paste0("NONE_FOR_", treat, "_", gsub("[^A-Za-z0-9_]", "_", cell_type), "_", gene, ".pdf")))
          ggsave(filename, plot = blank_plot, width = 8, height = 6, units = "in", dpi = 300)
        })
      }
    }
  }
}