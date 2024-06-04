make.merge.pdfs <- function(dir, genes, study_id, study_title) {

  library(pdftools)
  
  unlink(grep("pdf$", list.files("/tmp", full.names = TRUE), value = TRUE))

  study_id <- gsub("[^a-zA-Z0-9]", "", study_id)
  study_title <- gsub("[^a-zA-Z0-9]", " ", study_title)

  if (length(genes) == 1) {
    genes <- genes[[1]]
    genes_sym <- gsub("_", " and ", genes)

    # Generate title pdf
    txt <- paste(study_id, study_title, paste("Co-expression of", genes_sym), sep = "\n")
    page <- paste0("/tmp/", genes, ".pdf")
    pdf(page, width=8.27, height=11.69) # Default ggsave page size in inches (A4 format)
    par(mar=c(0.1, 0.1, 0.1, 0.1))
    plot.new()
    text(0.01, 0.98, txt, font=4, cex=1.5, col="#000000", xpd=TRUE, adj=c(0, 1))
    dev.off()

    subdirs <- list.dirs(path = file.path(dir, genes), full.names = TRUE, recursive = FALSE)
    subdir_name <- basename(subdirs)
    subdir_files <- list.files(path = subdirs, pattern = "\\.pdf$", full.names = TRUE)
    valid_files <- subdir_files[!grepl("/NONE_", subdir_files)]
    pdf_combine(c(page, valid_files), output = file.path(dir, basename(page)))

  } else {
    for (gene in genes) {
      genes_sym <- gsub("_", " and ", gene)

      # Generate title pdf
      txt <- paste(study_id, study_title, paste("Co-expression of", genes_sym), sep = "\n")
      page <- paste0("/tmp/", gsub(" ", "_", gene), ".pdf")
      pdf(page, width=8.27, height=11.69) # Default ggsave page size in inches (A4 format)
      par(mar=c(0.1, 0.1, 0.1, 0.1))
      plot.new()
      text(0.01, 0.98, txt, font=4, cex=1.5, col="#000000", xpd=TRUE, adj=c(0, 1))
      dev.off()

      subdirs <- list.dirs(path = file.path(dir, gene), full.names = TRUE, recursive = FALSE)
      subdir_name <- basename(subdirs)
      subdir_files <- list.files(path = subdirs, pattern = "\\.pdf$", full.names = TRUE)
      valid_files <- subdir_files[!grepl("/NONE_", subdir_files)]
      pdf_combine(c(page, valid_files), output = file.path(dir, basename(page)))
    }
  }
}