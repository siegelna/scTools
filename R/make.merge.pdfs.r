

make.merge.pdfs <- function(dir, genes, study_id, study_title) {

    library(pdftools)

    study_id <- gsub("[^a-zA-Z0-9]", "", study_id)
    study_title <- gsub("[^a-zA-Z0-9]", " ", study_title)
    
    txt <- paste(study_id, study_title,
    paste("Co-expression of", gsub("_", " and ", genes)),  sep = "\n")

    pdf("/tmp/page.pdf", paper="Letter");plot.new();text(0.5, 0.8, txt, font=4, cex=1.5, col="#000000");dev.off()

    subdirs <- list.dirs(path = file.path(dir, genes), full.names = TRUE, recursive = FALSE)
    subdir_name <- basename(subdirs)
    subdir_files <- list.files(path = subdirs, pattern = "\\.pdf$", full.names = TRUE)
    valid_files <- subdir_files[!grepl("/NONE_", subdir_files)]
    pdf_combine(c("/tmp/page.pdf", valid_files), output = "joined.pdf")

}