# scTools

## Dependencies:
```yaml
channels:
 - conda-forge
 - bioconda
dependencies:
 - r-seurat
 - r-stringr
 - r-data.table
 - r-gridextra
 - r-ggplot2
 - r-patchwork
 - r-pdftools
 - bioconductor-geoquery
```
## Functions

### generate.plot
```r
obj <- 1
treatment_col <- "stim"
cell_annotation <- "monaco.main"
generate.plot(
    seurat_object = seurat_objects[[obj]], 
    project = accession[[obj]], 
    genes = genes, 
    cell_annotation = cell_annotation,
    cell_types = unique(seurat_objects[[obj]]@meta.data[[cell_annotation]]), 
    treatment_col = treatment_col, 
    treatment = unique(seurat_objects[[obj]]@meta.data[[treatment_col]])
)
```

### make.merge.pdfs
```r
make.merge.pdfs(dir = "GSE135779_SC_Expression",
    genes = genes,
    study_id = "GSE135779",
    study_title = "test")
```
