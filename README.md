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

### remove_low_count_cells
```r
remove_low_count_cells(seurat_obj = obj, metadata_column = "monaco.main", threshold = 20)
```

### read_seurat_objects_from_files
```r
read_seurat_objects_from_files(file_paths = paths, default_assay = "RNA")
```


### frequency_table
```r
# Define inputs
seurat_object <- data2  # Assuming data2 is your Seurat object
meta.col <- "Type" # Treatment/stim column
group <- c("Inflammed", "Noninflammed", "Control") # stims
gene_name <- c("PD1", "PDL1") # Gene common names
filter_column <- c("PDCD1", "CD274") # Gene names present in data
ident <- "monaco.fine"

# Generate frequency table
result <- frequency_table(seurat_object = seurat_object,
                          meta.col = meta.col,
                          group = group,
                          gene_name = gene_name,
                          filter_column = filter_column,
                          ident = ident)
```
