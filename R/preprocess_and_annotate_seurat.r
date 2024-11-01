suppressPackageStartupMessages({
    library(Seurat)
    library(celldex)
})

# Function to preprocess and annotate Seurat object
preprocess_and_annotate_seurat <- function(obj, ncells = 3000, variable.features.n = 2000, vst.flavor = "v2", 
                                            npcs = 30, umap_dims = 1:20, resolution = 0.5) {
  
  # Preprocessing steps
  obj <- NormalizeData(obj)
  obj <- SCTransform(obj, ncells = ncells, variable.features.n = variable.features.n, 
                           vst.flavor = vst.flavor, method = 'glmGamPoi', conserve.memory = FALSE)
  obj <- RunPCA(obj, npcs = npcs, verbose = TRUE)
  obj <- RunUMAP(obj, reduction = "pca", dims = umap_dims)
  obj <- FindNeighbors(obj, reduction = "pca", dims = umap_dims)
  obj <- FindClusters(obj, resolution = resolution)
  
  # Load reference annotation data (replace with your actual loading logic)
  monaco.ref <- celldex::MonacoImmuneData()
  blueprint.ref <- celldex::BlueprintEncodeData()
  
  # Blueprint annotation (assuming LayerData and SingleR exist)
  sce <- LayerData(obj)
  blueprint.main <- SingleR(test = sce, assay.type.test = 1, ref = blueprint.ref, labels = blueprint.ref$label.main)
  obj@meta.data$blueprint.main <- blueprint.main$pruned.labels
  
  # Monaco annotation (assuming LayerData and SingleR exist)
  sce <- LayerData(obj)
  monaco.fine <- SingleR(test = sce, assay.type.test = 1, ref = monaco.ref, labels = monaco.ref$label.fine)
  obj@meta.data$monaco.fine <- monaco.fine$pruned.labels
  
  # Return the preprocessed and annotated Seurat object
  return(obj)
}
