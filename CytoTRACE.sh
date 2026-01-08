#### Load packages
library(tidydr)
library(ggplot2)
library(CytoTRACE)
library(Seurat)

#### Load data
seurat_obj = scRNA_harmony1

seurat_obj$celltype <- Idents(seurat_obj)

#### Extract phenotype information
phe <- seurat_obj$celltype
phe = as.character(phe)
names(phe) <- rownames(seurat_obj@meta.data)

#### Extract expression matrix
mat_3k <- as.matrix(seurat_obj[["RNA"]]$counts)
mat_3k[1:4, 1:4]

#### When the number of cells in the dataset exceeds 3,000,
#### CytoTRACE will automatically run in fast mode, which is a
#### subsampling-based method to reduce runtime and memory usage.
#### In addition, users can enable multithreading via ncores (default = 1),
#### or specify the subsampling size via subsamplingsize (default = 1000 cells).
#### Here, the dataset is run in fast mode using 8 cores with a subsampling size of 1000.
results <- CytoTRACE(mat = mat_3k)

plotCytoGenes(results, numOfGenes = 10)  # The number of genes can be adjusted

emb <- seurat_obj@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phe, emb = emb) 
# If emb is not provided, t-SNE will be automatically computed for visualization

#### Embed CytoTRACE results into the Seurat object
seurat_obj[["CytoTRACE"]] <- results$CytoTRACE

FeaturePlot(seurat_obj, features = "CytoTRACE", reduction = "umap") +
  ggtitle("CytoTRACE on UMAP")
