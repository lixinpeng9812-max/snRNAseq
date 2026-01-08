# Load required libraries
library(tidyverse)
library(Seurat)
library(harmony)

# Read Seurat object
pbmc <- readRDS("AF_DE.rds")

# Check the number of cells per sample
table(pbmc@meta.data$orig.ident)

# Inspect metadata columns
colnames(pbmc@meta.data)

# Visualize QC metrics (UMI counts, gene counts, mitochondrial percentage)
VlnPlot(
  pbmc,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mito"),
  group.by = "orig.ident",
  pt.size = 0.01,
  ncol = 3
)

# Quality control filtering
scRNA <- subset(
  pbmc,
  subset = nCount_RNA < 20000 &
           nCount_RNA > 1000 &
           nFeature_RNA > 500 &
           percent.mito < 5
)

# View Seurat objects before and after filtering
pbmc
scRNA

# Data normalization (log-normalization)
scRNA <- NormalizeData(
  scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identification of highly variable genes
scRNA <- FindVariableFeatures(
  scRNA,
  selection.method = "vst",
  nfeatures = 2000
)

# Data scaling based on highly variable genes
scale.genes <- VariableFeatures(scRNA)
scRNA <- ScaleData(
  scRNA,
  features = scale.genes
)

# Principal component analysis (PCA)
scRNA <- RunPCA(
  scRNA,
  features = VariableFeatures(scRNA),
  verbose = FALSE
)

# Batch effect correction using Harmony
scRNA_harmony <- RunHarmony(
  scRNA,
  group.by.vars = "orig.ident"
)

# Graph-based clustering
scRNA_harmony1 <- FindNeighbors(
  scRNA_harmony,
  reduction = "harmony",
  dims = 1:20
) %>%
  FindClusters(resolution = 0.8)

# Nonlinear dimensionality reduction
scRNA_harmony1 <- RunUMAP(
  scRNA_harmony1,
  reduction = "harmony",
  dims = 1:20
)

scRNA_harmony1 <- RunTSNE(
  scRNA_harmony1,
  reduction = "harmony",
  dims = 1:20
)

# Visualization of clustering results
DimPlot(scRNA_harmony1, reduction = "umap", label = TRUE)
DimPlot(scRNA_harmony1, reduction = "tsne", label = TRUE)

# Visualization by annotated cell types
DimPlot(
  scRNA_harmony1,
  reduction = "tsne",
  label = TRUE,
  group.by = "celltype"
)

# Visualization by sample origin
DimPlot(
  scRNA_harmony1,
  reduction = "umap",
  group.by = "orig.ident"
)

DimPlot(
  scRNA_harmony1,
  reduction = "umap",
  split.by = "orig.ident"
)
