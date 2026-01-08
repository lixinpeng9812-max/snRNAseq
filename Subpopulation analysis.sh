library(Seurat)
library(tidyverse)
library(reshape2)
library(tidyverse)
library(harmony)

scRNA_harmony1 = readRDS("../AF_anno.rds")
Idents(scRNA_harmony1) <- scRNA_harmony1@meta.data$celltype
Idents(scRNA_harmony1)

sub = subset(scRNA_harmony1, idents = c("Adipo"))
sub

scRNA <- NormalizeData(sub, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable genes
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# Based on highly variable genes
scale.genes <- VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA), verbose = FALSE)
scRNA <- FindNeighbors(scRNA, dims = 1:20) %>% FindClusters(resolution = 1)
scRNA <- RunUMAP(scRNA, dims = 1:20)
scRNA <- RunTSNE(scRNA, dims = 1:20)

DimPlot(scRNA, reduction = "umap", label = TRUE) |
  DimPlot(scRNA, reduction = "umap", group.by = "group")

# Harmony integration
library(harmony)
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA), verbose = FALSE)
scRNA_harmony <- RunHarmony(scRNA, group.by.vars = 'orig.ident')

scRNA_harmony1 <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1)

scRNA_harmony1 <- RunUMAP(scRNA_harmony1, reduction = "harmony", dims = 1:20)
scRNA_harmony1 <- RunTSNE(scRNA_harmony1, reduction = "harmony", dims = 1:20)

Idents(scRNA_harmony1) <- scRNA_harmony1@meta.data$seurat_clusters
DimPlot(scRNA_harmony1, reduction = "umap", label = TRUE)
DimPlot(scRNA_harmony1, reduction = "umap", label = TRUE, split.by = "group")
