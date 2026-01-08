## devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

# Load Seurat object
pbmc = readRDS("./TP_AF.rds")

# Check sample composition
table(pbmc@meta.data$orig.ident)

# Set cell identities based on sample origin
Idents(pbmc) = pbmc@meta.data$orig.ident

# Subset cells from abdominal fat of 6-month-old LTP pigs
scRNA_sub = subset(pbmc, idents = "LTP_6M_3_AbdominalFat")

# Quality control filtering
scRNA <- subset(
  scRNA_sub,
  subset =
    nCount_RNA < 20000 &
    nCount_RNA > 1000 &
    nFeature_RNA > 500 &
    percent.mito < 5
)

scRNA_sub
scRNA

# Data normalization
scRNA <- NormalizeData(
  scRNA,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify highly variable genes
scRNA <- FindVariableFeatures(
  scRNA,
  selection.method = "vst",
  nfeatures = 2000
)

# Data scaling based on highly variable genes
scale.genes <- VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

# Principal component analysis
scRNA <- RunPCA(
  scRNA,
  features = VariableFeatures(scRNA),
  verbose = FALSE
)

# Construct neighbor graph and perform clustering
scRNA1 <- FindNeighbors(scRNA, dims = 1:20) %>%
  FindClusters(resolution = 0.8)

# UMAP dimensionality reduction
scRNA1 <- RunUMAP(scRNA1, dims = 1:20)

# Check metadata columns
colnames(scRNA1@meta.data)

# Visualize clusters
DimPlot(scRNA1, reduction = "umap", group.by = "seurat_clusters")

# Parameter sweep to identify optimal pK for DoubletFinder
sweep.res.list <- paramSweep(scRNA1, PCs = 1:10, sct = FALSE)

# Since LogNormalize is used, set sct = FALSE
# If SCT normalization is used, set sct = TRUE
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

# Identify optimal pK value
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>%
  as.character() %>%
  as.numeric()

# Estimate expected doublet rate
scRNA
DoubletRate1 <- ncol(scRNA) * 8 * 1e-6  # General empirical estimate

# Estimate homotypic doublet proportion based on clustering
homotypic.prop <- modelHomotypic(scRNA1$seurat_clusters)

# Calculate expected number of doublets
nExp_poi <- round(DoubletRate1 * ncol(scRNA1))

# Adjust expected doublets by removing homotypic proportion
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

# Run DoubletFinder
scRNA2 <- doubletFinder(
  scRNA1,
  PCs = 1:10,
  pN = 0.25,
  pK = pK_bcmvn,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

# View metadata to identify DoubletFinder classification column
colnames(scRNA2@meta.data)

# Visualize singlets and doublets
DimPlot(
  scRNA2,
  reduction = "umap",
  group.by = "DF.classifications_0.25_0.005_424"
)

# Subset singlet cells only
cell_sub <- subset(
  scRNA2@meta.data,
  DF.classifications_0.25_0.005_424 == "Singlet"
)

scRNA_sub <- subset(scRNA2, cells = row.names(cell_sub))

# Save singlet-only Seurat object
saveRDS(scRNA_sub, "LTP_6M_3_AF_DF.rds")

scRNA_sub
