# Load required libraries
library(Seurat)
library(dplyr)
library(pheatmap)
library(cols4all)
library(ComplexHeatmap)
library(circlize)

# Read annotated Seurat object
seurat_object <- readRDS("./liver_anno.rds")

# Define marker genes for major liver cell types
genes <- c("HNF4A", "ALB", "PECAM1", "KDR", "RELN", "DCN",
           "AFF3", "PTPRC", "CD38", "VSIG4", "CD163", "CFTR")

# Calculate average expression of selected genes by cell type
aver_dt <- AverageExpression(
  seurat_object,
  features = genes,
  group.by = "celltype",
  slot = "data"
)

# Convert to data frame
aver_dt <- as.data.frame(aver_dt$RNA)
aver_dt[1:6, 1:6]

# Basic heatmap visualization
pheatmap(
  as.matrix(aver_dt),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE
)

# Row annotation: marker gene–associated cell types (Top marker genes)
gene_anno <- data.frame(
  gene_anno = markergene$cluster,
  row.names = markergene$gene
)

# Column annotation: cell types
cell_anno <- data.frame(
  cell_anno = colnames(aver_dt),
  row.names = colnames(aver_dt)
)

head(gene_anno)
head(cell_anno)

# Heatmap with row and column annotations
pheatmap(
  as.matrix(aver_dt),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = cell_anno,
  annotation_row = gene_anno
)

# Customize heatmap color palette
mycol <- colorRampPalette(c("#5E3C99", "white", "#E66101"))(50)

# Define annotation colors
celltype_col <- c4a("10", 6)
names(celltype_col) <- cell_anno$cell_anno

anno_col <- list(
  cell_anno = celltype_col,
  gene_anno = celltype_col
)

# Enhanced heatmap with customized colors
pheatmap(
  as.matrix(aver_dt),
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = cell_anno,
  annotation_row = gene_anno,
  annotation_colors = anno_col,
  color = mycol,
  border_color = "white"
)

# Define continuous color mapping for ComplexHeatmap
mycol2 <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))

# Z-score normalization by gene
aver_dtt <- t(scale(t(aver_dt)))

# Prepare colors for cell type annotations
cols <- c4a("classic10light", 6)
names(cols) <- unique(cell_anno$cell_anno)

# Column annotation (cell types)
cell <- data.frame(cell = cell_anno$cell_anno)
rownames(cell) <- colnames(aver_dtt)

col_anno <- HeatmapAnnotation(
  df = cell,
  show_annotation_name = FALSE,
  gp = gpar(col = "white", lwd = 2),
  col = list(cell = cols)
)

# Define gene-to-cell-type mapping for row annotation
gene_cell_mapping <- c(
  "HNF4A" = "Hep",
  "ALB" = "Hep",
  "PECAM1" = "LSEC",
  "KDR" = "LSEC",
  "RELN" = "HSC",
  "DCN" = "HSC",
  "AFF3" = "Immune",
  "PTPRC" = "Immune",
  "CD38" = "Immune",
  "VSIG4" = "Kuffer",
  "CD163" = "Kuffer",
  "CFTR" = "Chol"
)

row_cols <- cols[gene_cell_mapping[rownames(aver_dtt)]]

# Row annotation with colored gene labels
row_anno <- rowAnnotation(
  foo = anno_text(
    rownames(aver_dtt),
    location = 0,
    just = "left",
    gp = gpar(
      fill = row_cols,
      col = "black",
      fontface = "italic"
    ),
    width = max_text_width(rownames(aver_dtt)) * 1.25
  )
)

# Define column (cell type) order
column_order <- c("Hep", "LSEC", "HSC", "Immune", "Kuffer", "Chol")

# Final ComplexHeatmap visualization
Heatmap(
  aver_dtt,
  name = "expression",
  column_order = column_order,
  col = mycol2,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  column_names_side = "top",
  column_names_rot = 60,
  row_names_gp = gpar(fontsize = 12, fontface = "italic"),
  rect_gp = gpar(col = "white", lwd = 1.5),
  top_annotation = col_anno
) + row_anno
