# Load required libraries
library(stringr)
library(Seurat)
library(patchwork)
library(SummarizedExperiment)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(qs)

# Set up parallel processing
register(MulticoreParam(workers = 6, progressbar = TRUE))

# Load data
loom <- open_loom("./sample_SCENIC.loom")
sce <- readRDS("./6M_diffe_preadipo.rds")

# Filter cell types
table(sce$celltype)
sce <- sce[, sce$celltype %in% c('Differentiation inhibiting Preadipo', 'Differentiation promoting Preadipo')]

# Extract regulons
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons_incidMat[1:4, 1:4] # View first 4 rows and columns

# Convert regulons to gene lists
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

# Extract regulon AUC information
regulonAUC <- get_regulons_AUC(loom, column.attr.name = 'RegulonsAUC')

# Extract regulon thresholds
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

# Extract embeddings (though none were generated in pySCENIC analysis)
embeddings <- get_embeddings(loom)
embeddings

# Subset regulonAUC to match sce object
sub_regulonAUC <- regulonAUC[, match(colnames(sce), colnames(regulonAUC))]
dim(sub_regulonAUC)
sce # Check dimensions

# Verify consistency
identical(colnames(sub_regulonAUC), colnames(sce))

# Create cell type annotations
cellTypes <- data.frame(row.names = colnames(sce),
                        celltype = sce$celltype)
head(cellTypes)

sub_regulonAUC[1:4, 1:4]
table(sce$celltype)

# Split cells by group
selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), cellTypes[, selectedResolution])

# Keep only non-duplicated extended regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)), ]
dim(sub_regulonAUC)

# Calculate Regulon Specificity Scores (RSS)
selectedResolution <- "celltype"
rss <- calcRSS(AUC = getAUC(sub_regulonAUC),
               cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution])

# Remove NA values
rss <- na.omit(rss)

# Plot RSS heatmap
pdf('RSS_plot.pdf', width = 10, height = 20)
rssPlot <- plotRSS(
  rss,
  labelsToDiscard = NULL,  # Specify labels to exclude from heatmap
  zThreshold = 0.3,        # Regulon threshold (default = 1)
  cluster_columns = FALSE, # Whether to cluster columns
  order_rows = TRUE,       # Whether to sort rows
  thr = 0.01,             # Threshold for filtering RSS values (default = 0.01)
  varName = "cellType",
  col.low = '#330066',    
  col.mid = '#66CC66',    
  col.high = '#FFCC33',
  revCol = FALSE,
  verbose = TRUE
)
rssPlot$plot
dev.off()

# Plot RSS for specific cell type (Monocyte)
i <- 'Differentiation inhibiting Preadipo'
plotRSS_oneSet(rss, setName = i)


# Calculate mean regulon activity per cell group
# cellsPerGroup contains lists of samples for each cell population
# rowMeans(getAUC(sub_regulonAUC)[,x] calculates mean AUC for each regulon per cell group
regulonActivity_byGroup <- sapply(
  cellsPerGroup,
  function(x) rowMeans(getAUC(sub_regulonAUC)[, x])
)

# Check the range of regulon activities
range(regulonActivity_byGroup)

# Normalize the results (z-score scaling)
# Scaling is performed per regulon across different clusters
regulonActivity_byGroup_Scaled <- t(scale(
  t(regulonActivity_byGroup),
  center = TRUE,
  scale = TRUE
))

# Check dimensions and remove NA values
dim(regulonActivity_byGroup_Scaled)
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

# Plot heatmap using ComplexHeatmap
library(ComplexHeatmap)
library(circlize)

Heatmap(
  regulonActivity_byGroup_Scaled,
  name = "z-score",
  col = colorRamp2(
    seq(from = -2, to = 2, length = 11),
    rev(brewer.pal(11, "Spectral"))
  ),
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 12),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot = 0,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE
)

# Create Heatmap object from RSS data
p <- Heatmap(rss)

# Convert Heatmap matrix to data frame and calculate difference between groups
rss_df <- data.frame(p@matrix)
rss_df$diff <- rss_df[, 1] - rss_df[, 2]
head(rss_df)

# Select top 3 regulons with highest and lowest differences
top3_max <- rss_df %>% 
  top_n(3, diff) %>% 
  rownames()

top3_min <- rss_df %>% 
  top_n(-3, diff) %>% 
  rownames()

# Visualize RSS heatmaps
library(pheatmap)

# Full RSS heatmap
pheatmap(rss)

# Heatmap of top and bottom 3 differential regulons
pheatmap(rss[c(top3_max, top3_min), ])

# Add regulon AUC scores to Seurat object metadata
sce@meta.data <- cbind(sce@meta.data, t(sub_regulonAUC@assays@data@listData[["AUC"]]))
Idents(sce) <- sce$celltype

# Create dot plot of top regulons
DotPlot(sce, 
        features = c("FOXO1(+)", "PRRX2(+)", "ZNF200(+)"),
        cols = 'RdBu')

# Create clustered heatmap of top regulons using SCP package
library(SCP)
pdf('Top3_regulon_clustered_dotplot2.pdf', height = 8, width = 6)
GroupHeatmap(
  sce,
  features = c(top3_max, top3_min),
  group.by = "celltype",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = TRUE,
  cluster_column_slices = TRUE,
  add_dot = TRUE,
  add_reticle = TRUE,
  heatmap_palette = "viridis",
  nlabel = 0,
  show_row_names = TRUE,
  ht_params = list(
    row_gap = unit(0, "mm"),
    row_names_gp = gpar(fontsize = 10)
  )
)
dev.off()
