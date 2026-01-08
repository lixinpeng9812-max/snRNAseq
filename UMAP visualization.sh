# Load required R packages
library(tidydr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cols4all)
library(RColorBrewer)

# Read scRNA-seq data
## Preprocessed and annotated Seurat object
seurat_object <- readRDS("./AF_anno.rds")
dim(seurat_object)

# Extract cell type annotations
celltype <- Idents(seurat_object)
table(celltype)

# Extract UMAP embeddings
UMAP <- as.data.frame(seurat_object@reductions$umap@cell.embeddings)
UMAP <- cbind(UMAP, celltype)
head(UMAP)

# Define a customized theme
mytheme <- theme_void() + 
  theme(plot.margin = margin(5.5, 15, 5.5, 5.5))

# Calculate median UMAP coordinates for each cell type
## These coordinates are used for cluster labeling
label <- UMAP %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
head(label)

# Define color palette
display.brewer.all(type = "div")
mycol <- brewer.pal(8, "Set2")

# UMAP visualization with ellipses and cluster labels
ggplot(data = UMAP, aes(x = umap_1, y = umap_2)) +
  
  # Scatter plot of single cells
  geom_point(
    aes(color = celltype),
    size = 0.4,
    alpha = 0.8
  ) +
  
  # Dashed confidence ellipses (excluding "Unknow" cluster)
  stat_ellipse(
    data = ~ subset(.x, celltype != "Unknow"),
    aes(color = celltype),
    level = 0.95,
    linetype = 2,
    show.legend = FALSE
  ) +
  
  # Filled confidence ellipses (excluding "Unknow" cluster)
  stat_ellipse(
    data = ~ subset(.x, celltype != "Unknow"),
    aes(color = celltype, fill = celltype),
    level = 0.95,
    linetype = 1,
    geom = "polygon",
    alpha = 0.1,
    show.legend = FALSE
  ) +
  
  # Add customized axis arrows
  theme_dr(
    xlength = 0.2,
    ylength = 0.2,
    arrow = grid::arrow(
      length = unit(0.1, "inches"),
      ends = "last",
      type = "closed"
    )
  ) +
  
  # Add cluster labels
  geom_text(
    data = label,
    aes(x = umap_1, y = umap_2, label = celltype),
    fontface = "bold",
    color = "black",
    size = 4
  ) +
  
  # Legend and color settings
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol) +
  
  # Theme adjustments
  theme(
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
