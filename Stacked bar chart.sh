library(Seurat)

# Load Seurat object
scRNA <- readRDS("./AF_anno.rds")

# Set cell identities as cell type
Idents(scRNA) <- "celltype"

# Count the number of cells of each cell type in each group
table(Idents(scRNA), scRNA$group)

# Calculate the proportion of each cell type within each group
Cellratio <- prop.table(
  table(Idents(scRNA), scRNA$group),
  margin = 2
)

# Convert to data frame
Cellratio <- as.data.frame(Cellratio)

# Rename columns
colnames(Cellratio) <- c("celltype", "group", "ratio")

# Adjust the order of groups on the x-axis
Cellratio$group <- factor(Cellratio$group, levels = c())

library(ggplot2)
library(ggalluvial)
library(RColorBrewer)

# Choose color palette
mycol <- brewer.pal(8, "Set2")

# Plot stacked bar chart with alluvial flows
ggplot(
  Cellratio,
  aes(
    x = group,
    y = ratio,
    fill = celltype,
    stratum = celltype,
    alluvium = celltype
  )
) +
  geom_col(
    width = 0.7,
    color = "white",
    size = 0.9
  ) +
  geom_flow(
    width = 0.7,
    alpha = 0.7,
    knot.pos = 0,
    color = "white",
    size = 1
  ) +
  scale_fill_manual(values = mycol) +
  theme_classic() +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      size = 0.5,
      linetype = "solid"
    ),
    aspect.ratio = 0.8,
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  labs(
    x = "Group",
    y = "Cell proportion",
    fill = "Cell type"
  )
