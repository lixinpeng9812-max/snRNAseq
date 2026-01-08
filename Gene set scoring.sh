library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(clusterProfiler)
library(AUCell)

sc2 = readRDS("./AF_adipo_sub_anno.rds")
DimPlot(sc2, reduction = "umap", label = TRUE)

exp = GetAssayData(sc2, slot = "data")
cells_rankings <- AUCell_buildRankings(
  exp,
  nCores = 12,
  plotStats = TRUE
)

c2 <- read.csv("glycolytic process.csv")
geneSets <- lapply(unique(c2$term), function(x) {
  print(x)
  c2$gene[c2$term == x]
})
names(geneSets) <- unique(c2$term)

cells_AUC <- AUCell_calcAUC(
  geneSets,
  cells_rankings,
  nCores = 12,
  aucMaxRank = nrow(cells_rankings) * 0.1
)

geneSet <- "glycolytic process"

aucs <- as.numeric(getAUC(cells_AUC)[geneSet, ])
sc2$AUC <- aucs

df <- data.frame(
  sc2@meta.data,
  sc2@reductions$umap@cell.embeddings
)

class_avg <- df %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

colnames(df)

ggplot(df, aes(umap_1, umap_2)) +
  geom_point(aes(colour = AUC)) +
  viridis::scale_color_viridis(option = "H") +  # Change the viridis color palette option
  ggrepel::geom_label_repel(
    aes(label = celltype),
    data = class_avg,
    size = 3,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") +
  theme_bw() +
  facet_grid(. ~ celltype)

ggplot(df, aes(umap_1, umap_2)) +
  geom_point(aes(colour = AUC)) +
  viridis::scale_color_viridis(option = "H") +  # Change the viridis color palette option
  ggrepel::geom_label_repel(
    aes(label = celltype),
    data = class_avg,
    size = 3,
    label.size = 1,
    segment.color = NA
  ) +
  theme(legend.position = "none") +
  theme_bw()
