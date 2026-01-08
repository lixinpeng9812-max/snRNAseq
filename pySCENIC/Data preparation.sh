library(Seurat)
sce = readRDS("pbmc.rds")
table(sce@meta.data$celltype)
table(sce@meta.data$group)
write.csv(t(as.matrix(sce@assays$RNA@counts)),file = "sce_exp.csv")
cellInfo <- sce@meta.data[,c("celltype","nCount_RNA","nFeature_RNA")]
colnames(cellInfo) <- c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
write.csv(cellInfo, file = "cellInfo.csv")
