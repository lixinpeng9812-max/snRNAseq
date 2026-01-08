library(Seurat)
library(tidyverse)
library(reshape2)
library(tidyverse)
library(harmony)
library(monocle)


sc = readRDS("../preadipo_adipo_sub_anno.rds")
table(sc@meta.data$celltype)
sc1 = subset(sc, idents = c("Differentiation inhibiting Preadipo(AF)", "Differentiation promoting Preadipo(AF)"))
Idents(sc1) = "group"
sc2 = subset(sc1, idents = c("6M"))
sc2

scRNA.Osteoclastic= sc2
data <- as(as.matrix(scRNA.Osteoclastic@assays$RNA$counts), 'sparseMatrix')
scRNA.Osteoclastic
pd <- new('AnnotatedDataFrame', data = scRNA.Osteoclastic@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())



monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

HSMM=monocle_cds
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
HSMM <- setOrderingFilter(HSMM, disp.genes)
plot_ordering_genes(HSMM)

diff = differentialGeneTest(HSMM[disp.genes,], fullModelFormulaStr = "~celltype", cores = 12)
#diff = differentialGeneTest(HSMM, fullModelFormulaStr = "~celltype", cores = 8)
deg = subset(diff, qval < 0.05)
deg = deg[order(deg$qval, decreasing = F),]
ordergene = rownames(deg)

HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')

HSMM <- orderCells(HSMM)
HSMM <- orderCells(HSMM, root_state=3)
HSMM <- orderCells(HSMM, reverse = T)

plot_cell_trajectory(HSMM, color_by = "Pseudotime") |  plot_cell_trajectory(HSMM, color_by = "celltype")
plot_cell_trajectory(HSMM, color_by = "group")

Time_diff = differentialGeneTest(HSMM[ordergene,], cores = 12,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff = Time_diff[,c(5,2,3,4,1,6,7)]
Time_deg = subset(Time_diff, qval < 0.05 )
write.csv(Time_deg, "6M_preadipo_Time_diff.csv", quote = F, row.names = F)
Time_genes = Time_deg %>% pull(gene_short_name) %>% as.character()


p = plot_pseudotime_heatmap(HSMM[Time_genes,], 
                            num_clusters = 2, 
                            show_rownames = F, 
                            return_heatmap = T, 
                            cores = 12)
p
cluster = cutree(p$tree_row, k=3)
clustering = data.frame(cluster)
clustering[,1] = as.character(clustering[,1])
colnames(clustering) = "Gene_Cluster"
table(clustering)
write.csv(clustering, "6M_preadipo_Time_cluster.csv", quote = F)
