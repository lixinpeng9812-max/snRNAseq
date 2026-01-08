####################################################
#### PART 0: Data collection and preprocessing #####
####################################################

# Load DcjComm package and example data
library(DcjComm) ##load the DcjComm package
library('Matrix')
library('scRMD')
library(doParallel)
library(dplyr)
library(Seurat)

data("lr_human") ##ligand-receptor pairs databases (L-R)
LR_pairs_DB = lr_human[,c("ligand","receptor")]
data("TF_PPRhuman") ##Receptor-Transcription factor a-priori association(TF_PPR)
data("TF_TGhumandb") ##Transcriptional regulatory networks (TF_TG)
#data("lst_scrna")
#data <- lst_scrna$mat
data <- GetAssayData(seurat_obj, slot = "data")
data <- as.matrix(data)

# Select marker genes by Seurat
library(tidyverse)
library(Matrix)
library(Seurat)

#data<-log2(data+1)
#pbmc <- CreateSeuratObject(counts = data, project = "scrna", min.cells = 0, min.features = 0)
group<-Label
seurat_obj@meta.data[["dcjcomm"]]<-group
Idents(seurat_obj) <- "celltype"
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.25)
markers <- markers %>% dplyr::select(gene, everything()) %>% subset(p_val < 0.05)
row<-unique(markers$gene,markers$gene %in% rownames(data))
gene_expr <- data[row,] ##Infer cell-cell communications by marker genes

LR_db <- lr_human
TF_reg_DB <- TF_TGhumandb
R_TF_association <- TF_PPRhuman
Org <- "Homo sapiens" #choose the species source of gene, eg "Homo sapiens", "Mus musculus"

# Infer cell-cell communications
CCI <- Compute_CCI(
  data = gene_expr, Label = Label, cell_group = cell_group, LR_pairs_DB = LR_pairs_DB,
  TF_reg_DB = TF_reg_DB, R_TF_association = R_TF_association, LR_db = LR_db,
  N_cores = 16, DEmethod = "wilcoxon", n_bootstrap = 4, backend = "doParallel",
  other_inter_scores = NULL, cell_reference = NULL, Org = Org,
  use.type="median",probs = 0.9,method="weighted"
)
expr_l_r <-CCI$expr_l_r_log2_scale ##The results of cell-cell communications
