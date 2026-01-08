library(Seurat)
mt.genes=c('ENSSSCG00000018065','ENSSSCG00000018069',
           'ENSSSCG00000018075','ENSSSCG00000018078','ENSSSCG00000018080',
           'ENSSSCG00000018081','ENSSSCG00000018082','ENSSSCG00000018084',
           'ENSSSCG00000018086','ENSSSCG00000018087',
           'ENSSSCG00000018091','ENSSSCG00000018092','ENSSSCG00000018094')
mt.genes = c("ATP6", "ATP8", "COX1", "COX2", "COX3", "CYTB", "ND2", "ND3",
             "ND4", "ND4L", "ND5", "ND6")

pbmc.combined = readRDS("./pbmc.rds")
pbmc.combined@meta.data$percent.mito = "NA"
pbmc.combined[["percent.mito"]]= PercentageFeatureSet(
  pbmc.combined,
  features = c("ATP6", "COX1", "COX2", "COX3", "CYTB", "ND2", "ND3",
               "ND4", "ND4L", "ND5", "ND6")
)
