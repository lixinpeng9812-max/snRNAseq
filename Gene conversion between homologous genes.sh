library(biomaRt)
library(Seurat)

rds <- readRDS("../scRNA_sub/liver_celltype_anno.rds")

pig_features <- rownames(rds@assays$RNA@counts)

ensembl <- useMart("ensembl")
DBlist = listDatasets(ensembl)
dim(DBlist)
head(DBlist, 5)

ensembl.pig <- useMart("ensembl", dataset = "sscrofa_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

listAttributes(ensembl.pig)

pig_features[2400:2450]

pig_features_id_l <- grepl("ENSSSCG", pig_features)
pig_features_symbol <- pig_features[!pig_features_id_l]
pig_features_id <- pig_features[pig_features_id_l]

pig2human.gene_symbol <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name"),
                                filters = c("external_gene_name"),
                                values = pig_features_symbol,
                                mart = ensembl.pig,
                                attributesL = c("ensembl_gene_id", "external_gene_name"),
                                martL = ensembl.human,
                                uniqueRows = T)

pig2human.gene_id <- getLDS(attributes = c("ensembl_gene_id", "external_gene_name"),
                            filters = c("ensembl_gene_id"),
                            values = pig_features_id,
                            mart = ensembl.pig,
                            attributesL = c("ensembl_gene_id", "external_gene_name"),
                            martL = ensembl.human,
                            uniqueRows = T)

#去除未成功转换和重复的gene
pig2human.gene_id <- pig2human.gene_id[pig2human.gene_id$Gene.name.1!="",]
dim(pig2human.gene_id)
pig2human.gene_symbol <- pig2human.gene_symbol[pig2human.gene_symbol$Gene.name.1!="",]
dim(pig2human.gene_symbol)

pig2human.gene_id <- pig2human.gene_id[duplicated(pig2human.gene_id$Gene.stable.ID)==F,]
dim(pig2human.gene_id)
pig2human.gene_symbol <- pig2human.gene_symbol[duplicated(pig2human.gene_symbol$Gene.name)==F,]
dim(pig2human.gene_symbol)

pig2human.gene_id$Gene.name <- pig2human.gene_id$Gene.stable.ID
pig2human <- rbind(pig2human.gene_id, pig2human.gene_symbol)
dim(pig2human)

pig2human <- pig2human[duplicated(pig2human$Gene.name.1)==F,]
dim(pig2human)
saveRDS(pig2human, "pig2human.rds")

##################################################################################################################################
pig2human = readRDS("./pig2human.rds")
scRNA = scRNA_harmony1
pig_data_trans = subset(scRNA, features = pig2human$Gene.name)

RenameGenesSeurat_v2 <- function(obj,newnames,gene.use=NULL,de.assay="RNA") {
  print("Run this before integration. It only changes obj@assays$*@counts, @data and @scale.data, @var.features,@reductions$pca@feature.loadings")
  lassays <- Assays(obj)
  assay.use <- obj@reductions$pca@assay.used
  DefaultAssay(obj) <- de.assay
  if (is.null(gene.use)) {
    all_genenames <- rownames(obj)
  }else{
    all_genenames <- gene.use
    obj <- subset(obj,features=gene.use)
  }
  
  order_name <- function(v1,v2,ref){
    v2 <- make.names(v2,unique=T)
    df1 <- data.frame(v1,v2)
    rownames(df1) <- df1$v1
    df1 <- df1[ref,]
    return(df1)
  }
  
  df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
  all_genenames <- df1$v1
  newnames <- df1$v2
  
  if ('SCT' %in% lassays) {
    if ('SCTModel.list' %in%  slotNames(obj@assays$SCT)) {
      obj@assays$SCT@SCTModel.list$model1@feature.attributes <- obj@assays$SCT@SCTModel.list$model1@feature.attributes[all_genenames,]
      rownames(obj@assays$SCT@SCTModel.list$model1@feature.attributes) <- newnames
    }
  }
  change_assay <- function(a1=de.assay,obj,newnames=NULL,all_genenames=NULL){
    RNA <- obj@assays[a1][[1]]
    if (nrow(RNA) == length(newnames)) {
      if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
      if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
      if (length(RNA@var.features)) {
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=RNA@var.features)
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        RNA@var.features            <- newnames1
      }
      if (length(RNA@scale.data)){
        df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(RNA@scale.data))
        all_genenames1 <- df1$v1
        newnames1 <- df1$v2
        rownames(RNA@scale.data)    <- newnames1
      }
      
    } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
    obj@assays[a1][[1]] <- RNA
    return(obj)
  }
  
  for (a in lassays) {
    DefaultAssay(obj) <- a
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=rownames(obj))
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    obj <- change_assay(obj=obj,a1=a,newnames=newnames1,all_genenames=all_genenames1)
  }
  
  hvg <- VariableFeatures(obj,assay=assay.use)
  if (length(obj@reductions$pca)){
    df1 <- order_name(v1=all_genenames,v2=newnames,ref=hvg)
    all_genenames1 <- df1$v1
    newnames1 <- df1$v2
    rownames(obj@reductions$pca@feature.loadings) <- newnames1
  }
  
  return(obj)
}

pig_data_trans <- RenameGenesSeurat_v2(pig_data_trans, 
                                       newnames = pig2human$Gene.name.1,
                                       gene.use = pig2human$Gene.name,
                                       de.assay = 'RNA')
FeaturePlot(pig_data_trans, features = c("TBL1Y")) | FeaturePlot(scRNA, features = c("ENSSSCG00000012102"))

