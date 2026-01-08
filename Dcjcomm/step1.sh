library(Seurat)

seurat_obj = readRDS("../../6M_hep_adipo_CCC_trans.rds")

## cell_group
# Extract cell metadata
cell_metadata <- seurat_obj@meta.data
# Assume that cell type information is stored in the 'celltype' column
cell_group <- split(rownames(cell_metadata), cell_metadata$celltype)
# Inspect the grouped cell information
str(cell_group)

## true_labs
# Extract the true cell class labels
true_labs <- seurat_obj@meta.data$seurat_clusters
# Convert classes to factors and then to numeric values
# When converting factors to numeric, each class is encoded as a corresponding integer
true_labs <- as.numeric(as.factor(true_labs))
# Convert the numeric vector to matrix format
true_labs <- as.matrix(true_labs)
# Check the data type
class(true_labs)

## W
# Extract the neighbor graph and convert it to a standard matrix
graph_matrix <- as.matrix(seurat_obj@graphs$RNA_snn)
# If necessary, ensure the matrix is symmetric (in most cases it already is)
W <- (graph_matrix + t(graph_matrix)) / 2
# Inspect the matrix structure
str(W)

# Step 1: Extract the celltype column
celltype_info <- seurat_obj@meta.data$celltype

# Step 2: Convert the extracted column to matrix and array formats
Label <- matrix(celltype_info, ncol = 1)
