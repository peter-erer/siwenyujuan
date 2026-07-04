library(Seurat)
library(monocle3)
library(SingleCellExperiment)

pbmc <- readRDS("pbmc_annotated_for_downstream.rds")

expr_mat <- GetAssayData(
  pbmc,
  assay = "RNA",
  layer = "counts"
)

cell_metadata <- pbmc@meta.data
cell_metadata$cell_id <- rownames(cell_metadata)

gene_metadata <- data.frame(
  gene_short_name = rownames(expr_mat),
  row.names = rownames(expr_mat)
)

cds <- new_cell_data_set(
  expression_data = expr_mat,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(cds, num_dim = 20)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds)

plot_cells(
  cds,
  color_cells_by = "new.cluster.ids",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)