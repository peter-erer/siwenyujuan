options(stringsAsFactors = FALSE)

# =========================
# 0. 参数设置
# =========================
# 第一次运行时建议重点检查：
# 1. group_by：用于细胞分群的元数据列名
# 2. root_celltype：拟时序分析的起始细胞类型
# 3. species：人类填 "human"，小鼠填 "mouse"
config <- list(
  data_rds = "pbmc_annotated_for_downstream.rds/pbmc_annotated_for_downstream.rds",
  output_dir = "homework6_output",
  assay = NULL,
  species = "human",
  group_by = "new.cluster.ids",
  root_celltype = "Myeloid",
  root_cells_n = 200,
  monocle_num_dim = 50,
  cellchat_min_cells = 10,
  future_workers = 4,
  seed = 1234
)

# =========================
# 1. 加载所需 R 包
# =========================
required_pkgs <- c(
  "Seurat",
  "SeuratWrappers",
  "monocle3",
  "CellChat",
  "ComplexHeatmap",
  "ggplot2",
  "patchwork",
  "future"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall them first, then rerun this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(monocle3)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(future)
})

set.seed(config$seed)

# =========================
# 2. 创建输出目录
# =========================
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$output_dir, "monocle3"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$output_dir, "cellchat"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$output_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

# =========================
# 3. 定义辅助函数
# =========================
find_first_existing <- function(candidates, available) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) {
    return(NA_character_)
  }
  hit[[1]]
}

save_base_plot <- function(expr, filename, width = 8, height = 8) {
  pdf(filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

save_heatmap_plot <- function(ht, filename, width = 10, height = 8) {
  pdf(filename, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  ComplexHeatmap::draw(ht)
}

# =========================
# 4. 读取 Seurat 对象并检查基础信息
# =========================
message("读取 Seurat 对象：", config$data_rds)
obj <- readRDS(config$data_rds)

if (!inherits(obj, "Seurat")) {
  stop("输入文件不是 Seurat 对象，请先检查数据文件。")
}

if (!is.null(config$assay)) {
  DefaultAssay(obj) <- config$assay
}
assay_use <- DefaultAssay(obj)

metadata <- obj[[]]
available_meta <- colnames(metadata)

group_candidates <- c(
  config$group_by,
  "celltype",
  "cell_type",
  "celltype.l2",
  "celltype_l2",
  "celltype.l1",
  "celltype_l1",
  "annotation",
  "annot",
  "predicted.celltype.l2",
  "predicted.celltype.l1",
  "seurat_clusters"
)

group_by <- find_first_existing(group_candidates, available_meta)
if (is.na(group_by)) {
  obj$analysis_group <- as.character(Idents(obj))
  group_by <- "analysis_group"
  message("未识别到常见细胞注释列，改用 Idents(obj) 作为分组信息。")
} else {
  obj$analysis_group <- as.character(metadata[[group_by]])
}

obj$analysis_group <- ifelse(
  is.na(obj$analysis_group) | obj$analysis_group == "",
  "Unknown",
  obj$analysis_group
)
Idents(obj) <- factor(obj$analysis_group)

message("默认 assay：", assay_use)
message("分组列：", group_by)
message("元数据列：", paste(available_meta, collapse = ", "))

# 输出各细胞类型数量统计
group_table <- sort(table(obj$analysis_group), decreasing = TRUE)
write.csv(
  as.data.frame(group_table),
  file = file.path(config$output_dir, "tables", "cell_group_counts.csv"),
  row.names = FALSE
)

# 如果对象中已有 UMAP，则先保存一个基础分群图
if ("umap" %in% names(obj@reductions)) {
  p_umap <- DimPlot(obj, reduction = "umap", group.by = "analysis_group", label = TRUE) +
    ggtitle("UMAP of input Seurat object")
  ggsave(
    filename = file.path(config$output_dir, "tables", "input_umap_by_group.png"),
    plot = p_umap,
    width = 10,
    height = 8
  )
}

# =========================
# 5. Monocle3 拟时序分析
# =========================
message("开始进行 Monocle3 拟时序分析...")

# 将 Seurat 对象转换为 cell_data_set 对象
cds <- as.cell_data_set(obj)
colData(cds)$analysis_group <- obj$analysis_group

# 如果 Seurat 对象中已有 UMAP 坐标，则直接沿用
if ("umap" %in% names(obj@reductions)) {
  reducedDims(cds)$UMAP <- Embeddings(obj, reduction = "umap")
}

# 数据预处理、降维、聚类与轨迹学习
cds <- preprocess_cds(cds, num_dim = config$monocle_num_dim)
if (!"UMAP" %in% names(reducedDims(cds))) {
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
}
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = TRUE)

# 按细胞类型绘制轨迹图
trajectory_by_group <- plot_cells(
  cds,
  color_cells_by = "analysis_group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2
)

ggsave(
  filename = file.path(config$output_dir, "monocle3", "trajectory_by_group_before_ordering.png"),
  plot = trajectory_by_group,
  width = 10,
  height = 8
)

# 设定拟时序起点细胞
if (is.null(config$root_celltype)) {
  config$root_celltype <- names(group_table)[1]
  message(
    "未手动指定 root_celltype，暂时使用细胞数最多的群体：",
    config$root_celltype,
    "\n建议结合轨迹图和生物学意义再确认。"
  )
}

root_cells <- colnames(cds)[colData(cds)$analysis_group == config$root_celltype]
if (length(root_cells) == 0) {
  stop("未找到 root_celltype 对应的细胞，请检查 config$root_celltype。")
}
if (length(root_cells) > config$root_cells_n) {
  root_cells <- root_cells[seq_len(config$root_cells_n)]
}

# 对细胞进行拟时序排序
cds <- order_cells(cds, root_cells = root_cells)

# 按拟时序值绘图
trajectory_by_pseudotime <- plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2
)

ggsave(
  filename = file.path(config$output_dir, "monocle3", "trajectory_by_pseudotime.png"),
  plot = trajectory_by_pseudotime,
  width = 10,
  height = 8
)

# 再保存一张排序后的分组轨迹图
trajectory_by_group_after <- plot_cells(
  cds,
  color_cells_by = "analysis_group",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 2
)

ggsave(
  filename = file.path(config$output_dir, "monocle3", "trajectory_by_group_after_ordering.png"),
  plot = trajectory_by_group_after,
  width = 10,
  height = 8
)

# 导出每个细胞的拟时序值
pseudotime_df <- data.frame(
  cell = colnames(cds),
  analysis_group = as.character(colData(cds)$analysis_group),
  pseudotime = pseudotime(cds),
  stringsAsFactors = FALSE
)

write.csv(
  pseudotime_df,
  file = file.path(config$output_dir, "monocle3", "pseudotime_per_cell.csv"),
  row.names = FALSE
)

saveRDS(cds, file = file.path(config$output_dir, "monocle3", "monocle3_cds.rds"))

# =========================
# 6. CellChat 细胞通讯分析
# =========================
message("开始进行 CellChat 细胞通讯分析...")

future::plan("multisession", workers = config$future_workers)
on.exit(future::plan("sequential"), add = TRUE)

# 构建 CellChat 对象
cellchat <- createCellChat(object = obj, group.by = "analysis_group", assay = assay_use)

# 根据物种选择数据库
db_use <- switch(
  tolower(config$species),
  human = CellChatDB.human,
  mouse = CellChatDB.mouse,
  stop("config$species 只能填写 'human' 或 'mouse'。")
)

cellchat@DB <- db_use

# 标准 CellChat 分析流程
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = config$cellchat_min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

group_size <- as.numeric(table(cellchat@idents))

# 绘制细胞通讯数量网络图
save_base_plot(
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Number of interactions"
  ),
  filename = file.path(config$output_dir, "cellchat", "net_count_circle.pdf")
)

# 绘制细胞通讯强度网络图
save_base_plot(
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = "Interaction weights"
  ),
  filename = file.path(config$output_dir, "cellchat", "net_weight_circle.pdf")
)

# 绘制 outgoing / incoming signaling role 热图
outgoing_ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
incoming_ht <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

save_heatmap_plot(
  outgoing_ht,
  filename = file.path(config$output_dir, "cellchat", "outgoing_signaling_role_heatmap.pdf"),
  width = 10,
  height = 8
)

save_heatmap_plot(
  incoming_ht,
  filename = file.path(config$output_dir, "cellchat", "incoming_signaling_role_heatmap.pdf"),
  width = 10,
  height = 8
)

# 导出通讯结果表
communication_df <- subsetCommunication(cellchat)
write.csv(
  communication_df,
  file = file.path(config$output_dir, "cellchat", "all_communications.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(cellchat@net$count),
  file = file.path(config$output_dir, "cellchat", "interaction_count_matrix.csv"),
  row.names = TRUE
)

write.csv(
  as.data.frame(cellchat@net$weight),
  file = file.path(config$output_dir, "cellchat", "interaction_weight_matrix.csv"),
  row.names = TRUE
)

saveRDS(cellchat, file = file.path(config$output_dir, "cellchat", "cellchat_result.rds"))

# =========================
# 7. 保存环境信息
# =========================
writeLines(
  capture.output(sessionInfo()),
  con = file.path(config$output_dir, "sessionInfo.txt")
)
