library(Seurat)

# 1. 读取 10x 数据并构建 Seurat 对象
data_dir <- "filtered_feature_bc_matrix"
counts <- Read10X(data.dir = data_dir)

pbmc <- CreateSeuratObject(
  counts = counts,
  project = "scRNA",
  min.cells = 3,
  min.features = 200
)

cat("原始数据维度：\n")
print(dim(pbmc))
cat("原始细胞数：", ncol(pbmc), "\n")
cat("原始基因数：", nrow(pbmc), "\n")

# 查看 Seurat 对象中的原始表达矩阵
cat("\ncounts 矩阵预览：\n")
print(GetAssayData(pbmc, assay = "RNA", layer = "counts")[1:5, 1:5])

# 2. 质控（QC）
# 计算线粒体基因比例
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# 质控前可视化
png("01_QC_before_vlnplot.png", width = 1200, height = 800, res = 150)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

png("02_QC_scatter_before.png", width = 1200, height = 600, res = 150)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot1 + plot2)
dev.off()

# 根据课件条件过滤细胞
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 10
)

cat("\n过滤后数据维度：\n")
print(dim(pbmc))
cat("过滤后细胞数：", ncol(pbmc), "\n")
cat("过滤后基因数：", nrow(pbmc), "\n")

# 质控后可视化
png("03_QC_after_vlnplot.png", width = 1200, height = 800, res = 150)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# 3. 数据标准化与高变基因筛选

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
cat("\n前 10 个高变基因：\n")
print(top10)

png("04_variable_features.png", width = 1200, height = 800, res = 150)
plot_vf <- VariableFeaturePlot(pbmc)
plot_labeled <- LabelPoints(plot = plot_vf, points = top10, repel = TRUE)
print(plot_labeled)
dev.off()


# 4. 数据缩放、PCA 降维、聚类与 UMAP
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

cat("\nPCA 结果：\n")
print(pbmc[["pca"]])

png("05_elbowplot.png", width = 1000, height = 800, res = 150)
ElbowPlot(pbmc, ndims = 30)
dev.off()

# 根据课件建议取前 20 个主成分
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.3)
pbmc <- RunUMAP(pbmc, dims = 1:20)

cat("\n聚类结果（各 cluster 细胞数）：\n")
print(table(pbmc$seurat_clusters))

png("06_umap_clusters.png", width = 1000, height = 800, res = 150)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5)
dev.off()

# 5. Marker 基因表达可视化辅助注释
marker_genes <- c(
  "PTPRC", "CD3D", "CD3E", "CD8A", "CD4",
  "MS4A1", "CD79A", "GNLY", "NKG7",
  "LYZ", "CD68", "EPCAM", "KRT18", "KRT19",
  "PECAM1", "VWF", "DCN", "COL1A1",
  "CLDN18", "ABCA3", "AGER", "SFTPA1"
)

marker_genes_present <- marker_genes[marker_genes %in% rownames(pbmc)]
cat("\n数据中存在的 marker 基因：\n")
print(marker_genes_present)

if (length(marker_genes_present) > 0) {
  png("07_featureplot_markers.png", width = 1800, height = 1400, res = 150)
  print(FeaturePlot(
    object = pbmc,
    features = marker_genes_present,
    cols = c("gray", "blue"),
    min.cutoff = 0,
    max.cutoff = 2,
    ncol = 4
  ))
  dev.off()
}

# 6. SingleR 自动注释（可选）
if (requireNamespace("SingleR", quietly = TRUE) &&
    requireNamespace("celldex", quietly = TRUE) &&
    requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  library(SingleR)
  library(celldex)

  hpca.se <- celldex::HumanPrimaryCellAtlasData()
  clusters <- pbmc@meta.data$seurat_clusters
  expr_data <- GetAssayData(pbmc, assay = "RNA", layer = "data")

  pred.hesc <- SingleR(
    test = expr_data,
    ref = hpca.se,
    labels = hpca.se$label.main,
    clusters = clusters
  )

  cat("\nSingleR 注释结果：\n")
  print(table(pred.hesc$labels))

  celltype <- data.frame(
    ClusterID = rownames(pred.hesc),
    celltype = pred.hesc$labels,
    stringsAsFactors = FALSE
  )

  pbmc@meta.data$singleR <- celltype[
    match(clusters, celltype$ClusterID),
    "celltype"
  ]

  png("08_umap_singleR.png", width = 1200, height = 900, res = 150)
  DimPlot(
    pbmc,
    reduction = "umap",
    group.by = "singleR",
    label = TRUE,
    label.size = 4,
    pt.size = 0.5
  )
  dev.off()
} else {
  cat("\n未检测到 SingleR/celldex 相关包，已跳过自动注释步骤。\n")
}

# 7. 差异基因分析：寻找各 cluster 的 marker 基因
all_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(all_markers, "all_cluster_markers.csv", row.names = FALSE)

cat("\n各 cluster 前几行 marker 基因结果：\n")
print(head(all_markers))

top_markers <- all_markers |>
  dplyr::group_by(cluster) |>
  dplyr::slice_max(order_by = avg_log2FC, n = 10) |>
  dplyr::ungroup()

write.csv(top_markers, "top10_markers_each_cluster.csv", row.names = FALSE)

top10_genes <- unique(top_markers$gene)
top10_genes <- top10_genes[top10_genes %in% rownames(pbmc)]

if (length(top10_genes) > 1) {
  png("09_heatmap_top_markers.png", width = 1600, height = 1200, res = 150)
  print(DoHeatmap(pbmc, features = top10_genes, size = 3) + NoLegend())
  dev.off()
}


pbmc$new.cluster.ids <- as.character(pbmc$seurat_clusters)


png("10_umap_manual_annotation_template.png", width = 1000, height = 800, res = 150)
DimPlot(pbmc, reduction = "umap", group.by = "new.cluster.ids", label = TRUE, pt.size = 0.5)
dev.off()

# 9. 保存结果
saveRDS(pbmc, file = "pbmc_seurat_analysis.rds")

cat("\n分析完成。\n")
cat("已输出的主要结果包括：\n")
cat("1. 质控图、UMAP 图、marker 表达图、热图\n")
cat("2. all_cluster_markers.csv\n")
cat("3. top10_markers_each_cluster.csv\n")
cat("4. pbmc_seurat_analysis.rds\n")
