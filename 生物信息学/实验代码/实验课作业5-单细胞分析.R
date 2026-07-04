library(Seurat)
library(dplyr)
library(clustree)
library(ggplot2)
library(ggraph)

# 1. 读入 10X 数据并构建 Seurat 对象
pbmc.data <- Read10X(data.dir = "./filtered_feature_bc_matrix/")

pbmc <- CreateSeuratObject(
  counts = pbmc.data,
  project = "pbmc3k",
  min.cells = 3,
  min.features = 200
)

# 2.1 计算线粒体比例，先看过滤前 QC
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

p1 <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

p1

# 2.2 过滤低质量/异常细胞，再看过滤后 QC
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 3000 &
    percent.mt < 10
)

p2 <- VlnPlot(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

p2

# 3. 数据归一化 + 高变基因筛选
pbmc <- NormalizeData(
  pbmc,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

pbmc <- FindVariableFeatures(
  pbmc,
  selection.method = "vst",
  nfeatures = 2000
)

top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot1
plot2

# 4. 标准化 + PCA
all.genes <- rownames(pbmc)

pbmc <- ScaleData(
  pbmc,
  features = all.genes
)

pbmc <- RunPCA(
  pbmc,
  features = VariableFeatures(object = pbmc)
)

# 5.1 PCA 可视化与 ElbowPlot
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimHeatmap(
  pbmc,
  dims = 1:15,
  cells = 500,
  balanced = TRUE
)

ElbowPlot(pbmc, ndims = 30)

# 5.2 JackStraw 评估主成分
pbmc <- JackStraw(
  pbmc,
  num.replicate = 100,
  dims = 30
)

pbmc <- ScoreJackStraw(
  pbmc,
  dims = 1:30
)

JackStrawPlot(pbmc, dims = 1:30)

# 5.3 也可以粗略看累计方差贡献
pc_sd <- Stdev(pbmc, reduction = "pca")
cumulative_var <- cumsum(pc_sd^2) / sum(pc_sd^2)

cumulative_var[3]
cumulative_var[20]

# 实验里常直接选前 20 个 PC
dims_used <- 1:20

# 6. 聚类 + UMAP/TSNE
res.used <- seq(0.1, 1, by = 0.1)

pbmc <- pbmc %>%
  FindNeighbors(dims = dims_used) %>%
  FindClusters(resolution = res.used) %>%
  RunUMAP(dims = dims_used) %>%
  RunTSNE(dims = dims_used)

clas.tree.out <- clustree::clustree(pbmc) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Set1") +
  scale_edge_color_continuous(
    low = "grey80",
    high = "red"
  )

print(clas.tree.out)

sel.clust <- "RNA_snn_res.0.3"
pbmc <- SetIdent(pbmc, value = sel.clust)
pbmc$seurat_clusters <- pbmc@meta.data[[sel.clust]]

# 6.3 聚类结果可视化
DimPlot(
  pbmc,
  reduction = "umap",
  label = TRUE,
  label.size = 5
)

DimPlot(
  pbmc,
  reduction = "tsne",
  label = TRUE,
  label.size = 5
)

# 7. 用 marker 基因染色做细胞注释
FeaturePlot(
  object = pbmc,
  features = c("PTPRC"),
  cols = c("gray", "blue"),
  max.cutoff = 2,
  min.cutoff = 0
)

# 8. 找各 cluster 的 marker 基因
pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE
)

pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

top10 <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_head(n = 10)

DoHeatmap(
  pbmc,
  features = top10$gene
) +
  NoLegend()