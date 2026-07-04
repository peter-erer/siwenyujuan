library(DESeq2)
library(factoextra)
library(pheatmap)

countData <- read.csv("chose_TCGAcount.csv")

rownames(countData) <- countData$Gene
countData <- countData[, -1]
countData <- countData[rowSums(countData > 0) > 0, ]

meta <- data.frame(
  id = colnames(countData),
  stringsAsFactors = FALSE
)

meta$type <- sapply(
  meta$id,
  function(x) strsplit(x, "_")[[1]][4]
)

meta$label <- ifelse(meta$type == "01A", "tumor", "pang")
meta$label <- factor(meta$label)
rownames(meta) <- meta$id

dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = meta,
  design = ~ label
)

dds <- DESeq(dds)
res <- results(dds, contrast = c("label", "tumor", "pang"))
res_df <- as.data.frame(res)

diff_gene_up <- subset(
  res_df,
  padj < 0.05 & log2FoldChange > 1
)

diff_gene_down <- subset(
  res_df,
  padj < 0.05 & log2FoldChange < -1
)

DEGS <- c(rownames(diff_gene_up), rownames(diff_gene_down))
data_choose <- countData[DEGS, ]

data_choosescale <- scale(t(data_choose))

d <- dist(data_choosescale)
fit1 <- hclust(d, method = "ward.D2")

plot(fit1, hang = -1, cex = 0.3, main = "clustering")

groups <- cutree(fit1, k = 2)
rect.hclust(fit1, k = 2, border = "red")

annotation_row <- data.frame(
  label = as.character(meta$label),
  stringsAsFactors = FALSE
)
rownames(annotation_row) <- meta$id

pheatmap::pheatmap(
  data_choosescale,
  annotation_row = annotation_row,
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  show_rownames = FALSE,
  main = "Heatmap of DEGs (rows: samples, columns: genes)",
  fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6
)