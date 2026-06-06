if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "factoextra"))
install.packages("pheatmap")

library(DESeq2)
library(factoextra)
library(pheatmap)

countData <- read.csv("choose_TCGAcount.csv", stringsAsFactors = FALSE)
rownames(countData) <- countData$Gene
countData <- countData[ , -1] 
keep <- rowSums(countData > 0) >= 1
countData <- countData[keep, ]
sample_ids <- colnames(countData)
group_label <- ifelse(grepl("01A$", sample_ids), "tumor", 
                      ifelse(grepl("11A$", sample_ids), "pang", "other"))

meta <- data.frame(
  row.names = sample_ids,
  label = factor(group_label, levels = c("pang", "tumor")) 
)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = meta,
                              design = ~ label)

dds <- DESeq(dds)
res <- results(dds, contrast = c("label", "tumor", "pang"))
summary(res)

res_clean <- na.omit(res)
up_genes <- rownames(res_clean[res_clean$padj < 0.05 & res_clean$log2FoldChange > 1, ])
down_genes <- rownames(res_clean[res_clean$padj < 0.05 & res_clean$log2FoldChange < -1, ])

DEGs <- c(up_genes, down_genes)
cat("找到的上调基因数：", length(up_genes), "\n")
cat("找到的下调基因数：", length(down_genes), "\n")


deg_matrix <- countData[DEGs, ]
log_deg_matrix <- log2(deg_matrix + 1)
data_choosescale <- scale(t(log_deg_matrix))

d <- dist(data_choosescale, method = "euclidean")
fit1 <- hclust(d, method = "ward.D2")

plot(fit1, hang = -1, cex = 0.5, main = "Clustering of Samples based on DEGs")

annotation_row <- data.frame(Group = meta$label)
rownames(annotation_row) <- rownames(data_choosescale)

pheatmap::pheatmap(
  data_choosescale,
  annotation_row    = annotation_row,
  clustering_method = "ward.D2",
  show_colnames     = FALSE,
  show_rownames     = FALSE,
  main              = "Heatmap of DEGs (rows: samples, columns: genes)",
  fontsize          = 6,
  fontsize_row      = 6,
  fontsize_col      = 6
)