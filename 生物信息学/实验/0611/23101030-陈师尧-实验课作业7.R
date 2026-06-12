library(Biostrings)
library(msa)
library(ape)
library(ggmsa)
library(ggplot2)
library(pwalign)

options(stringsAsFactors = FALSE)

output_dir <- "作业7输出"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}



seq1 <- DNAString("ACGTC")
seq2 <- DNAString("AATC")

score_matrix <- nucleotideSubstitutionMatrix(
  match = 10,
  mismatch = -3,
  baseOnly = TRUE
)

pair_aln <- pairwiseAlignment(
  pattern = seq1,
  subject = seq2,
  type = "global",
  substitutionMatrix = score_matrix,
  gapOpening = 5,
  gapExtension = 0
)

task1_lines <- c(
  "Task 1：双序列全局比对（Needleman-Wunsch）",
  paste0("序列1：", as.character(seq1)),
  paste0("序列2：", as.character(seq2)),
  paste0("比对得分：", score(pair_aln)),
  paste0("相似度（PID）：", round(pid(pair_aln), 2), "%"),
  paste0("比对后序列1：", as.character(alignedPattern(pair_aln))),
  paste0("比对后序列2：", as.character(alignedSubject(pair_aln)))
)

writeLines(task1_lines, con = file.path(output_dir, "task1_pairwise_alignment.txt"))
cat(paste(task1_lines, collapse = "\n"), "\n")


dna_file <- "20231215.fas"
dna_raw <- readDNAStringSet(dna_file)
dna_widths <- width(dna_raw)
modal_len <- as.integer(names(sort(table(dna_widths), decreasing = TRUE))[1])
dna_chars <- as.character(dna_raw)
names(dna_chars) <- names(dna_raw)
dna_check_lines <- c(
  "Task 2 DNA数据检查",
  paste0("原始序列条数：", length(dna_raw)),
  paste0("原始长度范围：", min(dna_widths), "-", max(dna_widths), " bp"),
  paste0("主长度（出现次数最多的长度）：", modal_len, " bp")
)

outlier_idx <- which(dna_widths != modal_len)
if (length(outlier_idx) > 0) {
  for (idx in outlier_idx) {
    seq_name <- names(dna_raw)[idx]
    seq_len <- dna_widths[idx]
    seq_chr <- dna_chars[idx]
    if (
      seq_len == 2 * modal_len &&
      substr(seq_chr, 1, modal_len) == substr(seq_chr, modal_len + 1, seq_len)
    ) {
      dna_chars[idx] <- substr(seq_chr, 1, modal_len)
      dna_check_lines <- c(
        dna_check_lines,
        paste0(
          "已自动修正：", seq_name,
          " 长度为 ", seq_len, " bp，",
          "检测到其前后两半完全相同，已截取前 ", modal_len, " bp 参与后续分析。"
        )
      )
    } else {
      stop(
        paste0(
          "检测到异常长度序列且无法自动修正：",
          seq_name, "（", seq_len, " bp）。请检查原始 FASTA 文件后重试。"
        )
      )
    }
  }
} else {
  dna_check_lines <- c(dna_check_lines, "未检测到异常长度序列。")
}

dna_seqs <- DNAStringSet(dna_chars)
names(dna_seqs) <- names(dna_raw)

cat("DNA序列条数：", length(dna_seqs), "\n")
cat("DNA序列长度范围：", min(width(dna_seqs)), "-", max(width(dna_seqs)), "bp\n")

writeLines(dna_check_lines, con = file.path(output_dir, "task2_dna_data_check.txt"))

dna_msa <- msa(dna_seqs, method = "Muscle", order = "input")
dna_aligned <- as(dna_msa, "DNAStringSet")

writeXStringSet(
  dna_aligned,
  filepath = file.path(output_dir, "task2_dna_aligned.fasta"),
  format = "fasta"
)

dna_bin <- as.DNAbin(as.matrix(dna_aligned))
dna_dist <- dist.dna(dna_bin, model = "K80")
dna_tree <- nj(dna_dist)
dna_tree$edge.length[abs(dna_tree$edge.length) < 1e-12] <- 0
dna_tree$edge.length <- pmax(dna_tree$edge.length, 0)

write.tree(dna_tree, file = file.path(output_dir, "task2_dna_nj_tree.nwk"))
tree_text <- readLines(file.path(output_dir, "task2_dna_nj_tree.nwk"), warn = FALSE)
tree_text <- gsub(":\\-0([,\\)])", ":0\\1", tree_text)
tree_text <- gsub(":\\-0;", ":0;", tree_text)
writeLines(tree_text, con = file.path(output_dir, "task2_dna_nj_tree.nwk"))
write.csv(
  as.matrix(dna_dist),
  file = file.path(output_dir, "task2_dna_distance_matrix.csv"),
  row.names = TRUE
)

dna_tree_ppt <- ladderize(dna_tree, right = TRUE)
png(
  filename = file.path(output_dir, "task2_dna_tree.png"),
  width = 1600,
  height = 1200,
  res = 150
)
plot(
  dna_tree_ppt,
  use.edge.length = FALSE,
  main = "Task 2 DNA Neighbor-Joining Tree",
  cex = 0.9
)
dev.off()

# ggmsa 图展示前 100 个碱基位点，便于在作业中观察保守区和变异位点
dna_plot_end <- min(100, width(dna_aligned)[1])
dna_msa_plot <- ggmsa(
  file.path(output_dir, "task2_dna_aligned.fasta"),
  start = 1,
  end = dna_plot_end,
  seq_name = TRUE,
  show.legend = FALSE,
  disagreement = FALSE,
  use_dot = FALSE
) +
  ggtitle("Task 2 DNA Multiple Sequence Alignment (First 100 bp)")

ggsave(
  filename = file.path(output_dir, "task2_dna_ggmsa.png"),
  plot = dna_msa_plot,
  width = 14,
  height = 7,
  dpi = 150
)


protein_file <- "euk.Acetyltransf_1.HG1.seqs.fasta"
protein_seqs <- readAAStringSet(protein_file)

cat("蛋白质序列条数：", length(protein_seqs), "\n")
cat("蛋白质序列长度范围：", min(width(protein_seqs)), "-", max(width(protein_seqs)), "aa\n")

protein_msa <- msa(protein_seqs, method = "Muscle", order = "input")
protein_aligned <- as(protein_msa, "AAStringSet")

writeXStringSet(
  protein_aligned,
  filepath = file.path(output_dir, "task3_protein_aligned.fasta"),
  format = "fasta"
)

protein_dist <- stringDist(protein_aligned, method = "hamming")
protein_hclust <- hclust(protein_dist, method = "average")

write.csv(
  as.matrix(protein_dist),
  file = file.path(output_dir, "task3_protein_distance_matrix.csv"),
  row.names = TRUE
)

png(
  filename = file.path(output_dir, "task3_protein_cluster.png"),
  width = 1600,
  height = 1200,
  res = 150
)
plot(
  protein_hclust,
  main = "Task 3 Protein Clustering Dendrogram",
  xlab = "Protein sequence",
  sub = "Distance method: Hamming distance on aligned amino acids",
  cex = 0.9
)
dev.off()

protein_plot_end <- width(protein_aligned)[1]
protein_msa_plot <- ggmsa(
  file.path(output_dir, "task3_protein_aligned.fasta"),
  start = 1,
  end = protein_plot_end,
  color = "Chemistry_AA",
  seq_name = TRUE,
  show.legend = TRUE
) +
  ggtitle("Task 3 Protein Multiple Sequence Alignment")

ggsave(
  filename = file.path(output_dir, "task3_protein_ggmsa.png"),
  plot = protein_msa_plot,
  width = 14,
  height = 6,
  dpi = 150
)
