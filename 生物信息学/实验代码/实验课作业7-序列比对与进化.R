library(Biostrings)
library(msa)
library(ape)
library(ips)
library(muscle)
library(ggmsa)
library(pwalign)

# Task 1：Needleman-Wunsch 双序列全局比对
substitution_matrix <- data.frame(
  A = c(10, -1, -3, -4),
  G = c(-1, 7, -5, -3),
  C = c(-3, -5, 9, 0),
  T = c(-4, -3, 0, 8)
)

rownames(substitution_matrix) <- c("A", "G", "C", "T")
colnames(substitution_matrix) <- c("A", "G", "C", "T")

seq1 <- c("A", "C", "G", "T", "C")
seq2 <- c("A", "A", "T", "C")

M <- length(seq1)
N <- length(seq2)
gap <- -5

score_matrix <- matrix(0, nrow = N + 1, ncol = M + 1)
rownames(score_matrix) <- c("0", seq2)
colnames(score_matrix) <- c("0", seq1)

score_matrix[1, ] <- seq(0, gap * M, by = gap)
score_matrix[, 1] <- seq(0, gap * N, by = gap)

for (i in 2:(N + 1)) {
  for (j in 2:(M + 1)) {
    base_seq2 <- rownames(score_matrix)[i]
    base_seq1 <- colnames(score_matrix)[j]

    diagonal_score <- score_matrix[i - 1, j - 1] +
      substitution_matrix[base_seq2, base_seq1]
    up_score <- score_matrix[i - 1, j] + gap
    left_score <- score_matrix[i, j - 1] + gap

    score_matrix[i, j] <- max(
      diagonal_score,
      up_score,
      left_score
    )
  }
}

score_matrix

dna_file <- "20231215.fas"

dna_seq <- readDNAStringSet(dna_file)
dna_aln <- msa(dna_seq)

dna_aln_seqinr <- msaConvert(
  dna_aln,
  type = "seqinr::alignment"
)

dna_aligned_seqs <- DNAStringSet(dna_aln_seqinr$seq)
names(dna_aligned_seqs) <- dna_aln_seqinr$nam

writeXStringSet(
  dna_aligned_seqs,
  filepath = "dna_aligned.fas"
)

dna_bin_raw <- fasta2DNAbin("dna_aligned.fas")
dna_bin_trimmed <- trimEnds(dna_bin_raw)

dist_raw <- dist.dna(dna_bin_raw, model = "K80")
dist_trimmed <- dist.dna(dna_bin_trimmed, model = "K80")

all.equal(
  as.matrix(dist_raw),
  as.matrix(dist_trimmed)
)

tree_trimmed <- nj(dist_trimmed)

plot(
  tree_trimmed,
  type = "phylogram",
  direction = "rightwards",
  use.edge.length = FALSE,
  cex = 0.7,
  label.offset = 0.3,
  no.margin = TRUE
)

# Task 2：DNA 多序列比对与系统发育树
dna_file <- "20231215.fas"

dna_seq <- readDNAStringSet(dna_file)
dna_aln <- msa(dna_seq)

dna_aln_seqinr <- msaConvert(
  dna_aln,
  type = "seqinr::alignment"
)

dna_aligned_seqs <- DNAStringSet(dna_aln_seqinr$seq)
names(dna_aligned_seqs) <- dna_aln_seqinr$nam

writeXStringSet(
  dna_aligned_seqs,
  filepath = "dna_aligned.fas"
)

dna_bin_raw <- fasta2DNAbin("dna_aligned.fas")
dna_bin_trimmed <- trimEnds(dna_bin_raw)

dist_raw <- dist.dna(dna_bin_raw, model = "K80")
dist_trimmed <- dist.dna(dna_bin_trimmed, model = "K80")

all.equal(
  as.matrix(dist_raw),
  as.matrix(dist_trimmed)
)

tree_trimmed <- nj(dist_trimmed)

plot(
  tree_trimmed,
  type = "phylogram",
  direction = "rightwards",
  use.edge.length = FALSE,
  cex = 0.7,
  label.offset = 0.3,
  no.margin = TRUE
)

# Task 3：蛋白质多序列比对与聚类/进化树
protein_file <- "outAAseq_trimmsf-1_NG1_seqs.fasta"

ggmsa(protein_file) +
  geom_seqlogo() +
  geom_msaBar()

protein_seq <- readAAStringSet(
  protein_file,
  format = "fasta"
)

protein_aln <- muscle(protein_seq)

protein_aln_trim <- maskGaps(
  protein_aln,
  min.fraction = 0.5,
  min.block.width = 4
)

protein_dist <- pwalign::stringDist(
  as(protein_aln_trim, "AAStringSet"),
  method = "hamming"
)

protein_clust <- hclust(
  protein_dist,
  method = "single"
)

plot(
  as.phylo(protein_clust),
  type = "phylogram",
  direction = "rightwards",
  use.edge.length = FALSE,
  cex = 0.8,
  label.offset = 0.2,
  no.margin = TRUE
)