if (!require("Biostrings")) {
  install.packages("BiocManager")
  BiocManager::install("Biostrings")
}
library(Biostrings)
file_path <- "DNA_sequence_new.fasta"
lines <- readLines(file_path)
seq_lines <- lines[!grepl("^>",lines)]
dna_seq <- toupper(paste(seq_lines, collapse = ""))
cat("序列长度：",nchar(dna_seq),"\n\n")

expand_inr_motifs <- function() {
  pos1 <- c("C", "T")
  pos2 <- c("C", "T")
  pos3 <- c("A")
  pos4 <- c("A", "T", "C", "G")
  pos5 <- c("A", "T")
  pos6 <- c("C", "T")
  pos7 <- c("C", "T")
  combos <- expand.grid(
    pos1, pos2, pos3, pos4, pos5, pos6, pos7, 
    stringsAsFactors = FALSE
  )
  apply(combos, 1, paste0, collapse="")
}

inr_128 <- expand_inr_motifs()
print(length(inr_128))  # 结果会是 128
print(head(inr_128))    # 看看前几个长什么样

predict_TSS <- function(seq){
  s <- DNAString(seq)
  tss_candidates <- data.frame()
  motifs <- list(
    TATA = c("TATAAAA", "TATAAAT", "TATATAA", "TATATAT"),
    Inr  = expand_inr_motifs(),
    BRE  = c("GGAGCC", "GGTGCC", "GCAGCC", "GCTGCC"),
    CAAT = c("GGCCAATCT"),
    GC   = c("GGGCGG")
  )
  for (name in names(motifs)) {
    for (mot in motifs[[name]]) {
      matches <- matchPattern(DNAString(mot), s)
      if (length(matches) > 0) {
        tss_candidates <- rbind(
          tss_candidates,
          data.frame(
            Motif = name,
            Motif_seq = mot,
            Start = start(matches),
            End = end(matches),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  return(tss_candidates)
}

predict_TSE <- function(seq){
  s <- DNAString(seq)
  tse_candidates <- data.frame()
  polyA__motifs <- c("AATAAA","ATTAAA","AGTAAA","TATAAA")
  for (mot in polyA__motifs) {
    matches <- matchPattern(DNAString(mot), s)
    if (length(matches) > 0) {
      tse_candidates <- rbind(
        tse_candidates,
        data.frame(
          Motif = "polyA",
          Motif_seq = mot,
          Start = start(matches),
          End = end(matches),
          stringsAsFactors = FALSE
      )
    )
    }
  }
  return(tse_candidates)
}
tss_all <- predict_TSS(dna_seq)
tse_all <- predict_TSE(dna_seq)

true_tss <- 2018
true_tse <- 194629

if (!is.null(tss_all) && nrow(tss_all) > 0) {
  tss_all$Predicted_TSS <- NA
  tss_all$Predicted_TSS[tss_all$Motif == "TATA"] <- tss_all$End[tss_all$Motif == "TATA"] + 25
  tss_all$Predicted_TSS[tss_all$Motif == "BRE"]  <- tss_all$End[tss_all$Motif == "BRE"] + 35
  tss_all$Predicted_TSS[tss_all$Motif == "CAAT"] <- tss_all$End[tss_all$Motif == "CAAT"] + 100
  tss_all$Predicted_TSS[tss_all$Motif == "GC"]   <- tss_all$End[tss_all$Motif == "GC"] + 80
  tss_all$Predicted_TSS[tss_all$Motif == "Inr"]  <- floor((tss_all$Start[tss_all$Motif == "Inr"] + tss_all$End[tss_all$Motif == "Inr"]) / 2)
  
  tss_all$Distance <- abs(tss_all$Predicted_TSS - true_tss)
  tss_sorted <- tss_all[order(tss_all$Distance), ]
  
  cat("\n=== TSS (转录起始位点) 预测结果 ===\n")
  cat("TSS 候选位点总数:", nrow(tss_all), "\n")
  cat("距离真实 TSS (2018) 最近的 Top 10 候选位点:\n")
  print(head(tss_sorted, 10))
}

if (!is.null(tse_all) && nrow(tse_all) > 0) {
  tse_all$Predicted_TSE <- tse_all$End + 20
  tse_all$Distance <- abs(tse_all$Predicted_TSE - true_tse)

  tse_sorted <- tse_all[order(tse_all$Distance), ]
  
  cat("\n=== TSE (转录终止位点) 预测结果 ===\n")
  cat("TSE 候选位点总数:", nrow(tse_all), "\n")
  cat("距离真实 TSE (194629) 最近的 Top 10 候选位点:\n")
  print(head(tse_sorted, 10))
}