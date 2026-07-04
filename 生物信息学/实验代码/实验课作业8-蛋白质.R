library(Biostrings)
library(Peptides)
library(protr)
library(bio3d)
library(r3dmol)
library(magrittr)
library(rstudioapi)

# Task1：蛋白质序列特征分析
protein_file <- "uniprotkb_accession_Q7W...fasta"

protein_sequence <- readAAStringSet(protein_file)
sequence_choose <- as.character(protein_sequence[[1]])

hydrophobicity(
  sequence_choose,
  scale = "KyteDoolittle"
)

mw(sequence_choose)

aIndex(sequence_choose)

pI(
  sequence_choose,
  pKscale = "EMBOSS"
)

extractAAC(sequence_choose)

dc <- extractDC(sequence_choose)
head(dc, n = 5)

data(AAdata)

autoCorrelation(
  sequence = sequence_choose,
  lag = 1,
  property = AAdata$Hydrophobicity$KyteDoolittle
)

protFP(sequence_choose)

blosumIndices(sequence_choose)

mswhimScores(sequence_choose)

vhseScales(sequence_choose)

# Task2：蛋白质结构分析
pdb <- read.pdb("1hel")

print(pdb)

modes <- nma(pdb)
plot(modes)

plot(
  modes,
  sse = pdb
)

plot.bio3d(
  pdb$atom$b[pdb$calpha],
  sse = pdb,
  typ = "l",
  ylab = "B-factor"
)

# Task3：蛋白质三维结构可视化
data(
  pdb_6zsl,
  package = "r3dmol"
)

r3dmol() %>%
  m_add_model(
    data = pdb_6zsl,
    format = "pdb"
  ) %>%
  m_zoom_to() %>%
  m_set_style(
    style = m_style_cartoon(
      color = "spectrum"
    )
  )