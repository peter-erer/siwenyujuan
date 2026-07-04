# 实验课7 R序列比对与进化 R代码整理

对应材料：

- `实验/课件/实验课7-R序列比对与进化-根据代码重新修改版-202606(1).pptx`
- `回放/生物信息学（协和班）第15周星期4第8,9节_笔记.txt`
- 老师板书整理

## 说明

- 这份整理以老师板书为主，结合课件和第15周实验课回放做了必要补全。
- 这节实验课分成三块：
  - `Needleman-Wunsch` 双序列全局比对
  - DNA 多序列比对与系统发育树
  - 蛋白质多序列比对与聚类/进化树
- 老师在回放里反复强调：这节课“特别适合理解考试逻辑”，尤其是第一部分的经典算法思想和双重循环填表过程。

## 一、老师对这节实验课的相关要求

老师在回放里对这节课的定位说得比较明确：

- 这节课属于“进化分析”和“序列比对”非常基础、也非常关键的一部分。
- 第一部分 `Needleman-Wunsch` 是经典算法，老师明确说“特别适合考试”。
- 虽然老师又补了一句“不一定直接这样考代码”，但明显是在提醒大家：
  - 这部分逻辑一定要理解
  - 公式、循环、打分矩阵、`max` 取值思路要会
- 后两部分 DNA/蛋白多序列比对，老师更强调“流程理解”和“格式转换 + 算距离 + 建树”的主链条。

老师还特别提醒了两类复习策略：

### 1. 不能只抄代码，要把过程理解清楚

回放里老师明确说：

- 有些代码可能有抄的成分没关系
- 但过程一定要理解清楚
- 特别是第一部分，如果不懂算法就写不出来

### 2. 考试应用题和简答题一定要写满

回放里老师专门提醒：

- 简答题要尽量写满
- 应用题分值高，能写的流程、逻辑、步骤都要写出来

所以这节课的复习重点不是“死背每一行函数名”，而是：

- 输入是什么
- 打分矩阵怎么用
- 动态规划怎么填表
- 多序列比对后怎样转格式、算距离、建树
- 蛋白和 DNA 流程的异同是什么

## 二、实验所需包：老师板书和课件附录

老师板书左上角和课件附录提到的包，核心是这些：

```r
library(Biostrings)
library(msa)
library(ape)
library(ips)
library(muscle)
library(ggmsa)
library(pwalign)
```

补充说明：

- 老师板书里还写了 `rstudioapi`，但这不是本实验主流程的核心包。
- 课件附录专门提醒：
  - 新版 `Biostrings` 里 `stringDist()` 已移到 `pwalign`
  - 所以应写成 `pwalign::stringDist()`
- 课件还提醒：
  - 在部分新环境里 `ggtree` 兼容性不好
  - 本实验统一用 `ape::plot()` 画树

## 三、Task 1：Needleman-Wunsch 双序列全局比对

### 1. 板书整理后的主干 R 代码

```r
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
```

### 2. 这段代码对应老师板书的哪几块

老师第一张板书几乎完整写了这个算法的骨架：

- 先构造 `substitution_matrix`
- 再定义：
  - `seq1`
  - `seq2`
  - `M`
  - `N`
  - `gap`
- 再初始化 `score_matrix`
- 再用两层循环逐格填表

板书右侧的核心公式就是：

```r
diagonal_score <- score_matrix[i - 1, j - 1] + substitution_matrix[base_seq2, base_seq1]
up_score <- score_matrix[i - 1, j] + gap
left_score <- score_matrix[i, j - 1] + gap

score_matrix[i, j] <- max(diagonal_score, up_score, left_score)
```

也就是课件第 7 页反复演示的那一套：

- 左上对角线：配对得分
- 从左：加 gap
- 从上：加 gap
- 取三者最大值

### 3. 这一部分最该掌握什么

老师回放里对这一段讲得最重：

- 这是最经典的算法
- 特别适合考试
- 核心就是：
  - 替换矩阵
  - gap 罚分
  - 初始化第一行第一列
  - 两层循环
  - 每格取 `max`

也就是说，考试时你最起码要会解释：

```text
为什么第一行第一列要按gap累加初始化
为什么每个格子要比较左上、上、左三个方向
为什么左上代表配对，上/左代表插入gap
```

### 4. 课件里讲了但板书没完全展开的：回溯

课件第 8 页讲了回溯思想，但老师板书这次主要写到“填分数矩阵”。

为了复习完整性，可以把回溯逻辑概括成这样：

```text
从右下角开始回溯：
1. 如果当前值来自左上 + 匹配/错配得分，就走左上
2. 如果来自上方 + gap，就向上走
3. 如果来自左方 + gap，就向左走
4. 一直走到左上角，得到最终比对结果
```

老师这节课更强调的是“动态规划填表”这一步，因为它最适合基础编程题和算法逻辑题。

## 四、Task 2：DNA 多序列比对与系统发育树

### 1. 板书整理后的主干 R 代码

```r
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
```

### 2. 这段代码对应老师板书的哪几块

老师第二张板书的主线非常清楚：

1. 读入 DNA `fasta`
2. 用 `msa()` 做多序列比对
3. `msaConvert()` 做格式转换
4. 写出对齐后的 `fasta`
5. 转成 `DNAbin`
6. `trimEnds()`
7. `dist.dna()`
8. `nj()`
9. `plot()`

也就是说，这部分的最短主链条就是：

```text
读fasta
-> msa多序列比对
-> 格式转换
-> 写出对齐结果
-> 转成DNAbin
-> trimEnds
-> dist.dna
-> nj建树
-> plot画树
```

### 3. 老师强调的“核心代码”在哪里

回放里老师对这部分说得很直接：

- 前面很多代码都是“格式转换”
- 真正核心的是：
  - 算距离
  - 建树
  - 可视化

所以你复习时别被代码长度吓住，要抓住：

- 多序列比对只是前处理
- 真正进化树的主线是：

```text
alignment -> distance -> tree
```

### 4. 为什么要 `trimEnds()`

板书里专门写了 `trimEnds()`，这个点值得记：

- 多序列比对后，两端往往会有很多 gap
- 这些区域可能对距离计算有干扰
- 所以先裁掉两端 gap 较多的部分，再算距离更稳妥

## 五、Task 3：蛋白质多序列比对与聚类/进化树

### 1. 板书整理后的主干 R 代码

```r
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
```

### 2. 这段代码对应老师板书的哪几块

老师第三张板书主线也很清楚：

1. `ggmsa()` 先看多序列比对的整体情况
2. `readAAStringSet()` 读蛋白序列
3. `muscle()` 做蛋白多序列比对
4. `maskGaps()` 去掉 gap 太多的区域
5. `pwalign::stringDist(..., method = "hamming")`
6. `hclust()`
7. `plot(as.phylo(...))`

### 3. 和 DNA 多序列比对的异同

这一部分最适合考试对比题：

相同点：

- 都是先做多序列比对
- 再算距离
- 再聚类/建树
- 最后可视化

不同点：

- DNA 这里主要用：
  - `readDNAStringSet`
  - `dist.dna`
  - `nj`
- 蛋白这里主要用：
  - `readAAStringSet`
  - `muscle`
  - `pwalign::stringDist`
  - `hclust`

老师回放里也提到：

- 蛋白这一部分其实也比较简单
- 多重比对后，把 gap 很多的区域去掉
- 有了距离就可以聚类

## 六、老师回放里对考试最有价值的提醒

### 1. 第一部分最适合考试

老师回放里明显重复了几次：

- `Needleman-Wunsch` 很经典
- 特别适合考试
- 核心就是公式 + 两层循环 + 打分矩阵

虽然老师又说“不一定直接考写这段代码”，但这其实是在提醒你：

- 至少会解释
- 至少会写伪代码
- 至少会把流程写出来

### 2. 简答题和应用题要写满

老师明确提醒：

- 分值高的题，能写多少写多少
- 流程、原理、步骤都尽量展开

所以这节课如果出应用题，最稳妥的写法不是一句话带过，而是要写出：

- 输入文件格式
- 用了什么包/函数
- 中间做了什么转换
- 输出是什么

### 3. 不要只盯着代码，要抓“核心步骤”和“格式转换”

老师反复强调：

- 大量代码只是为了做格式转换
- 真正核心的那几步其实不多

这节课可以直接总结成：

```text
Task1:
替换矩阵 -> 初始化 -> 双重循环填score matrix

Task2:
读DNA -> 多序列比对 -> 格式转换 -> 算距离 -> 建树

Task3:
读蛋白 -> 多序列比对 -> 去gap多区域 -> 算距离 -> 聚类/建树
```

## 七、这节实验课最该怎么复习

### 1. `Needleman-Wunsch`

最该背：

- `substitution_matrix`
- `gap`
- 第一行第一列初始化
- 两层循环
- `diagonal / up / left`
- `max()`

最该会解释：

- 为什么这是动态规划
- 为什么可以逐格累积最优子结构

### 2. DNA 多序列比对

最该背：

- `readDNAStringSet`
- `msa`
- `msaConvert`
- `fasta2DNAbin`
- `trimEnds`
- `dist.dna`
- `nj`

### 3. 蛋白多序列比对

最该背：

- `readAAStringSet`
- `muscle`
- `maskGaps`
- `pwalign::stringDist`
- `hclust`

## 八、如果考试让你写流程图或伪代码，最稳妥的写法

### 1. 双序列全局比对

```text
输入：两个序列seq1和seq2、替换矩阵、gap罚分

1. 初始化score matrix
2. 第一行和第一列按gap累加
3. 对矩阵每个内部位置：
   计算左上得分、上方得分、左方得分
4. 取三者最大值填入当前位置
5. 得到完整score matrix
6. 如需最终比对结果，再从右下角回溯
```

### 2. DNA 多序列比对与系统发育树

```text
输入：DNA fasta文件

1. 读取DNA序列
2. 做多序列比对
3. 将比对结果转换成合适格式
4. 去除两端gap较多区域
5. 计算序列间距离
6. 根据距离矩阵构建系统发育树
7. 绘制进化树
```

### 3. 蛋白多序列比对与聚类

```text
输入：蛋白质fasta文件

1. 读取蛋白质序列
2. 做多序列比对
3. 去除gap较多区域
4. 计算蛋白序列间距离
5. 根据距离做聚类
6. 将聚类结果画成树状图/进化树
```

## 九、一句提醒

这节实验课最重要的不是把所有函数名硬背下来，而是你要会解释：

- `Needleman-Wunsch` 为什么用双重循环和 `max`
- 多序列比对之后为什么还要做格式转换
- 为什么算距离之后才能建树
- DNA 流程和蛋白流程的异同是什么
- 哪些步骤是“格式准备”，哪些步骤才是真正核心计算
