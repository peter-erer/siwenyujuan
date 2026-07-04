# 实验课8 R蛋白 R代码整理

对应材料：

- `实验/课件/实验课8-R蛋白-更新.pptx`
- `回放/生物信息学（协和班）第16周星期4第8,9节_笔记.txt`
- 老师板书整理

## 说明

- 这份整理以老师板书为主，结合课件和第16周实验课回放做了必要补全。
- 这节实验课主要分成两块：
  - `Task1` 蛋白质序列特征分析
  - `Task2` 蛋白质结构分析与三维结构可视化
- 老师在回放里明确说，这节课整体是“最容易的一堂实验课”之一，但不代表可以不复习；它对应的是蛋白序列和结构分析的常用基础方法。

## 一、老师对这节实验课的相关要求

老师在回放里强调了几件很关键的事：

- 八堂实验课学过的步骤、流程、代码思路都要特别熟。
- 不只是会跑代码，还要会：
  - 画流程图
  - 用语言描述完整步骤
  - 解释每一步输入、处理和输出
- 这节课的代码本身不算难，重点是理解“一个蛋白序列读进来之后，可以做哪些常见特征分析”和“拿到 PDB 后可以怎样做结构分析和可视化”。

老师对这节课的定位大致是：

- 很实用
- 偏基础
- 适合做批处理
- 适合作为实验课方法题、流程题或应用题素材

另外老师还特别讲了两点：

### 1. 为什么还要写代码，不直接用在线工具

老师回放里明确说：

- 如果只分析一两条蛋白序列，在线工具也能做
- 但代码的优势是方便批处理
- 如果以后要处理成百上千甚至更多序列，写代码更实用

### 2. 这节课的重点不是“背全所有高级特征名”

老师说得很清楚：

- 课上列的是常用例子
- 真正做研究时，序列特征可以扩展出上千、上万种
- 所以考试和复习更重要的是知道：
  - 先读序列
  - 再看组成
  - 再看理化性质
  - 再看高级描述符
  - 再进入结构分析和 3D 可视化

## 二、实验所需包：老师板书和课件附录

老师板书最上面写了这些包，课件附录也对应列出来了：

```r
library(Biostrings)
library(Peptides)
library(protr)
library(bio3d)
library(r3dmol)
library(magrittr)
library(rstudioapi)
```

课件附录还专门提醒：

- `readAAStringSet()` 来自 `Biostrings`
- `%>%` 来自 `magrittr`
- `pdb_6zsl` 是 `r3dmol` 的内置示例数据，不需要自己下载

## 三、Task1：蛋白质序列特征分析

### 1. 板书整理后的主干 R 代码

```r
library(Biostrings)
library(Peptides)
library(protr)

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
```

## 四、Task1 代码对应老师板书的哪几块

### 1. 先读蛋白序列，再选一条序列出来

对应第一张板书左上：

```r
protein_file <- "uniprotkb_accession_...fasta"
protein_sequence <- readAAStringSet(protein_file)
sequence_choose <- as.character(protein_sequence[[1]])
```

老师回放里讲得很清楚：

- 第一步就是把序列读进去
- 然后选一条出来分析

也就是说，这节课最前面的输入非常简单：

- 一个蛋白质 `fasta` 文件
- 读入后取其中一条序列做演示

### 2. 常见理化性质

对应第一张板书左下：

```r
hydrophobicity(sequence_choose, scale = "KyteDoolittle")
mw(sequence_choose)
aIndex(sequence_choose)
pI(sequence_choose, pKscale = "EMBOSS")
```

这些分别对应老师回放里提到的常见指标：

- 疏水性指数
- 分子量 `MW`
- 脂肪指数 `A index`
- 等电点 `pI`

这部分是最基础、最容易理解的一块。

### 3. 氨基酸组成和二肽组成

对应第一张板书右上：

```r
extractAAC(sequence_choose)

dc <- extractDC(sequence_choose)
head(dc, n = 5)
```

这里的逻辑是：

- `extractAAC()`：氨基酸组成
- `extractDC()`：二肽组成

老师回放里说得很直接：

- 先看组成
- 再看理化性质

这块很适合出“某类函数是干什么的”这种基础题。

### 4. 更高级一点的描述符

对应第一张板书中下和右下：

```r
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
```

老师回放里对这部分的定位是：

- 这只是举例
- 真正做研究时，高级特征远不止这些
- 课上主要是让大家知道“蛋白序列特征可以往很深很广扩展”

所以考试复习时，这块不一定要求你死背每一个名字，但要会说：

- 这类函数是在提取更高级的数值化特征
- 后续可用于机器学习、聚类、分类或批量分析

## 五、Task1 最该掌握什么

如果按考试思路复习，这部分最该抓的是这条主线：

```text
读入蛋白序列
-> 选一条序列
-> 看氨基酸组成
-> 看理化性质
-> 看高级描述符
-> 输出一组可以继续分析的数值特征
```

老师回放里实际上强调的是：

- 输入一个蛋白序列
- 不同函数只是从不同角度提取特征
- 这些特征以后可以批量算

## 六、Task2：蛋白质结构分析

### 1. 板书整理后的主干 R 代码

```r
library(bio3d)

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
```

## 七、Task2 代码对应老师板书的哪几块

### 1. 读入 PDB 结构

对应第二张板书左上：

```r
pdb <- read.pdb("1hel")
print(pdb)
```

老师回放里说得很清楚：

- 结构分析一般需要一个 `PDB` 文件
- 如果知道 `PDB ID`，可以直接读入

这里的 `1hel` 就是老师举的示例。

### 2. 正常模态分析

对应第二张板书左中：

```r
modes <- nma(pdb)
plot(modes)
plot(modes, sse = pdb)
```

这部分是结构层面的分析：

- `nma()` 做正常模态分析
- 再通过 `plot()` 看结果

老师这节课没有把这一块讲得特别深，更多是在展示“PDB 读进来以后可以做什么分析”。

### 3. B-factor 波动图

对应第二张板书左下：

```r
plot.bio3d(
  pdb$atom$b[pdb$calpha],
  sse = pdb,
  typ = "l",
  ylab = "B-factor"
)
```

老师回放里明确点了这个：

- 结构分析部分要看一个典型的 `B-Factor` 波动图

所以这块最该会解释的是：

- `B-Factor` 反映局部原子/残基的波动或柔性
- 曲线不同位置高低不同，说明不同区域稳定性不同

## 八、Task3：蛋白质三维结构可视化

### 1. 板书整理后的主干 R 代码

```r
library(r3dmol)
library(magrittr)

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
```

## 九、Task3 代码对应老师板书的哪几块

对应第二张板书右侧：

```r
data(
  pdb_6zsl,
  package = "r3dmol"
)

r3dmol() %>%
  m_add_model(data = pdb_6zsl, format = "pdb") %>%
  m_zoom_to() %>%
  m_set_style(style = m_style_cartoon(color = "spectrum"))
```

这块逻辑很直观：

1. 先加载内置 `PDB` 数据
2. 创建 3D 可视化对象
3. 加入模型
4. 缩放到合适视角
5. 用 cartoon 样式和 `spectrum` 配色显示

老师回放里提到：

- `R` 里可以直接做 3D 可视化
- 还能用鼠标拖动、缩放、旋转
- 但如果做更专业的结构展示，他自己更常用 `PyMOL`

所以这节课对考试最重要的，不是 3D 界面怎么操作，而是你要知道：

- `PDB` 结构读进来后可以在 `R` 里直接做基础三维展示

## 十、老师回放里对考试最有价值的提醒

### 1. 实验课所有流程都要特别熟

老师这节回放里其实是在最后总复习，所以他说得很直接：

- 八堂实验课的步骤、流程、代码思路都要熟
- 不只是会写，还要会讲
- 也要会画流程图

### 2. 应用题很可能要求“文字描述 + 流程图”

老师明确说：

- 应用题里要写完整描述性的语言
- 还要画标准流程图

所以这节课如果出应用题，最稳妥的写法就是：

- 输入什么
- 先做什么
- 再做什么
- 最后输出什么结果

### 3. 这节课代码很适合说明“为什么实验课用代码”

老师强调：

- 单条序列在线做也行
- 但写代码适合批处理

这其实很适合写进应用题答案里，因为它能解释“为什么不用网站，而要写 `R`”。

## 十一、这节实验课最该怎么复习

### 1. 蛋白序列特征分析

最该背的不是所有函数名，而是这条主线：

```text
读入蛋白序列
-> 选定一条序列
-> 分析氨基酸组成
-> 分析理化性质
-> 提取高级描述符
```

### 2. 蛋白结构分析

最该背的是：

```text
读入PDB
-> 做结构分析
-> 画B-factor等结构图
-> 进行三维可视化
```

### 3. 这节课和考试的关系

这节课本身不一定是最核心的大题来源，但它和老师总复习要求高度一致的地方在于：

- 实验课代码思路
- 输入到输出的完整流程
- 用语言描述每一步
- 流程图规范表达

## 十二、如果考试让你写流程图或伪代码，最稳妥的写法

### 1. 蛋白序列特征分析

```text
输入：蛋白质fasta序列

1. 读取蛋白质序列
2. 选取需要分析的目标序列
3. 计算氨基酸组成与二肽组成
4. 计算理化性质，如疏水性、分子量、脂肪指数、等电点
5. 进一步提取高级描述符
6. 输出一组蛋白质序列特征
```

### 2. 蛋白质结构分析与三维可视化

```text
输入：PDB文件或PDB编号

1. 读取蛋白质三维结构
2. 对结构做基础分析
3. 计算或展示B-factor等结构特征
4. 构建三维可视化对象
5. 设置显示样式并进行可视化
6. 输出结构分析结果和三维图像
```

## 十三、一句提醒

这节实验课最重要的不是把所有特征名都背全，而是你要会解释：

- 蛋白序列读进来后可以提取哪些常见特征
- 为什么写代码比在线工具更适合批处理
- `PDB` 文件能支持哪些基础结构分析
- `B-factor` 图和三维可视化各自在看什么
