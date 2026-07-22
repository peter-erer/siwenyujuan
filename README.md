# 🐟 四文鱼卷：一个知识花园

> 后之览者，亦将有感于**斯文**。
> ——王羲之《兰亭集序》

## 🌊 从“三文鱼卷”开始的旅途

**“三文鱼卷”**是由北航10系学长创立的微信公众号（[最后一篇招新推送中含有简要介绍](https://mp.weixin.qq.com/s/q_XouH0q9dekRLP-RRbD4g)），曾经是覆盖航类、计类、医工班等学院的学习资源分享平台，帮助过无数伙伴（包括笔者），但是在2024年受不可抗力停止更新。（[最后一篇推送](https://mp.weixin.qq.com/s/eyuAGtSxbU5LPHF02o3x_Q)）我们希望能够把学习与分享的火种传递下去，于是便有了这个仓库——“四文鱼卷”。~~（东晋著名书法家王羲之曾在《兰亭集序》中公开表达其对**四文**鱼卷的赞美之情，令人动容。）~~ 

它主要面向**北航医工交叉学科群**（生物与医学工程学院10系、医学科学与工程学院45系，包括生物医学工程、智能医学工程等专业），提供一个分享学习经验、保留学习资料的平台。

它保持三文鱼卷的初心不变，但是转移了阵地。为什么不继续使用微信公众号？因为相比于公众号只能由少数人管理， [GitHub](https://github.com/) 提供了**开放协作、去中心化、开源共享**的天然土壤。这意味着，即使我们这群**老登**毕业之后，任何**小登**都可以继续为这个仓库贡献资料，传递下去。

## ❤ 贡献指南

欢迎贡献！无论是一份考试回忆、一处错别字，还是一门课程的完整资料，都可以成为这座知识花园的一部分。提交前请先阅读下方的 [仓库原则与收录指南](#📜-仓库原则与收录指南)，确认资料可以公开分享。

### 资料怎样放置

- 补充已有课程：把资料放入对应的课程文件夹，并尽量使用清晰的文件名，例如 `2025年期末考试回忆.md`、`第三章复习笔记.pdf`。
- 新增课程：在仓库根目录新建以课程名称命名的文件夹，并在其中创建 `README.md`，简要介绍课程和所收录的资料。
- Markdown 文档中的图片建议放在同一课程目录的 `assets` 或 `images` 文件夹中，并使用相对路径引用。
- 如果新增了一门课程，请同时把它补充到本 README 的“内容概览”表格中。
- 不需要修改 `site-config`、GitHub Actions 或生成的网站文件；资料合并到 `main` 分支后，网站会自动重新构建并部署。

### 方法一：提交 Issue（最简单）

如果还不熟悉 Git，可以直接[新建 Issue](https://github.com/peter-erer/siwenyujuan/issues/new)，说明课程名称、资料类型和内容来源，并上传附件或提供可访问的文件链接。仓库维护者会协助检查和整理。

### 方法二：在 GitHub 网页上提交 Pull Request

1. 打开[四文鱼卷仓库](https://github.com/peter-erer/siwenyujuan)，点击右上角的 **Fork**，将仓库复制到自己的账号下。
2. 进入自己 Fork 后的仓库，打开对应课程目录；如需新建课程，可通过 **Add file → Create new file** 创建 `课程名称/README.md`。
3. 使用 **Add file → Upload files** 上传资料，并按需编辑课程目录中的 `README.md`。
4. 点击 **Commit changes** 保存修改，提交说明可写成“补充××课程考试回忆”或“修正××课程资料链接”。
5. 返回仓库页面，点击 **Contribute → Open pull request**，简要说明新增或修改的内容，然后提交 Pull Request。

建议一次 Pull Request 只处理一门课程或一类相关修改，这样更便于审核。如果维护者提出调整建议，直接在同一分支继续修改即可，Pull Request 会自动更新。

### 方法三：使用 Git 在本地提交

先在 GitHub 上 Fork 本仓库，再运行以下命令；请把 `<你的用户名>` 和分支名替换成自己的信息：

```bash
git clone https://github.com/<你的用户名>/siwenyujuan.git
cd siwenyujuan
git switch -c add-course-materials

# 将资料复制到对应课程目录并完成必要的 README 修改后
git add -- "课程名称/" README.md
git commit -m "补充××课程资料"
git push -u origin add-course-materials
```

推送完成后，按照 GitHub 页面提示创建 Pull Request 即可。关于 [Git](https://git-scm.com/) 与 [GitHub](https://github.com/) 的进一步用法，推荐阅读[廖雪峰老师的 Git 教程](https://liaoxuefeng.com/books/git/)。当然，也可以把仓库链接、准备贡献的资料和具体需求交给 Codex 等 AI Agent，请它协助整理文件并创建 Pull Request。

## 📂 内容概览

在网站中可使用顶部搜索框检索全部 Markdown 笔记；PDF、Word、TXT 与相关代码可从各课程首页直接打开或下载。目前内容主要是协和医班的课程，欢迎生医、智医、大类同学积极补充。

| 学期 | 已收录课程 |
| --- | --- |
| **大一下** | [中国近现代史纲要](./中国近现代史纲要/README.md) |
| **大二下** | [计算机组成与编程基础](./计算机组成与编程基础/README.md) · [细胞生物学](./细胞生物学/README.md) · [生物医学电子](./生物医学电子/README.md) · [军事理论](./军事理论/README.md)  |
| **大三上** | [生物力学设计与仿真](./生物力学设计与仿真/README.md) · [生物力学](./生物力学/README.md) · [模式识别与机器学习](./模式识别与机器学习/README.md) · [生物统计学](./生物统计学/README.md) · [分子生物学](./分子生物学/README.md) |
| **大三下** | [生物材料](./生物材料/README.md) · [生物信息学](./生物信息学/README.md) · [医疗器械前沿](./医疗器械前沿/README.md) · [遗传学](./遗传学/README.md)  · [医学伦理学](./医学伦理学/README.md) |

## 📜 仓库原则与收录指南

我们不希望步三文鱼卷的后尘。所以一点收录规则还是有必要的：

### ✅ 我们欢迎并收集：
- 📝 **考试回忆**：历年真题回忆与期末重点梳理。（最有含金量的财富）
- 💡 **经验分享**：学习方法、备考经验与心得体会。
- 📚 **学习笔记**：优质的个人课堂笔记、总结与思维导图等。
- ⚠️ **关于格式**：推荐 Markdown 格式。当然，你也可以上传其他格式，交给笔者整理。

### 🚫 我们不接受：
- ⚠️ 涉及版权保护的官方教材、原版专著及电子书等。
- ⚠️ 未经授课教师允许外传的课件、PPT 、课程资料等。

欢迎大家通过提交 `Pull Request` 或 `Issues` 的方式参与共建，让知识在这里自由流动与共享！

## 🌐 本地构建网页

在仓库根目录运行：

```powershell
.\serve-site.ps1
```

然后访问 `http://127.0.0.1:8000/siwenyujuan/`。如需生成静态网站，运行 `.\build-site.ps1`，输出位于仓库旁的 `siwenyujuan-site` 目录。

## 💐 特别鸣谢

本项目在建设过程中得到了以下优秀项目的启发与参考：

- **百航全书**：[disk.stuhelper.com](https://disk.stuhelper.com/)

- **浙江大学课程攻略共享计划**：[QSCTech/zju-icicles](https://github.com/QSCTech/zju-icicles)
- **清华大学计算机系课程攻略**：[PKUanonym/REKCARC-TSC-UHT](https://github.com/PKUanonym/REKCARC-TSC-UHT)

网站基于 [MkDocs](https://www.mkdocs.org/) 构建。

## 待办
- [ ] 大三上课程资料及readme
