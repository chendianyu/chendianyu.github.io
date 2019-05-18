---
title: 利用Bioconductor实现scRNA-seq数据的分析
description: Bioconductor提供了完成单细胞转录组分析所需的大量软件包，本综述向我们展示了如何利用Bioconductor实现整个分析流程
categories:
 - Single Cell
tags:
 - paper
 - single cell
 - review
---

# Introduction
Bioconductor 是一个基于 R 语言的生物学开源软件包仓库，吸引了来自不同领域，包括基因组学，蛋白质组学，代谢组学，转录组学等的大量开发者和用户。随着 scRNA-seq 技术的广泛应用，其产生的数据集的两大特性需要我们有效地解决：一是**样本量规模的大幅增加**，另一个则是**数据稀疏度的增加**。数据复杂度以及稀疏度的增加需要我们在数据获取，管理，分析等方面进行根本性地改变，并开发出新的方法来分析这些数据，而 Bioconductor 恰好提供了最新的分析软件和流程，从而实现：  
1. 高效的数据导入和表示  
2. 统一的数据容器，用于存储和软件包之间的交互  
3. 快速且强健的方法，将原始的单细胞数据转换成可用于后续分析的加工后数据  
4. 交互式数据可视化  
5. 下游分析，注释以及生物学解释  
  
本文后续将对整个分析流程以及用到的软件进行说明。此外，作者还提供了[在线电子书](https://osca.bioconductor.org)，囊括了更详细的说明和更多的相关资源  
  
# 测序数据预处理  
得到测序数据后，通常我们会将 reads 比对到参考基因组，量化，构建表达矩阵。R/Bioconductor 上提供了 `scPipe` 包，能够利用 `Rsubread` 包在 R 中实现比对。对于基于液滴的测序方法如 10X Genomics 等的数据，`DropletUtils` 包能够读入通过 `CellRanger` 流程产生的 UMI 计数矩阵  
此外，一些新的 `pseudoaligners` 使得比对能够在 PC 上完成，如 `Salmon` 和 `Kallisto`。`tximport` 包可以将这些比对软件的结果导入到 R 会话中，并在 `tximeta` 包的配合下，构建后续分析所需要的实例  
  
# 数据结构  
Bioconductor 分析流程的一大优势就是使用了统一的数据结构，从而实现模块化以及软件包之间的交互。Bioconductor 利用 S4 面向对象编程风格，对数据的存储和获取进行标准化，并在单个数据容器中纳入丰富的元信息注释，便于分析。下面我们将着重介绍 `SingleCellExperiment` 类  
`SingleCellExperiment` 类衍生自 `SummarizedExperiment` 类，因此继承了其适于大规模数据分析的优势，并增加了单细胞数据分析特需的方法和结构。不同信息部分称为 `slot`，以下是常见的一些 slot：  
* `assay` slot，对应最核心的数据，即表达矩阵。可以是 read 或 UMI 计数，且可以有一至多个矩阵，比如分别表示原始以及标准化后的计数。通常矩阵的行对应基因（或转录本等 feature），列对应细胞（sample）  
* `rowData` slot，行和列都会有丰富的注释信息，该 slot 中放置的即为行相关的元信息，比如基因的 ID，GC 含量等，特别地，如果行对应的是基因组特征，那么会有 `rowRanges` slot 存储基因组坐标信息  
* `colData` slot，列元信息，即样本层面的信息，可以自行增加列元信息，比如每个细胞中测到的 reads 数量等，适用于质控时的筛选  
* `sizeFactors` slot，指向列（细胞），包含标准化所需的信息  
* `reducedDims` slot，包含数据在低维度空间中的表示，可以有多种降维方法，如 PCA，t-SNE，UMAP  
  
另外，最近涌现了许多同时测量细胞中基因，表观遗传以及转录信息的方法，`MultiAssayExperiment` 包能够整合异质的数据类型，比如 `SingleCellExperiment`，`DelayedArray`等，以便进行后续分析  
  
# 质控和标准化
构建好 `SingleCellExperiment` 对象后，第一步需要做的就是鉴定、移除、校正低质量的 feature 和 sample，整个步骤包括：**过滤低质量的细胞**，**选择能提供有用信息的特征**，**细胞和基因层面的标准化**以移除细胞和基因特异性偏好性，以及**针对已知协变量和隐变量进行调整**  
  
## 细胞和基因质控 
对于基于液滴的测序方法，`DropletUtils` 包可用于鉴定空液滴以及减弱 *barcode swapping* 带来的影响。`scater` 包能够自动计算一系列质控指标，其中**文库大小**（library size，单个细胞中 read 或 UMI 总量）、**spike-ins 或线粒体基因上计数所占比例**、**检测到的基因数**等是常用的细胞筛选指标；而基因的**平均表达量**或**检测到的频率**可用于移除那些不提供信息的基因。`simpleSingleCell` 流程中展示了如何利用这些包进行质控  
  
## 基因表达标准化
细胞和基因水平的偏好可以通过**比例因子（scaling factor）**，也称**量化因子（size factor）**进行校正，从而使得不同的细胞能够进行比较。`SCnorm`，`scran` 和 `scater` 包均提供了标准化的方法，其中 `SCnorm` 包能够估算并消除由于测序深度，GC含量以及基因长度等引起的基因水平差异，`scran` 首先将表达谱相似的细胞聚集在一起，然后计算细胞特异性量化因子用于标准化，`scater` 包则是基因文库大小计算量化因子。需要注意的是，**以上方法只能计算能直接从数据中推导出来的内部因子，无法针对实验水平比如批次效应进行校正**  
