---
title: SummarizedExperiment 介绍
description: SummarizedExperiment 对象是存储转录组测序等实验结果的结构，也是 DESeq 等包进行分析的基础，在此对其进行一个简单的介绍。
categories:
 - RNA-seq
tags:
 - RNA-seq
 - bioconductor
---

# Introduction  
`SummarizedExperiment` 类用于存储测序或者芯片试验的结果。 每个对象中包含一个或者多个样本的观察值，以及描述观察值（observations），即特征（features） 和样本（samples），即表型（phenotypes）的元数据（meta-data）。  
`SummarizedExperiment` 类的一个重要特征是**元数据与实验数据的协调**，例如想要剔除某个样本，那么单个操作即可完成对元数据和试验的处理，确保元数据和观察值保持同步。  
`SummarizedExperiment` 在其行信息更加灵活，可以是基于 `GRanges` 或者任意的 `DataFrame`。  
  
# 解析 SummarizedExperiment  
`SummarizedExperiment` 包中包含两个类：`SummarizedExperiment` 和 `RangedSummarizedExperiment`  
  
* SummarizedExperiment
  `SummarizedExperiment` 是一个矩阵样的容器，每一行代表感兴趣的**特征**（feature），如基因，转录本，外显子等，每一列则代表**样本**（sample）。每个对象可以包含一个或者多个**实验**（assay），而它们每一个都用类矩阵的对象表示。  
  特征相关的信息存储在一个数据框（DataFrame）对象中，可以通过函数 `rowData()` 来获取。数据框的每一行提供了 `SummarizedExperiment` 对象中对应每一行特征的信息，数据框中的每一列则提供了感兴趣的特征属性，例如基因或者转录本的 ID 等。  

* RangedSummarizedExperiment  
  `RangedSummarizedExperiment` 是 `SummarizedExperiment` 的子类，所以后者所有的方法在 `RangedSummarizedExperiment` 中也适用。  
  两者最主要的区别在于 `RangedSummarizedExperiment` 中的每一行均为 `genomic ranges`，而不是 包含特征的数据框. 这些 genomic ranges 通过 `GRanges` 或 `GRangesList` 对象描述，可通过 `rowRanges()` 函数来获取。  
  
下图是 `SummarizedExperiment` 对象的示意图，并对行列之间的关系进行了高亮。  
![SummarizedExperiment](/img/2018-09-07-SummarizedExperiment/SummarizedExperiment.svg)  
  
## 实验（Assay）  
可通过 `assays()` 访问器从 `SummarizedExperiment` 对象中获取实验数据；由于单个对象中可以包含多个实验数据集，因此可以通过 `$` 操作符获取具体某个数据集。  
```R  
> assays(se)
List of length 1
assays(1): counts  
> assays(se)$counts  
...  
```  
  
## 'Row' (regions-of-interest) data  
通过 `rowRanges()` 访问器获取 `RangedSummarizedExperiment` 中的 range information (`rowData()` for `SummarizedExperiment`)。数据存储在 `GRangesList` 对象中。  
  
## 'Column' (sample) data  
通过 `colData()` 访问器后去样本的元数据，将会返回数据框  
样本的元数据的信息可以通过 `$` 操作符来获取，因此我们可以基于此对实验数据取子集， e.g. `se[, se$dex == "trt"]`，注意**元数据和观察值将同步变化**  
  
## Experiment-wide metadata  
通过 `metadata()` 来获取描述整个实验设计以及参考文献等元数据信息  
`metadata()` 得到的结果是一个简单的列表，因此可以添加任意希望存档的数据，比如模型公式，e.g. `metadata(se)$formula <- counts ~ dex + albut`  
  
# SummarizedExperiment 常用操作  
## Subsetting  
* Two dimensional subsetting  
```R
se[1:5,1:3] 
se[1:5,] 
```  
* $ operates on `colData()` columns  
```R  
se[, se$cell == "N61311"]  
```  

## Getters and setters  
* `assays()` and `assay()`  
  `assays()` 将所有实验数据作为一个整体进行操作，而 `assay()` 每次仅作用于单个实验数据。 `assay(x, i)` 等价于 `assays(x)[[i]]` (`assay` defaults to the first assay if no `i` is given)  
* `rowRanges() / (rowData())`  
* `colData()`  
* `metadata()`  

# Construct DESeqDataSet object from RangedSummarizedExperiment object  
```R  
dds <- DESeqDataSet(se, design = ~ cell + dex)  
```  
  
# REF  
1. https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html
