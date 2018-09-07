---
title: SummarizedExperiment 介绍
description: SummarizedExperiment 类的介绍
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
  `SummarizedExperiment` 是一个矩阵样的容器，每一行代表感兴趣的**特征**（feature），如基因，转录本，外显子等，每一列则代表**样本**（sample）。每个对象可以包含一个或者多个**实验**（assay），而它们每一个则用类矩阵的对象表示。  
  特征相关的信息存储在一个数据框（DataFrame）对象中，可以通过函数 `rowData()` 来获取。数据框的每一行提供了 `SummarizedExperiment` 对象中对应每一行特征的信息，数据框中的每一列则提供了感兴趣的特征属性，例如基因或者转录本的 ID 等。  

* RangedSummarizedExperiment  
  `RangedSummarizedExperiment` 是 `SummarizedExperiment` 的子类，所以后者所有的方法在 `RangedSummarizedExperiment` 中也适用。  
  两者最主要的区别在于 `RangedSummarizedExperiment` 中的每一行均为 `genomic ranges`，而不是 包含特征的数据框. 这些 genomic ranges 通过 `GRanges` 或 `GRangesList` 对象描述，可通过 `rowRanges()` 函数来获取。  
  
下图是 `SummarizedExperiment` 对象的示意图，并对行列之间的关系进行了高亮。  
![SummarizedExperiment](/img/2018-09-07-SummarizedExperiment/SummarizedExperiment.svg)  
  
## 实验（Assay）  
可通过函数To retrieve the experiment data from a SummarizedExperiment object one can use the `assays()` accessor. An object can have multiple assay datasets, each of which can be accessed using the `$` operator  
```R  
> assays(se)
List of length 1
names(1): counts  
> assays(se)$counts  
...  
```  
  
## 'Row' (regions-of-interest) data  
Using `rowRanges()` accessor to view the range information for a `RangedSummarizedExperiment`. (`rowData()` for `SummarizedExperiment`)  
  
## 'Column' (sample) data  
Using `colData()` accessor to view sample meta-data describing the samples, and it returns a DataFrame that can store any number of descriptive columns for each sample row  
The sample metadata can be accessed using the `$` accessor, e.g. `se$dex`. So it's easy to subset the entire object by a given phenotype  
  
## Experiment-wide metadata  
Using `metadata()` to get Meta-data describing the experimental methods and publication references  
`metadata()` is just a simple `list`, so user can add any metadata they wish to save, such as storing model formulas  
  
# Common operations on SummarizedExperiment  
## Subsetting  
* Two dimensional subsetting  
```R
se[1:5,1:3] 
se[1:5,] 
```  
* $ operates on colData() columns  
```R  
se[, se$cell == "N61311"]  
```  

## Getters and setters  
There are two accessor functions for extracting the assay data: `assays()` and `assay()`  
`assays()` operates on the entire list of assay data as a whole, while `assay()` operates on only one assay at a time. `assay(x, i)` is equivalent to `assays(x)[[i]]` (`assay` defaults to the first assay if no `i` is given)  
  
# Construct DESeqDataSet object from RangedSummarizedExperiment object  
```R  
dds <- DESeqDataSet(se, design = ~ cell + dex)  
```  
