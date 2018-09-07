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
`SummarizedExperiment` class is used to store rectangular matrices of experimental results. Each object stores observations of one or more samples, along with additional meta-data describing both the observations (features) and samples (phenotypes)
A key aspect of the `SummarizedExperiment` class is the **coordination of the meta-data and assays** when subsetting. For example, if you want to exclude a given sample you can do for both the meta-data and assay in one operation, which ensures the meta-data and observed data will remain in sync
`SummarizedExperiment` is more flexible in it’s row information, allowing both `GRanges` based as well as those described by `arbitrary DataFrames`. This makes it ideally suited to a variety of experiments, particularly sequencing based experiments such as RNA-Seq and ChIp-Seq  
  
# Anatomy of a SummarizedExperiment  
The SummarizedExperiment package contains two classes: `SummarizedExperiment` and `RangedSummarizedExperiment`  
  
* SummarizedExperiment
  `SummarizedExperiment` is a matrix-like container where rows represent `features` of interest (e.g. genes, transcripts, exons, etc.) and columns represent `samples`. The objects contain one or more `assays`, each represented by a matrix-like object of numeric or other mode  
  Information about features is stored in a `DataFrame` object, accessible using the function `rowData()`. Each row of the DataFrame provides information on the feature in the corresponding row of the SummarizedExperiment object. Columns of the DataFrame represent different attributes of the features of interest, e.g., gene or transcript IDs, etc.
* RangedSummarizedExperiment  
  `RangedSummarizedExperiment` is the child of the SummarizedExperiment class, so all the methods on SummarizedExperiment also work on a RangedSummarizedExperiment  
  The fundamental difference between the two classes is that the rows of a RangedSummarizedExperiment object represent `genomic ranges` of interest instead of a DataFrame of features. The RangedSummarizedExperiment ranges are described by a `GRanges` or a `GRangesList` object, accessible using the `rowRanges()` function  
  
## Assay  
To retrieve the experiment data from a SummarizedExperiment object one can use the `assays()` accessor. An object can have multiple assay datasets, each of which can be accessed using the `$` operator  
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
