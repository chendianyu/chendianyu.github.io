---
title: Introduction to Human Reference Genome and Annotation
categories:
 - NGS
tags:
 - NGS
 - reference
 - annotation
---

当我们进行NGS数据分析时，通常我们需要将测序得到的reads比对至参考基因组上（reference genome）。目前有多家机构提供不同版本的参考基因组，针对我们
的需求，我们应当选择合适的参考基因组。此外，选择哪个参考基因组将会影响到后续对基因组注释的选择。接下来，让我们来了解一下不同的参考基因组及注释。

<!-- more -->

# Reference Genome
> A reference genome is a digital nucleic acid sequence database, assembled by scientists as a representative example of a species' 
set of genes. As they are often assembled from the sequencing of DNA from a number of donors, reference genomes do not accurately 
represent the set of genes of any single person. Instead a reference provides a haploid mosaic of different DNA sequences from 
each donor              \-- from wikipedia
  
目前存在多个版本的参考基因组，例如 `hg19`, `hg38`, 以及有多个不同机构维护其各自版本，包括 `USCS`, `NCBI`, `ENSEMBL`等。总体而言，不同机构维护
的版本是相似的，但是仍存在一定的差异。此外，还有一些项目会针对这些版本进行一些定制化或者改进，例如千人基因组计划等  
  
首先让我们看一下不同机构之间参考基因组的对应关系：  
| NCBI | UCSC | ENSEMBL |  
| --- | --- | --- |  
| GRCh37 | hg19 | ENSEMBL release_59/61/64/68/69/75 |  
| GRCh38 | hg38 | ENSEMBL release_76/77/78/80/81/82 |
  
