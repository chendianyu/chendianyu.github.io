---
title: pybedtools 使用
description: 
categories:
 - software
tags:
 - bedtools
 - tutorial
---

bed 文件是我们在生物信息学中常会碰到的一种格式，一般我们会通过 bedtools 等工具对这类文件进行处理。Pybedtools 对 bedtools 的各种功能进行了封装和拓展，从而实现在 python 中对基因组区间进行各种类型的操作。本文中我们将对 pybedtools 进行介绍。  

<!-- more -->

# Installation  
  
# BedTool Object
一般来说，单个 `BedTool` 指向单个区间文件，可以是 BED, GFF, GTF, VCF, SAM, or BAM format 以及对应的压缩格式等。`BedTool` 对象封装了所用可用的 bedtools 程序，使之可以在 python 中进行调用。
