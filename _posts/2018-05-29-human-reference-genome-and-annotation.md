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

> A reference genome is a digital nucleic acid sequence database, assembled by scientists as a representative example of a species' 
set of genes. As they are often assembled from the sequencing of DNA from a number of donors, reference genomes do not accurately 
represent the set of genes of any single person. Instead a reference provides a haploid mosaic of different DNA sequences from 
each donor              \-- from wikipedia
  
目前存在多个版本的参考基因组，例如 `hg19`, `hg38`, 以及有多个不同机构维护其各自版本，包括 `USCS`, `NCBI`, `ENSEMBL`等。总体而言，不同机构维护
的版本基本是相同的，但是仍存在一定的差异。此外，还有一些项目会针对这些版本进行一些定制化或者改进，例如千人基因组计划等  
  
以下是不同机构之间参考基因组的对应关系：  

| NCBI | UCSC | ENSEMBL |  
| --- | --- | --- |  
| GRCh37 | hg19 | ENSEMBL release_59/61/64/68/69/75 |  
| GRCh38 | hg38 | ENSEMBL release_76/77/78/80/81/82 |  

# 提供参考基因组的机构  
## The Genome Reference Consortium (GRC)
NCBI 旗下的 GRC 尽最大努力为我们提供最完善的人类参考基因组。除了基础的版本外，它还提供了 `alternate loci` 以表示那些过于复杂而无法用单条通路表示的区域。
此外，它还以 `patch` 的形式修正某些区域，例如你会看到 `GRCh37.p7` 这样的版本号。通过这种方式，既可以保证整个染色体坐标的稳定性，同时提供更加准确的碱基组成。  
  
## UCSC  
UCSC本身并不会进行基因组的测序及组装，它只是将GRC的版本进行了一定的调整，例如在染色体编号之前加上 `chr` 等，具体的我们会在接下来各版本参考基因组中进行更详细的说明。
  
## ENSEMBL  

# Reference Genome
## hg19/GRCh37/b37/ENSEMBL release_59/61/64/68/69/75
* GRCh37 中包含  
  * 24条基本完整的染色体序列，即 1-22, X, Y  
  * 完整的线粒体序列  
  * `unlocalized sequences`, 指知道具体来自哪条染色体，但具体坐标未知的序列  
  * `unplaced sequences`, 指无法确定来自于那条染色体的序列  
  * `alternate loci`, 这些序列包含人类基因组特定区域的 alternate representations
由于GRC给出这些序列时并没有给出标准的命名，也没有指定顺序，导致其他不同的组织采用各自的方式进行标记  
* hg19  
  * 针对24条基本完整的染色体序列，使用 `chr1`-`chr22`, `chrX` and `chrY` 的命名  
  * 没有使用来自GRCh37的线粒体 `NC_012920`, 而是使用了老版本的 `NC_001807`，命名为 `chrM`  
  * 针对 `unlocalized sequences`, `unplaced sequences`, `alternate loci` 序列，给出了自定义的名称  
  * 用小写字母标记 `repeats` and `low complexity regions`, 即进行 `soft-masked`  
尽管hg19版本存在不少的问题，但通过 `UCSC genome browser` 带来的大量曝光，hg19事实上应用广泛，包括不少外显子捕获试剂制造商利用该版本进行坐标定位  
* b37 (1000 Genomes Project Phase I)  
  * 24条基本完整的染色体序列，命名为 `1`-`22`, `X`, `Y`  
  * 使用来自GRCh37的线粒体，命名为 `MT`  
  * `unlocalized sequences`及`unplaced sequences`使用accession number命名，未包含 `alternate loci` 序列  
这些约定在之后的 `ENSEMBL genome browser`, `the NCBI dbSNP` (in VCF files), `the Sanger COSMIC` (in VCF files) 等中得以继承。后续多数新的项目更倾向于该标准  
* b37+decoy / hs37d5 (1000 Genomes Project Phase II)  
可以理解为b37的升级版，做了一些调整：  
  * A human herpesvirus (疱疹病毒) 4 type 1 sequence，命名为 `NC_007605`  
  * "decoy" sequence (诱饵序列)（命名为 `hs37d5`） 来自HuRef、BAC和质粒克隆和NA12878，可以提高序列比对的正确率  
  * Y染色体上 `pseudo-autosomal regions (PAR)` 被masked (用 `N`代替)，所以X染色体上对应区域可能被视为二倍体  
以上这些变化该套序列能够减少假阳性，且与b37兼容，从而成为序列比对和变异识别的最佳之选
  
# Ref
1. https://genestack.com/blog/2016/07/12/choosing-a-reference-genome/  
2. http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use  
3. https://wiki.dnanexus.com/Scientific-Notes/human-genome  

