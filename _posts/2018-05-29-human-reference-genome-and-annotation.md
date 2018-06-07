---
title: Introduction to Human Reference Genome
categories:
 - NGS
tags:
 - NGS
 - reference
 - Genome
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
| GRCh37 | hg19 | ENSEMBL 59/61/64/68/69/75 |  
| GRCh38 | hg38 | ENSEMBL 76/77/78/80/81/82 |  

# 提供参考基因组的机构  
## The Genome Reference Consortium (GRC)
NCBI 旗下的 GRC 尽最大努力为我们提供最完善的人类参考基因组。除了基础的版本外，它还提供了 `alternate loci` 以表示那些过于复杂而无法用单条通路表示的区域。GRC 版本的基因组是 `one-based coordinate system`。
此外，它还以 `patch` 的形式修正某些区域，例如你会看到 `GRCh37.p7` 这样的版本号。通过这种方式，既可以保证整个染色体坐标的稳定性，同时提供更加准确的碱基组成。  
  
## UCSC  
UCSC 本身并不会进行基因组的测序及组装，它只是将 GRC 的版本进行了一定的调整，例如在染色体编号之前加上 `chr` 等。需要注意的是 UCSC 给出的基因组版本是 `zero-based coordinate system`。 其他具体的调整我们会在接下来各版本参考基因组中进行更详细的说明。  
  
## ENSEMBL
ENSEMBLE 同样是使用来自 GRC 给出的参考基因组。与 UCSC 不同的是，其使用的是 `one-based coordinate system`。

# Reference Genome
## GRCh37/hg19/b37/ENSEMBL 59/61/64/68/69/75
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
  * `unlocalized sequences`及`unplaced sequences`使用 `accession number` 命名，未包含 `alternate loci` 序列  
这些约定在之后的 `ENSEMBL genome browser`, `the NCBI dbSNP` (in VCF files), `the Sanger COSMIC` (in VCF files) 等中得以继承。后续多数新的项目更倾向于该标准  
* b37+decoy / hs37d5 (1000 Genomes Project Phase II)  
可以理解为b37的升级版，做了一些调整：  
  * A human herpesvirus (疱疹病毒) 4 type 1 sequence，命名为 `NC_007605`  
  * "decoy" sequence (诱饵序列)（命名为 `hs37d5`） 来自HuRef、BAC和质粒克隆和NA12878，可以提高序列比对的正确率  
  * Y染色体上 `pseudo-autosomal regions (PAR)` 被masked (用 `N`代替)，所以X染色体上对应区域可能被视为二倍体  
  以上这些变化使得该套序列能够减少假阳性，且与b37兼容，从而成为序列比对和变异识别的最佳之选  
  
## GRCh38/hg38/ENSEMBL 76/77/78/80/81/82
相较于 GRCh37，GRCh38 的碱基组成以及坐标位置都有所变化，此外添加了更多的 ALT contig，因此两者之间存在不小的差异。从2013年12月以来，38版本的基因组
及相关的注释其实都已经很完善了，因此会建议新的项目都是用该版本的参考基因  
  
# Which one to choose?  
BWA 作者 Heng Li 在其[博客](http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)中给出了他的建议：  
* 如果使用37版本的参考基因组，则选择 `hs37-1kg`：  
`ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz`  
* 如果使用37版本的参考基因组，且认为诱饵序列有助于之后的分析，则选择 `hs37d5`：  
`ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz`  
* 如果使用的是38版本的参考基因组，则选择：  
`ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`  
  
他列举了不同版本参考基因组中存在的问题，包括：  
1. 包含 `ALT` contig。这些序列中有很大部分侧翼序列与 `primary human assembly` 中的序列基本一致，大部分序列比对软件给到这些区域的 map quality 值为0，这就导致了后续变异识别以及其他分析的敏感性下降。可以通过选择 `ALT-aware` 的软件解决这个问题，但目前主流的软件并不支持  
2. `Padding ALT contigs`，当中有大量的 `N` 碱基，类似于1带来的问题  
3. 包含 `multi-placed` 的序列。例如 X 染色体上的 PAR 也放在 Y 染色体上。如果使用同时包含这两部分的参考基因组，就无法使用标准流程对这块区域进行变异识别。还有像在 GRCh38 中，部分 `alpha satellites` 被放置多次。正确的解决方法时间 Y 染色体上的 PAR 区域以及 alpha repeats 多出来的拷贝 hard maseked  
4. [rCRS](https://en.wikipedia.org/wiki/Cambridge_Reference_Sequence) 线粒体序列的使用。`rCRS` 在群体遗传性中大量使用，而 GRCh37 所使用的线粒体序列要比它长2bp，这可能导致在进行线粒体系统进化分析时出现问题。GRCh38 使用的则是 rCRS  
5. 将半模糊的 [IUB codes](http://biocorp.ca/IUB.php) 转换成了 `N`，虽然总体而言影响不大，因为人类染色体序列中只有极少部分这样的碱基  
6. 用 `accession number` 代替染色体名字进行命名，不够直观  
7. 没有包含 `unplaced and unlocalized contigs`  
  
具体到个版本参考基因组，问题就是这样的：  
* hg19/chromFa.tar.gz from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/): 1, 3, 4 and 5  
* hg38/hg38.fa.gz from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/): 1, 3 and 5  
* GCA_000001405.15_GRCh38_genomic.fna.gz from [NCBI](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/): 1, 3, 5 and 6  
* Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz from [EnsEMBL](http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/): 3  
* Homo_sapiens.GRCh38.dna.toplevel.fa.gz from [EnsEMBL](http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/): 1, 2 and 3.
  
**当然，怎么选择参考基因组还是需要根据具体的情况来判定，包括实验室或者公司之前使用的情况，不同版本对应的注释文件情况等等**  

# Ref
1. https://genestack.com/blog/2016/07/12/choosing-a-reference-genome/  
2. http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use  
3. https://wiki.dnanexus.com/Scientific-Notes/human-genome  

