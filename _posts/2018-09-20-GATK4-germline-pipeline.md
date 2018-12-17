---
title: GATK4 WES/WGS 分析流程
description: GATK4 检测生殖系突变流程
categories:
 - GATK
tags:
 - GATK4
 - software
 - pipeline
 - WES
 - WGS
---

# Introduction  
本流程基于 GATK4 进行 WES/WGS 生殖系突变，主要是 SNP 和 Indel 的检测。流程使用 `GRCh38` 作为参考基因组，以双端测序 FastQ 文件起始，最终得到包含 SNP 和 Indel 的 (g)VCF 文件  

# Reference  
GRCh38 参考基因组 FASTA 文件包含 `alternate contigs`, 此外还需要用到 `.fai` index 文件， `.dict` dictionary 文件以及5个 BWA-specific 索引文件 `.sa`, `.amb`, `.bwt`, `.ann` and `.pac`    
1. a `.dict` dictionary of the contig names and sizes  
我们通过 `Picard` 的 `CreateSequenceDictionary`，基于参考基因组的 FASTA 文件构造出一个后缀名为 `.dict` 的序列字典文件。该文件实质为一个 SAM 格式的文件，**包含 header，但无 SAMRecords，且 header 部分仅包含序列记录，包括序列的名称和长度等信息**。参考基因组文件可以被压缩，即支持 `.fasta.gz` 格式    
Usage:  
```shell
java -jar picard.jar CreateSequenceDictionary \ 
      R=<ref.fasta> \ 
      O=<ref.dict>  
```  
e.g.:  
```shell  
java -jar picard.jar CreateSequenceDictionary \ 
      R=/home/refer/ucsc.hg38.fasta \ 
      O=/home/refer/ucsc.hg38.dict  
```  

2. a `.fai` fasta index file to allow efficient random access to the reference bases  
我们通过 `Samtools` 的 `faidx` 构建参考基因组的索引文件 `<ref.fasta>.fai`。输入文件可以是 BGZF 格式的压缩文件。根据 `.fai` 文件和原始的`fasta` 文件，我们能够快速的提取任意区域的序列  
Usage : `samtools faidx <ref.fasta>`  
e.g. : `samtools faidx /home/refer/ucsc.hg38.fasta`  
  
3. BWA 比对也需要其对应的索引文件，索引是**算法特异性**。我们通过 `BWA` 的 `index` 构建索引，得到后缀名为 `.amb`, `.ann`, `.bwt`, `.pac` and `.sa` 的文件，之后工具将会默认这些索引文件与参考基因组位于**同一文件夹内**  
Usage : `bwa index -a bwtsw <ref.fasta>`  
e.g. : `bwa index -a bwtsw /home/refer/ucsc.hg38.fasta`    

当然也可以进入 GATK resource bundle 下载已准备好的文件  
```shell
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi &
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz &
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_omni2.5.hg38.vcf.gz.tbi
```  

# Expected input
数据最初是按照不同的 `readgroups` 分成不同的子集，对应 `libraries`（DNA 取自不同的生物学样本，以及在文库制备过程中的片段化和 barcode 标记过程）与  `lanes` (测序仪的物理分隔单元) 通过 `multiplexing`（混合多个文库，并在多条泳道上进行测序，以减少风险和人为误差）交叉得到的不同组合  
Read groups 在 SAM/BAM/CRAM 中通过一系列标签进行定义。在 SAM 等文件的 header 中，以 `@RG` 起始，后面跟着以下几个 tag：  
* ID：每个 read group 的 ID 是唯一的。对于 Illumina 数据，一般就是 flowcell + lane name and number。对于 BQSR 等处理，每一个 read group 被认为是独立的，被认为共享相同的错误模型  
* PU：Platform Unit，包含 {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} 信息。{FLOWCELL_BARCODE} 是 flowcell 的唯一标识，{LANE} 表示 flow cell 中的 lane，而 {SAMPLE_BARCODE} 则是 sample/library-specific identifier。PU 不是 GATK 必需的，但如果出现，将优先于 ID 作为 BQSR 时分组的依据  
* SM：该 read group 中测序的样本的名称。GATK 将所有具有相同 SM 的数据视为来源于同一样本，这也将会是 VCF 文件中 sample column 所使用的名称  
* PL：平台或使用的技术名称。Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO  
* LB：文库标识符。当同一 DNA 文库在多条 lane 进行测序时，`MarkDuplicates` 会基于 LB 来确定哪些 read group 可能含有分子重复  

正常情况下我们会在 BWA 比对时候加上 read group 信息，如下所示。但假如忘了，可以通过 Picard `AddOrReplaceReadGroups` 补上：  
```shell
java -jar picard.jar AddOrReplaceReadGroups \
    I= <reads_without_RG.bam> \
    O= <reads_with_RG.bam> \
    SORT_ORDER=coordinate \
    RGID=<group1> \
    RGLB=<lib1> \
    RGPL=<illumina> \
    RGSM=<sample1> \
    CREATE_INDEX=True
```

# Map to Reference  
**Tools involved:** `BWA`, Picard 的`MergeBamAlignments`
第一步按照每个 read group，将每对 read 比对到参考基因组上。由于比对算法能够独立处理每对 read，因此可以通过并行提高速度  

```shell
bwa mem -R <read_group> \  # e.g. : '@RG\tID:group1\tSM:sample1\tPL:illumina\tLB:lib1\tPU:unit1'
    -M \ 
    -t <threads_num> \ 
    <db.prefix> \           # 存储ref.fa及相关index等的目录
    <Read1> <Read2> > \
    <input.sam>
```

* `-M` : Mark shorter split hits as secondary (for Picard compatibility)
  
# Coordinate sort and index  
Tools involved: Picard's `SortSam` and `BuildBamIndex`  
后续变异识别和可视化比对情况需要 SAM/BAM 文件先进行排序。`SortSam` 能够根据 **coordinate**, **queryname (QNAME)**, or **some other property** 对 SAM/BAM 文件排序，排序依据存放在 `@HD` tag 的 `SO` 字段中  
按照坐标排序时，read 首先参考序列字典（`@SQ` tag）中的参考序列名称（`RNAME` 字段）排序，然后是最左边的比对位置（`POS`），在之后就随机排序  
按照 query name（`QNAME`）排序时，按照 query name 聚合，但各组内的序列并不一定进行了排序。具有相同 query name 的 reads 来源于同一模版  
排序后得到的文件为 `reads_sorted.bam`，使用 `CREATE_INDEX` 参数能同时对应的索引文件，后缀为 `.bai`，e.g. `reads_sorted.bai`。 索引文件用于快速检索 BAM 文件。需要注意该参数**不能用于 SAM 文件**，且 **BAM 文件必须按照坐标轴排序**  
Usage :  
```shell
java -jar picard.jar SortSam \
    INPUT=<input.sam> \
    OUTPUT=<output_sorted.bam> \ 
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
```  
  
可以用 `BuildBamIndex` 对排完序的 BAM 文件构建索引  
```shell  
java -jar picard.jar BuildBamIndex \
      I=<input.bam>  
```  
  
或者 `samtools index`  
Usage: `samtools index <input.bam>` 产生的文件为 `<input.bam.bai>` 只有这个与 Picard 有区别，文件内容本质上应该是一致的  

## Mark Duplicates
Tools involved: `MarkDuplicates`  
重复可以是在样本准备过程中发生，如通过 PCR 构建文库，称为 `PCR duplicates`；也可以是单个扩增簇被测序仪的光学传感系统误认为是多个簇导致，称为 `optical duplicates`  
重复标记过程是按照每个样本（`per-sample`）进行的，标记重复序列，从而在后续变异识别过程中忽略这些 reads。该步骤需要对每个 sample 内所有 reads 进行两两比较，因此是主要的限速步骤  
输出为 **a new SAM or BAM file**, 重复 reads bitwise flag 标记为十六进制值 `0x0400`, 对应十进制值为1024；另外为了标记重复的类型，最近在 SAM/BAM 文件的 'optional field' section 引入了一个新的 tag。通过 `TAGGING_POLICY` 选项，可以控制程序标记所有重复（All），仅光学重复（OpticalOnly）或者不标记重复（DontTag，**默认**）。输出 SAM/BAM 文件中会带上 `DT` tag，值为 `library/PCR-generated duplicates (LB)`, or `sequencing-platform artifact duplicates (SQ)`  
此外还会输出一个 metrics file，标记 reads 的重复数；使用 `CREATE_INDEX` 能够得到索引文件  
  
Usage : 
```shell
java -jar picard.jar MarkDuplicates \
    INPUT=<input_sorted.bam> \
    OUTPUT=<output_sorted_duplicates.bam> \ 
    METRICS_FILE=<output_markduplicates_metrics.txt> \ 
    CREATE_INDEX=true
```  
  
可以通过 `REMOVE_DUPLICATE`（剔除所有重复序列） 和 `REMOVE_SEQUENCING_DUPLICATES`（剔除所有因测序而不是样本制备导致的重复）参数来剔除重复序列   
(另外可以使用 Picard's `FixMateInformation` 确认成对 reads 之间的信息是否一致；如有必要，进行修复)  
Usage:  
```shell
java -jar picard.jar FixMateInformation \ 
    I=<input_sorted_duplicates.bam> \ 
    O=<output_sorted_duplicates_fixed.bam> \ 
    ADD_MATE_CIGAR=true  
```
  
# Base (Quality Score) Recalibration
Tools involved: `BaseRecalibrator`, `Apply Recalibration`, `AnalyzeCovariates` (optional)
变异识别算法主要依赖碱基质量分数，但是机器得到的分值受各种来源的误差干扰，比如文库制备的生化过程，测序芯片或者测序仪等的缺陷，并不准确。`Base quality score recalibration (BQSR)` 是基于机器学习方法，对这些错误进行建模，然后校正质量分数。这一步是按照每个 read group 进行的
BQSR 包括两个主要的步骤：  
* `BaseRecalibrator` 基于输入数据和一组已知变异（一般是 dbSNP）构建协变模型，生成重校准文件。已知变异的作用是**将实际（预期）发生变异位点的碱基 mask，避免将真实的变异视为错误**，除了这些位点，所有不匹配的位置均视为错误  
* `ApplyBQSR` 基于重校准文件对质量分数进行校正，生成一个新的 BAM 文件    

为了构建重校准模型，BaseRecalibrator 会遍历所有 reads，按照碱基的下列特征制作表格：  
* read 所属的 read group  
* 碱基质量分数  
* 得到该碱基时的机器循环次数，即碱基在 read 中所处的位置  
* 当前位置的碱基和前一个位置的碱基（dinucleotide，考虑测序时化学试剂的影响）  
  
对于每个 bin 中的碱基，排除已知变异稳点后，计算碱基数量和与参考序列不匹配的碱基数量，此信息将以 GATKReport 格式输出到重校准文件  
然后 ApplyBQSR 根据每个碱基所处的 bin 对其进行校准，新的质量分数值为：  
* 报告的碱基质量分数和实际值之间的总体差异  
* 加上每个 quality bin 特异的偏移值  
* 加上循环和 dinucleotide 的影响值  
  

## Steps  
1. Analyze patterns of covariation in the sequence dataset  
```shell
gatk BaseRecalibrator \ 
    -I <input_sorted_markduplicates.bam> \
    -R <ref.fasta> \ 
    -knownSites <dbsnp.vcf> \ 
    -knownSites <1000g.vcf> \ 
    -O <recal_data.table>  
```

2. Apply the recalibration to your sequence data  
```shell
gatk ApplyBQSR \  
    -R <ref.fasta> \ 
    -I <input_sorted_markduplicates.bam> \ 
    --bqsr-recal-file <recal_data.table> \ 
    -O <output_sorted_markduplicates_recal.bam> 
```
* By default, the original quality scores are discarded in order to keep the file size down
  
# Variants calling  

# REF
1. https://gatkforums.broadinstitute.org/gatk/discussion/11165/data-pre-processing-for-variant-discovery
