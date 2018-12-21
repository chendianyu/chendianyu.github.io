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
本流程基于 GATK4 进行 WES/WGS 生殖系突变，主要是 SNP 和 Indel 的检测。流程使用 `GRCh38` 作为参考基因组，以双端测序 FastQ 文件起始，最终得到包含 SNPs 和 Indels 的 (g)VCF 文件  

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
Tools involved: Picard's `MarkDuplicates`  
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
Tools involved: GATK `BaseRecalibrator`, `Apply Recalibration`  
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
Tools involved: GATK `HaplotypeCaller`  
HaplotypeCaller 通过对活跃区域局部重组装实现变异识别：  
1. 定义活跃区域。首先是计算基因组各位置的活跃值，得到 raw activity profile，然后通过平滑算法得到 actual activity profile，最后根据活跃度曲线找出局部极大值，并以此确定出准确的区间  
2. 对于每个活跃区间，通过 De Bruijn-like graph 进行从头组装，确定可能的单体型，然后通过 Smith-Waterman 算法将单体型与参考单体型比对，确定可能存在变异的位点  
3. 通过 PairHMM 算法将活跃区内的 reads 成对比对到单体型上，会得到单体型似然度矩阵，将似然度边缘化,得到给定 reads 情况下可能存在变异位点各等位基因的似然度  
4. 对每个可能的变异位点，根据贝叶斯法则和各等位基因的似然度确定基因型  
  
## 单个样本  
Usage：  
```shell
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R <ucsc.hg38.fasta> \
   -I <input_sorted_markduplicates_recal.bam> \
   -O <output.vcf.gz> 
```  
  
## 多样本 GVCF  
GVCF 代表 genomic VCF，相较常规的 VCF 文件，包含更多的信息，适用于多个样本的变异识别。GVCF 包含基因组（或指定区间）**所有的位点**，而不管样本在这个位点是不是存在变异，用于后续 joint analysis  
GVCF 有两种模式，分别通过 `-ERC GVCF` 和 `-ERC BP_RESOLUTION` 得到。其中前者会将连续且基因型质量分数（genotype quality，GQ）在一定区间内的非变异位点整合成 block，并在 header 部分标注 `#GVCFBlock` 的信息；后者则是每个位点占据一行。GVCF 模式文件较小，一般采用该种模式  
Usage：  
```shell
# 单样本生成中间 GVCF 文件
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R <ucsc.hg38.fasta> \
   -I <input_sorted_markduplicates_recal.bam> \
   -ERC GVCF \
   -O <output.g.vcf.gz>
# GVCF 文件联合变异识别只接受单个输入，所以想得把样本中间结果合并
gatk CombineGVCFs \
   -R <ucsc.hg38.fasta> \
   --variant <sample1.g.vcf.gz> \
   --variant <sample2.g.vcf.gz> \
   -O <cohort.g.vcf.gz>
# 联合变异识别
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R <ucsc.hg38.fasta> \
   -V <cohort.g.vcf.gz> \
   -O <output.vcf.gz>
```
  
# Select Variants  
不同类型的变异具有不同的特征，因此对应的过滤条件也有区别，所以我们需要先将变异按照 SNP 和 Indel 区分  
Usage：  
```shell
# select SNPs
gatk SelectVariants \
     -R <ucsc.hg38.fasta> \
     -V <input.vcf.gz> \
     --select-type-to-include SNP \
     -O <raw_snps.vcf>
# select Indel
gatk SelectVariants \
     -R <ucsc.hg38.fasta> \
     -V <input.vcf.gz> \
     --select-type-to-include INDEL \
     -O <raw_indels.vcf>
```
  
# Variant Filter  
## hard filter
基于 INFO (`--filter-expression` 参数) 和 FORMAT (`--genotype-filter-expression` 参数) 注释对 VCF 文件进行过滤，达到标准的变异 FILTER 列显示为 PASS，未达到的则是其他值，过滤掉的变异默认也会输出  
Usage:  
```shell
# SNP 过滤
gatk VariantFiltration \
   -R <ucsc.hg38.fasta> \
   -V <raw_snps.vcf> \
   -O <filtered_snps.vcf> \
   --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
   --filter-name "my_snp_filters"    
# 满足表达式的变异将被过滤
# filter expressions 和 filter names 之间必须一一对应
# 列出多个过滤表达式和对应的名称
# --filter-name One --filter-expression "X < 1" --filter-name Two --filter-expression "X > 2"
#
# Indel 过滤
gatk VariantFiltration \
   -R <ucsc.hg38.fasta> \
   -V <raw_indels.vcf> \
   -O <filtered_indels.vcf> \
   --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
   --filter-name "my_indel_filters"
```
## VQSR (Variant Quality Score Recalibration)
基于机器学习方法对训练集中变异进行建模，并将其用于实际变异集的过滤。VQSR 不会对原来的质量分数值进行校正，而是将 QUAL 列中的分值没有考虑的各项与变异相关的属性纳入其中，重新计算一个新的质量分数 VQSLOD (for variant quality score log-odds)，尽可能平衡好敏感性和准确性。VQSLOD 分值，是真实变异与假阳性之间的 log odds ratio，添加至 INFO 列中  
VQSR 适用于单个全基因组，或者30个全外显子组，少于这个数据量可能导致效果不佳，尤其是 indel 重校准
`VariantRecalibrator` 找出训练/真实集与实际变异集之间的重叠部分，然后对变异与所指定的注释项之间的分布建模，期望将它们聚集成簇，然后根据这些簇对变异赋上 VQSLOD 值，越接近各分支中心的变异的分值要高于位于边缘的变异。这一步会得到 VCF 格式的重校准文件和一些图表附件  
`ApplyRecalibration` 会进行过滤，在输出的 VCF 文件中标记哪些变异通过了筛选条件，哪些则没有。在前一步中，变异获得了 VQSLOD 值，同时真实集中的变异会按照分支的大小排序。当通过 `ApplyRecalibration` 指定 tranche sensitivity （以百分比表示，如99.9%）后，程序将会查找高于那个分值时，真实集中99.9%的变异都会囊括在内，然后用这个值作为阈值筛选自己数据中的变异，超过这个值的变异 FILTER 标记为 PASS，没超过的则标记为如 VQSRTrancheSNP99.90to100.00，意味着其 VQSLOD 值处在真实集后0.1%的范围内  
Usage：
```shell
# SNPs VQSR
gatk VariantRecalibrator \
   -R <ucsc.hg3838.fasta> \
   -V <raw_snps.vcf> \
   --resource hapmap,known=false,training=true,truth=true,prior=15.0:hapmap_3.3.hg38.sites.vcf.gz \
   --resource omni,known=false,training=true,truth=false,prior=12.0:1000G_omni2.5.hg38.sites.vcf.gz \
   --resource 1000G,known=false,training=true,truth=false,prior=10.0:1000G_phase1.snps.high_confidence.hg38.vcf.gz \
   --resource dbsnp,known=true,training=false,truth=false,prior=2.0:Homo_sapiens_assembly38.dbsnp138.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O <snps_recalibrate.recal> \
   --tranches-file <snps_recalibrate.tranches> \ 
   --rscript-file <snps_recalibrate.plots.R>
gatk ApplyVQSR \
   -R <ucsc.hg38.fasta> \
   -V <raw_snps.vcf> \
   -O <snps_recalibrate.vcf> \
   --truth-sensitivity-filter-level 99.5 \
   --tranches-file <snps_recalibrate.tranches> \
   --recal-file <snps_recalibrate.recal> \
   -mode SNP
##
# Indels VQSR
gatk VariantRecalibrator \
   -R <ucsc.hg38.fasta> \
   -V <raw_indels.vcf> \
   --maxGaussians 4 \
   -resource:mills,known=false,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf  \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf\
   -an QD -an DP -an FS -an SOR -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff \
   -mode INDEL \
   -O <indels_recalibrate.recal> \
   --tranches-file <indels_recalibrate.tranches> \
   --rscript-file <indels_recalibrate.plots.R>
gatk ApplyVQSR \
   -R <ucsc.hg38.fasta> \
   -V <raw_indels.vcf> \
   -O <indels_recalibrate.vcf> \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file <indels_recalibrate.tranches> \
   --recal-file <indels_recalibrate.recal> \
   -mode INDEL
```
  
# Annotation  
Tools involved: `ANNOVAR`  
ANNOVAR 能够利用公共数据库中最新的信息对基因组变异进行注释。当给定染色体，起始位置，终止位置，参考碱基和替换碱基后，ANNOVAR 能够：  
* Gene-based annotation，确定 SNPs 或 CNVs 是否造成蛋白质编码的改变，以及所影响的氨基酸。用户可以选择 RefSeq，UCSC，ENSEMBL，GENCODE，AceView 或者其他的基因定义系统  
* Region-based annotation，确定变异是否存在于特定基因组区域，如不同物种之间的保守区域，预测的转录因子结合位点等  
* Filter-based annotation，确认变异是否在 dbSNP 数据库中，在千人基因组中的频率是多少，变异功能预测分值等  
* Other functionalities，其他自定义的一些数据库  

Usage：  
```shell
perl table_annovar.pl \ 
    <clean_variants.vcf> \
    <humandb/> \ 
    -buildver <hg38> \ 
    -outfile <variants_annovar> \ 
    -remove \     # 移除所有中间文件
    -protocol <refGene,cytoBand,exac03,avsnp147,dbnsfp30a> \ 
    -operation <g,r,f,f,f> \ 
    -nastring <.> \ 
    -vcfinput
```  
* `-operation`：告知 ANNOVAR 对各 protocol 进行何种操作  
    * `g`：gene-based  
    * `gx`：gene-based with cross-reference annotation (需要通过 `-xref` 参数指定 xref 文件)  
    * `r`：region-based  
    * `f`：filter-based  
* `--outfile`：输出文件名的前缀  
* `--otherinfo`：invoke `--includeinfo` in `convert2annovar.pl` (default `-includeinfo -allsample -withfreq`)  
* `--argument <string>`：各操作的可选参数，用逗号分隔。e.g. : `-arg '-splicing 5',,,`  
* 如果 ANNOVAR 遇到不符合标准的输入行，就会将改行输出至文件 `<outfile>.invalid_input`，`<outfile>` 是 `--outfile` 参数指定。如果都符合标准，该文件不存在  
* 上述程序将会输出文件 `<variants_annovar>.<hg38>_multianno.vcf`，`<variants_annovar>.<hg38>_multianno.txt` (tab-delimited) 和 `<variants_annovar>.avinput`  
  
# Others
## VQSR 使用的训练资源和注释
### For SNPs
* Truth resource：认证过的高可信度变异集，认为当中的变异都是真实的（true=true），可用于训练重校准模型（training=ture），之后还会将其用于根据敏感度确定阈值，如 Hapmap 数据库 hapmap_3.3.hg38.sites.vcf.gz 和 Omni 芯片数据库 1000G_omni2.5.hg38.sites.vcf.gz  
* Training resource：具有一定可信度的变异集，认为当中既有真实变异，也存在假阳性（true=false），会将其用于训练（training=ture），如千人基因组数据库 1000G_phase1.snps.high_confidence.hg38.vcf.gz  
* Known sites resource：未经过验证的变异集（true=false），不会用于训练（training=false），但会根据变异是否存在于其中作为 Ti/Tv 比率等输出指标的分层依据，如 dbSNP 数据库 dbsnp_146.vcf.gz  
  
### For Indels
* Truth resource：Mills 的 Mills_and_1000G_gold_standard.indels.hg38.vcf  
* Known sites resource：dbSNP 的 dbsnp_146.hg38.vcf  

### 使用的注释  
* 外显子数据不宜使用测序深度 `DP`，因为捕获的靶序列在测序深度上存在极大的差异，not indicative for error  
* `InbreedingCoeff` 是群体水平的统计量，至少需要10个样本，因此不适用于样本数较少或者样本之间关系较近（如家系）的情况  
  
## Genotype likelihoods and genotype quality  
VCF 文件中存在多种质量分数：  
* `QUAL` field  
计算公式：$ -10log_{10} prob(call in ALT is wrong) $ 。如果 ALT 为 `.` (即没有变异)，那么该值为 $ -10log_{10} prob(variant) $；如果 ALT 不为 `.`，那么该值为 $ -10log_{10} prob(no variant) $；如果未知，则给出缺失值  
QUAL 告诉我们的是给定位点存在某种变异的可靠性，该变异可能出现在一个或者多个样本中。它有一个标准化后的值 `QualByDepth (QD)`，有变异可信度（QUAL）除以非参考样本的未过滤深度（DP）得到，以避免深度测序带来的膨胀  
QUAL（或者 QD）在多样本情境中更有用。当对一组数据进行重校准时，看的是位点水平的注释，因为希望是了解变异的整体情况  
* `PL` in `FORMAT` field  
“标准化”后的每种基因型的 Phred-scaled 似然度，计算公式为 $ PL = −10 log_{10} P(Genotype|Data) $，得到每种基因型的 PL 之后，都减去当中最小的那个值，那么最有可能的那种基因型对应的 PL 值为0  
* `GQ` in `FORMAT` field  
给出的基因型（GT）的可靠程度。GQ 等于第二可能的基因型的 PL 值减去最可能基因型的 PL 值，所以即为第二小的 PL 值。**GATK 中该值上限设为99**  
  
## 深度 
* `AD (DepthPerAlleleBySample)` 输出的是每个样本中每个等位基因的未经过滤的测序深度  
* `DP (Coverage)` 输出的则是每个样本经过过滤的测序深度和所有样本未经过滤的测序深度  

## Allele-specific filtering （测试中版本） 
传统的 VQSR 重校正针对的是每个位点，对于一些含有多个等位基因的位点可能会导致假阴性。Allele-Specific filtering 会将每个位点的等位基因独立对待，适用于多等位基因的位点，一般在大样本中表现较好  
与传统的 VQSR 流程相比，不需要其他的资源，不过需要从样本的 BAM 文件开始操作。当运行完 Allele-Specific filtering 后，会在 INFO 一列添加一些新的注释，之后 VQSR 会针对这些新的注释进行校正  
该流程必须按照 GVCF 模式进行  
workflow：  
* 生成 GVCF 文件时加上 `-G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation` 参数  
* 合并 GVCF 文件时加上 `-G StandardAnnotation -G AS_StandardAnnotation` 参数  
* 联合识别时加上 `-G StandardAnnotation -G AS_StandardAnnotation` 参数  
* VQSR 时加上 `-AS` 参数
  
## 评估变异集的质量  
实验层面的方法包括 Sanger 测序和芯片等，这里主要考虑计算层面的。通过比较金标准数据和自己的变异集的多个度量，来评估变异的可靠性，以及可能找出错误的来源。使用该方法的要求是（1）真实集中的数据都是真实可靠的；（2）自己变异集样本的基因组与真实集之间应当类似  
两个衡量指标：  
* SNPs/Indels 的数量，WGS 为4.4M，WES 为41k  
* TiTv Ratio，transition (Ti)，转换，嘌呤和嘌呤，或者嘧啶和嘧啶之间的转换；transversion (Tv)，颠换，嘌呤和碱基之间的转换。如果不考虑生物学影响，这个比例应该是0.5，但实际中甲基化的胞嘧啶容易脱氨基转换成胸腺嘧啶；另外像 CpG island，多位于引物区域，具有大量甲基化胞嘧啶，所以包含这些区域的 WES ti/tv 的值会更高。一般 WGS 为2.0-2.1，WES 为3.0-3.3，另外就是对于外显子组测序，该比例受旁侧序列长度影响  
* Insertions/Deletions，按照目的进行过滤之后，如果寻找常见变异，值约为1，如果是罕见变异，则为0.2-0.5  

# REF
1. https://gatkforums.broadinstitute.org/gatk/discussion/11165/data-pre-processing-for-variant-discovery  

## BQSR  
1. https://software.broadinstitute.org/gatk/documentation/article?id=11081  

## HaplotypeCaller 
1. https://software.broadinstitute.org/gatk/documentation/article?id=11077  
2. https://software.broadinstitute.org/gatk/documentation/article?id=3893  


## Hard filter
1. https://software.broadinstitute.org/gatk/documentation/article?id=11069  

## VQSR
1. https://software.broadinstitute.org/gatk/documentation/article?id=11084  
2. https://software.broadinstitute.org/gatk/documentation/article.php?id=1259  
3. https://software.broadinstitute.org/gatk/documentation/article.php?id=2805  
4. https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_vqsr_ApplyVQSR.php  
5. https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php  

## 质量分数  
1. https://software.broadinstitute.org/gatk/documentation/article?id=11024  
2. https://software.broadinstitute.org/gatk/documentation/article?id=11075  
3. https://software.broadinstitute.org/gatk/documentation/article?id=11079  
  
## 深度  
1. https://software.broadinstitute.org/gatk/documentation/article?id=11072  
  
## 变异集好坏的评估  
1. https://software.broadinstitute.org/gatk/documentation/article?id=11071  
