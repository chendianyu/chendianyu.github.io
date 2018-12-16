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

# Data pre-processing for variant discovery  
## Reference  
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

## Expected input
数据最初是按照不同的 `readgroups` 分成不同的子集，对应 `libraries`（DNA 取自不同的生物学样本，以及在文库制备过程中的片段化和 barcode 标记过程）与  `lanes` (测序仪的物理分隔单元) 通过 `multiplexing`（混合多个文库，并在多条泳道上进行测序，以减少风险和人为误差）交叉得到的不同组合  
Read groups 在 SAM/BAM/CRAM 中通过一系列标签进行定义。在 SAM 等文件的 header 中，以 `@RG` 起始，后面跟着以下几个 tag：  
* ID：每个 read group 的 ID 是唯一的。对于 Illumina 数据，一般就是 flowcell + lane name and number。对于 BQSR 等处理，每一个 read group 被认为是独立的，被认为共享相同的错误模型  
* PU：Platform Unit，包含 {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE} 信息。{FLOWCELL_BARCODE} 是 flowcell 的唯一标识，{LANE} 表示 flow cell 中的 lane，而 {SAMPLE_BARCODE} 则是 sample/library-specific identifier。PU 不是 GATK 必需的，但如果出现，将优先于 ID 进行 BQSR 时的分组  
* SM：该 read group 中测序的样本的名称。GATK 将所有具有相同 SM 的数据视为来源于同一样本，这也将会是 VCF 文件中 sample column 所使用的名称  
* PL：平台或使用的技术名称。Valid values: ILLUMINA, SOLID, LS454, HELICOS and PACBIO  
* LB：文库标识符。当同一 DNA 文库在多条 lane 进行测序时，`MarkDuplicates` 会基于 LB 来确定哪些 read group 可能含有分子重复  


# Main steps  
## Map to Reference  
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
  
## Coordinate sort and index  
Tools involved: Picard's `SortSam` and `BuildBamIndex`  
`SortSam` can sort the input SAM or BAM file by **coordinate**, **queryname (QNAME)**, or **some other property** of the SAM record. The SortOrder of a SAM/BAM file is found in the SAM file header tag `@HD` in the field labeled `SO`  
(For a coordinate sorted SAM/BAM file, read alignments are sorted first by the reference sequence name (RNAME) field using the reference sequence dictionary (@SQ tag). Alignments within these subgroups are secondarily sorted using the left-most mapping position of the read (POS). Subsequent to this sorting scheme, alignments are listed arbitrarily  
For queryname-sorted alignments, all alignments are grouped using the queryname field but the alignments are not necessarily sorted within these groups. Reads having the same queryname are derived from the same template)  
  
This creates a file called `reads_sorted.bam` containing reads sorted by genomic location, aka coordinate, and a `.bai` index file with the same prefix as the output, e.g. `reads_sorted.bai`, within the same directory. Index file for the input BAM that allows fast look-up of data in a BAM file, lke an index on a database. Note that this tool **cannot be run on SAM files**, and that the input BAM file must be sorted in coordinate order  
Usage :  
```shell
java -jar picard.jar SortSam \
    INPUT=<input.sam> \
    OUTPUT=<output_sorted.bam> \ 
    SORT_ORDER=coordinate \
    CREATE_INDEX=true
```  
  
You can also use `BuildBamIndex` to create index  
```shell  
java -jar picard.jar BuildBamIndex \
      I=<input.bam>  
```  
  
OR `samtools index`  
Usage: `samtools index <input.bam>` 产生的文件为 `<input.bam.bai>` 只有这个与picard有区别，文件内容本质上应该是一致的  

## Mark Duplicates
Tools involved: `MarkDuplicates`  
Duplicates can arise during sample preparation e.g. library construction using PCR. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as `optical duplicates`  
This processing step is performed `per-sample` and consists of identifying read pairs that are likely to have originated from duplicates of the same original DNA fragments through some artifactual processes. The program tags all but of the read pairs within each set of duplicates, causing them to be ignored by default during the variant discovery process. This step constitutes a major bottleneck since it involves making a large number of comparisons between all the read pairs belonging to the sample, across all of its readgroups  
The tool's main output is **a new SAM or BAM file**, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of `0x0400`, which corresponds to a decimal value of 1024  
(To identify the type of duplicate, a new tag called the `duplicate type (DT)` tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. Invoking the `TAGGING_POLICY` option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag, as either `library/PCR-generated duplicates (LB)`, or `sequencing-platform artifact duplicates (SQ)`)  
  
The following commands take a coordinate-sorted and indexed BAM and return:  
(i) a BAM with the same records in coordinate order and with duplicates marked by the 1024 flag,  
(ii) a duplication metrics file,   
(iii) an optional matching BAI index  
  
Usage : 
```shell
java -jar picard.jar MarkDuplicates \
    INPUT=<input_sorted.bam> \
    OUTPUT=<output_sorted_duplicates.bam> \ 
    METRICS_FILE=<output_markduplicates_metrics.txt> \ 
    CREATE_INDEX=true
```  
  
(Maybe you need Picard's `FixMateInformation` to verify mate-pair information between mates and fix if needed)  
Usage:  
```shell
java -jar picard.jar FixMateInformation \ 
    I=<input_sorted_duplicates.bam> \ 
    O=<output_sorted_duplicates_fixed.bam> \ 
    ADD_MATE_CIGAR=true  
```
  
## Base (Quality Score) Recalibration
Tools involved: `BaseRecalibrator`, `Apply Recalibration`, `AnalyzeCovariates` (optional)
This third processing step is performed `per-sample` and consists of applying machine learning to detect and correct for patterns of systematic errors in the base quality scores, which are confidence scores emitted by the sequencer for each base  
  
Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. `Base quality score recalibration (BQSR)` is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls  
The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. The known variants are used to mask out bases at sites of real (expected) variation, to avoid counting real variants as errors. Outside of the masked sites, every mismatch is counted as an error. The rest is mostly accounting  
This process is accomplished by analyzing the covariation among several features of a base. For example:  
* Reported quality score  
* The position within the read  
* The preceding and current nucleotide (sequencing chemistry effect) observed by the sequencing machine  

## Steps  
1. Analyze patterns of covariation in the sequence dataset  
```shell
gatk BaseRecalibrator \ 
    -I <input_sorted_markduplicates.bam> \
    -R <ref.fasta> \ 
    -knownSites <dbsnp.vcf> \ 
    -knownSites <1000g.vcf> \ 
    -o <recal_data.table>  
```
2. Do a second pass to analyze covariation remaining after recalibration (Optional)  
```shell
java -jar GenomeAnalysisTK.jar \ 
    -T BaseRecalibrator \ 
    -R <ref.fasta> \ 
    -I <input_sorted_markduplicates.bam> \
    -knownSites <dbsnp.vcf> \ 
    -knownSites <1000g.vcf> \ 
    -BQSR <recal_data.table> \ 
    -o <post_recal_data.table> 
```
* the `-BQSR` flag tells the GATK engine to perform on-the-fly recalibration based on the first recalibration data table  
3. Generate before/after plots (Optional)  
```shell
java -jar GenomeAnalysisTK.jar \ 
    -T AnalyzeCovariates \ 
    -R <ref.fasta> \ 
    -before <recal_data.table> \
    -after <post_recal_data.table> \
    -plots <recalibration_plots.pdf>
```
4. Apply the recalibration to your sequence data  
```shell
java -jar GenomeAnalysisTK.jar \ 
    -T PrintReads \ 
    -R <ref.fasta> \ 
    -I <input_sorted_markduplicates.bam> \ 
    -BQSR <recal_data.table> \ 
    -o <output_sorted_markduplicates_recal.bam> 
```
* By default, the original quality scores are discarded in order to keep the file size down
  
# REF
1. https://gatkforums.broadinstitute.org/gatk/discussion/11165/data-pre-processing-for-variant-discovery
