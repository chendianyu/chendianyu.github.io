---
title: GATK4 WES/WGS 分析流程
description: GATK4 检测细胞系突变流程
categories:
 - GATK
tags:
 - GATK4
 - software
 - pipeline
 - SNP
 - Indel
---

# Introduction  
This pipeline uses `GRCh38` as the reference genome. It begins with unaligned paired reads in BAM format and results in a sample-level SNP and INDEL variant callset in GVCF format  

# Data pre-processing for variant discovery  
## Reference  
GRCh38 reference genome FASTA including `alternate contigs`, corresponding `.fai` index and `.dict` dictionary and five BWA-specific index files `.sa`, `.amb`, `.bwt`, `.ann` and `.pac`  

GATK uses two files to access and safety check access to the reference files:  
1. a `.dict` dictionary of the contig names and sizes  
You can use `Picard` `CreateSequenceDictionary` to create a `SAM/BAM file` from a fasta containing reference sequence. The reference sequence can be gzipped (both `.fasta` and `.fasta.gz` are supported). The output SAM file contains **a header but no SAMRecords**, and the header contains only sequence records  
  
Usage:  
```shell
java -jar picard.jar CreateSequenceDictionary \ 
      R=<ref.fasta> \ 
      O=<ref.dict>  
```  
e.g.:  
```shell  
java -jar picard.jar CreateSequenceDictionary \ 
      R=ucsc.hg38.fasta \ 
      O=ucsc.hg38.dict  
```  

2. a `.fai` fasta index file to allow efficient random access to the reference bases  
You can use `Samtools` `faidx` to index the reference sequence file and create `<ref.fasta>.fai` on the disk. The input file can be compressed in the BGZF format  
  
Usage : `samtools faidx <ref.fasta>`  
e.g. : `samtools faidx ucsc.hg38.fasta`  
  
BWA alignment requires an indexed reference genome file. Indexing is **specific to algorithms**. To index the human genome for BWA, we apply BWA's `index function` on the reference genome file. This produces five index files with the extensions `amb`, `ann`, `bwt`, `pac` and `sa`  
Usage : `bwa index -a bwtsw <ref.fasta>`  
e.g. : `bwa index -a bwtsw /home/refer/ucsc.hg38.fasta`  
The tool automatically locates the index files within the `same folder` as the reference FASTA file  

当然也可以进入GATK FTP下载已有的文件  
```shell
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.fai & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.dict & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz & 
nohup wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi   
```  

## Expected input
The data is initially organized in distinct subsets called `readgroups`. These correspond to the intersection of `libraries` (the DNA product extracted from biological samples and prepared for sequencing, which includes fragmenting and tagging with identifying barcodes) and `lanes` (units of physical separation on the DNA sequencing chips) generated through `multiplexing` (the process of mixing multiple libraries and sequencing them on multiple lanes, for risk and artifact mitigation purposes)  

# Main steps  
## Map to Reference  
Tools involved: `BWA`, Picard's `MergeBamAlignments`
This first processing step is performed `per-read group` and consists of mapping `each individual read pair` to the reference genome, which is a synthetic single-stranded representation of common genome sequence that is intended to provide a common coordinate framework for all genomic analysis. Because the mapping algorithm **processes each read pair in isolation**, this can be massively parallelized to increase throughput as desired  

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
