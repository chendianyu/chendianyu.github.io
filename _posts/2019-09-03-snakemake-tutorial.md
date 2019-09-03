---
title: Snakemake 学习笔记
description: 
categories:
 - Workflow
tags:
 - workflow
 - pipeline
 - tutorial
---

# Setup
利用`Conda`安装`Snakemake`  

```shell
conda install -c bioconda -c conda-forge snakemake
```
  
# Basic
`Snakemake workflow`可以通过在`Snakefile`中指定规则来定义。规则指明了如何将基于输入文件构造输出文件，将整个流程划分成一个个小步骤。`Snakemake`能够基于文件名的匹配情况自动确定`rule`之间的依赖关系。`Snakemake`语言在`Python`的基础上进行了拓展，添加了规则定义以及其他用于控制的语句结构。**所有添加的语法结构都是以一个关键字起始，后面跟着代码块**  
接下来将以变异检测流程为例，介绍Snakemake的语法  
  
## Step 1: Mapping reads
在 Snakefile 中定义以下规则： 

```shell
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/A.fastq"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

* 每个rule会有一个名字（本例中为`bwa_map`）和一组指令（本例中包括`input`, `output`和`shell`）  
* `input`和`output`指令之后会跟着一组文件名，**注意多个文件时要用逗号分隔**。最简单的情况就是像本例中直接为显式的python字符串  
* `shell`指令后的字符串包含了要执行的shell命令。在命令中如果需要指向规则规则内的元素，那么就用花括号括起来即可。本例中由于有多个输入文件，Snakemake会将它们连接在一起，用空格分隔，即用`data/genome.fa data/samples/A.fastq`替换`{input}`

构建好Snakefile之后就可以执行  

```shell
snakemake -np mapped_reads/A.bam
## 在包含Snakefile的工作目录中执行该命令
## -n --dryrun 表示伪执行，仅展示执行计划
## -p 表示同时打印shell命令便于理解
## mapped_reads/A.bam 指定要生成的 target files
```

为了生成目标文件，Snakemake以从上往下的方式执行Snakefile中的规则，规则的执行称为`job`  
  
## Step 2: Generalizing the read mapping rule
Snakemake能够通过使用通配符实现更加普遍性的规则，如下：  

```shell
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
```

本例中用`{sample}`取代了上例中的`A`，Snakemake 会将所有符合`*.fastq` 的文件处理成`*.bam`。在文件路径名称中可以有多个通配符，不过为了避免与同一规则内其他的job发生冲突，要求同一规则内所有输出文件含有相同的通配符  
执行命令  

```shell
snakemake -np mapped_reads/{A,B,C}.bam
## 等价于 snakemake -np mapped_reads/A.bam mapped_reads/B.bam mapped_reads/C.bam
## 告诉snakemake要生成的target files，就会用A、B、C去替换{sample}
## 本例中不会对A在进行运行，因为mapped_reads/A.bam文件已经存在了
## 除非存在有输入文件比某个输出文件更新或者有输入文件被其他的job更新才会对其重新运行一遍流程
```
  
## Step 3: Sorting read alignments
生成 BAM 文件之后我们会对其进行排序，在Snakefile中添加如下语句：  

```shell
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

Snakemake会在执行job之前自动创建缺失的目录。本例中因为用了`-T`参数，需要指定中间临时文件的前缀名，即我们需要用到通配符{sample}对应的值，Snakemake允许通过`wildcards`对象在shell命令中调用这些值。该对象中会有对应的属性指向每一个通配符的值  
通过执行`snakemake -np sorted_reads/B.bam`，Snakemake会先执行bwa_map，然后执行samtools_sort。如之前所说，Snakemake会根据文件名自动解决依赖关系  
  
## Step 4: Indexing read alignments and visualizing the DAG of jobs
接下来需要对排序后的BAM文件构建索引：  

```shell
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
```

如果需要构建流程的有向非循环图（directed acyclic grap, DAG），可以在调用snakemake命令的时候指定，将会调用`Graphviz`中的`dot`命令来画图：  

```shell
snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg
## snakemake用dot语言指明要怎么画图，传递给dot命令，然后被渲染成SVG格式
```

DAG图中，每个点代表一个job，每条边表示依赖关系。如果job没有运行（比如已经存在结果文件了），那么会用虚线表示  
  
## Step 5: Calling genomic variants
这一步中我们需要用到samtools和bcftools 并整合所有样本中比对上的reads来call variants。当然可以一个个显示地在input中定义出来：  

```shell
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bamA="sorted_reads/A.bam"
        bamB="sorted_reads/B.bam"
        baiA="sorted_reads/A.bam.bai"
        baiB="sorted_reads/B.bam.bai"
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bamA} {input.bamB} | "
        "bcftools call -mv - > {output}"
```

但这显然非常麻烦，所以考虑用列表推导式的方式来实现.在规则外定义`SAMPLES=["A","B"]`，然后在规则内使用`bam=["sorted_reads/{}.bam".format(sample) for sample in samples]`。由于列表推导式比较常用，所以Snakemake定义了`expand` 来进行简化，最后结果就是：  

```shell
SAMPLES=["A","B"]
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```

对于较长的shell命令，建议像代码中那样将字符串分成多行，python会自动将它们合为一行  
  
## Step 6: Using custom scripts
有时候流程中会包含自己写的脚本，比如代码较长时我们会把它写入独立的脚本文件中然后去调用这个脚本。Snakemake提供了`script`指令来运行这些脚本  

```python
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

脚本的**路径总是相对于Snakefile而言**。在脚本中，规则中的所有内容，比如input，output以及wildcards等都会是全局`snakemake`对象的属性，例如脚本可以这样写：  

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysam import VariantFile

quals = [record.qual for record in VariantFile(snakemake.input[0])]
plt.hist(quals)

plt.savefig(snakemake.output[0])
```

虽然我们可以通过比如shell命令来调用独立的脚本，但是像上面代码那样操作可以使脚本的逻辑独立于流程的逻辑  
除了python，我们还可以使用R脚本。类似的，会生成一个名为`snakemake`的S4对象，比如可以通过`snakemake@input[[1]]`来访问第一个输入文件  
  
## Step 7: Adding a target rule
之前执行流程的时候都会在命令行指定目标文件。除了文件名，当指向的规则不存在通配符时，Snakemake 也接受规则名作为目标，这样我们可以通过目标规则收集想要的结果。如果没有在命令行指定目标的话，Snakemake会将Snakefile的第一条规则定位目标，所以我们可以将名为`all`的规则放在最上部，该规则将所有想要的结果文件作为输入：  

```python
rule all:
    input:
        "plots/quals.svg"
```

之后只需在命令行中执行`snakemake -n`即可。此外，如果在Snakefile顶部添加了多条命令规则，可以通过`snakemake -n mytarget`指定输出哪些结果  
  
# Advanced
## Step 1: Specifying the number of used threads
有些工具建议使用多个线程进行加速，Snakemake 提供了`threads`指令：  

```python
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

执行流程时，Snakemake通过调度程序来调控任务需要的线程数，并确保正在运行的任务使用的总线程数不超过给定的可用核心数，后者可在命令行中指定（默认为1）`snakemake --cores 10`  
像上面的代码中，跑bwa_map时需要8个线程，所以这样的任务一次只能跑一个。但是Snakemake调度程序会尝试去跑其他的任务来让剩余可用核心数也利用起来，比如在一个比对跑完之后，同时跑下一个比对和samtools_sort，保证效率最大化  
  
## Step 2: Config files
在上面的代码中，我们是将要跑的样本写在Snakefile，这样的话对于不同的项目，我们每次需要对Snakefile进行修改。对于这种情况，Snakemake提供了配置文件机制。配置文件可以用`JSON`或`YAML`语法，例如下面这个`config.yaml`文件：  

```shell
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

而Snakefile文件就可以写成：  

```python
configfile: "config.yaml"
## Snakemake会加载配置文件并将文件中的内容放入名为config的全局字典中
## 形式为{'samples': {'A': 'data/samples/A.fastq', 'B': 'data/samples/B.fastq'}}

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["smaples])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
## 虽然sample是个字典，但在展开时只会使用key的部分，即A和B
```

## Step 3: Input functions
既然将FASTQ文件的路径也放在了配置文件中，那么我们也可以一般化bwa_map规则来使用这些路径。为了搞清楚与上一步之间的差异，首先我们需要理解Snakemake流程执行的3个阶段：  
1. initialization阶段：解析流程并将所有规则实例化  
2. DAG阶段：用具体的值填补通配符，输入文件和输出文件对上号，构建任务的DAG图  
3. scheduling阶段：执行任务的DAG  

所以上一步中expand函数是在第一阶段执行的。在这一阶段，对于任务是什么，通配符的值以及规则之间的依赖关系我们都不知道，所以也就没办法根据配置文件来决定bwa_map规则里FASTQ文件的路径，因为我们甚至不知道该规则会生成什么样的任务。我们需要将输入文件的确定推迟至DAG阶段，这可以通过在输入中指定一个函数而不是直接使用字符串来实现：  

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"
```

这里通过匿名函数（也可以用其他常规的函数）将wildcards对象作为单变量输入，然后通过属性（wildcards.sample）获取通配符对应的具体值。返回的结果为字符串或字符串列表  
  
## Step 4: Rule parameters
有时候shell命令不知包含输入、输出以及静态的标记，尤其是有些参数需要根据任务的通配符来确定。对此，Snakemake可通过`params`指令实现在规则中定义任何参数：  

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"
```

## Step 5: Logging
在执行大型流程时，往往需要将各步骤的输出放到一个文件中，对此可通过`log`指令实现：  

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        "mapped_reads/{sample}.bam"
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

shell命令会收集bwa和samtools的STDERR结果输出至log文件中。log文件和输出文件需要用相同的通配符，以避免冲突  
  
## Step 6: Temporary and protected files
Snakemake允许将输出文件标记为临时文件，这样对应的任务执行完之后会将他们杀出，以节省硬盘空间：  

```python
rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"
```

这样就会在执行完后续的samtools_sort任务后将BAM文件进行删除  
不过由于比对耗费资源大，为了避免最终的BAM文件被意外删除或修改，我们可以对其进行保护：  

```python
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"
```

在执行完该任务之后，Snakemake会对文件系统中的输出文件进行写入保护，使得它无法被覆盖或删除  
  
## Summary
```python
configfile: "config.yaml"

rule all:
    input:
        "plots/quals.svg"

rule bwa_map:
    input:
        "data/genome.fa",
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    params:
        rg=r"@RG\tID:{sample}\tSM:{sample}"
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 8
    shell:
        "(bwa mem -R '{params.rg}' -t {threads} {input} | "
        "samtools view -Sb - > {output}) 2> {log}"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        protected("sorted_reads/{sample}.bam")
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```
