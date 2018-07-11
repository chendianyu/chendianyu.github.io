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
pybedtools 的安装要求有以下软件：  
1. Python 2.7 及更高的版本  
2. bedtools，最好是最新的，且位于 PATH 中  
3. C/C++ 编译器，Linux 中就是 gcc，通常已经安装好了  

最简单的安装方法是通过 Anaconda 进行安装：  
```shell
conda install -c bioconda pybedtools
conda install -c bioconda bedtools htslib
```  
另外也可通过 pip 进行安装：  
```shell
pip install pybedtools  
```  
  
# Convention  
* integer values 一直都是 0-based，而不考虑文件格式，这就意味着所有的 `Interval` 对象的起止位置值都可以统一对待  
* string values 则根据对应的文件格式的要求进行调整，原文件是什么样就什么样，哪怕表示的是起止坐标  
  
# Design Principles 
## 自动创建和删除临时文件  
由于 BedTool 对象必须指向硬盘上的某个文件，因此每当我们对 BedTool 对象进行处理得到新的对象时，都会创建出一个临时文件（通过 python 的 tempfile 模块处理），其路径可通过 `.fn` 属性得到。  
临时文件默认存放于 `/tmp` 目录下（由 tempfile 模块指定），命名模式为 `/tmp/pybedtools.*.tmp`。当退出本次操作时，这些临时文件也会被自动删除。但如果非正常退出，就可能导致这些文件未被清理，可进行手动清理。  
如果想在会话中删除之前所创建的临时文件，可通过 `pybedtools.cleanup()`；如果使用 `pybedtools.cleanup(remove_all=True)`，则不管是不是本次会话中创建的，只要文件名符合 `<tempdir>/pybedtools.*.tmp`模式，均会被删除，其中 `<tempdir>` 是 `pybedtools.get_tempdir()` 的值。  
如果想要自己指定临时文件存放的目录，可通过 `pybedtools.set_tempdir(<path>)` 进行设置。  
  
## 取决于 bedtools 版本
由于 BedTool 各种方法是对 bedtools 程序的封装，因此所能使用的各种方法以及参数取决于你当前安装的 bedtools 版本  
  
## 默认参数的设置  
```python
>>> # 1. 单个输入文件时  
>>> result1 = a.merge(d=100, s=True)
>>> result2 = a.merge(i=a.fn, d=100, s=True)  
>>> result1 == result2
True  
>>> 
>>> # 2. 多个输入文件时  
>>> result3 = a.intersect(b)
>>> result4 = a.intersect(a=a.fn, b=b.fn)    # 可以是文件名，也可以是 BedTool 对象
>>> result3 == result4
True
```
除了指向 Bed 文件的参数（`-i`, `-a`, `-b`）会根据输入智能化地确定对应的值外，其余的参数均为关键字参数，不设定默认值，不显式指定的话取决于所安装的 beedtools 所设定的默认值。  
  
## 串联命令  
大部分的方法返回的结果是 BedTool 对象，可以直接继续调用方法，从而实现多个步骤串联起来，相当于 shell 中使用管道命令  
```python  
>>> a.intersect(b).merge().saveas('shared_merged.bed')  
>>> # 等价于 intersectBed -a a.bed -b b.bed | merge -i stdin > shared_merged.bed  
```  

# BedTool Object
一般来说，单个 `BedTool` 指向单个区间文件，可以是 BED, GFF, GTF, VCF, SAM, or BAM format 或者对应的压缩格式等。`BedTool` 对象封装了所用可用的 bedtools 程序，使之可以在 python 中进行调用。  
  
## BedTool 对象创建
```python 
# 1. 从文件中读取
import pybedtools
bed_file = pybedtools.BedTool(<bed_file>）     

# 2. 直接读取字符串
>>> # using a longer string to make a bed file.  Note that
>>> # newlines don't matter, and one or more consecutive
>>> # spaces will be converted to a tab character.
>>> larger_string = """
... chrX 1    100   feature1  0 +
... chrX 50   350   feature2  0 -
... chr2 5000 10000 another_feature 0 +
... """

>>> fromscratch = pybedtools.BedTool(larger_string, from_string=True)
>>> print(fromscratch)
chrX    1   100 feature1    0   +
chrX    50  350 feature2    0   -
chr2    5000    10000   another_feature 0   +
```  
  
## BedTool 结果保存   
```python 
>>> # 1. .saveas() 方法除了保存结果外，还会对结果进行拷贝，即返回指向所保存的文件的 BedTool 对象，
>>> # 因此对于大型文件会比较耗时间和内存，但适合数据流操作，常用于在流程中间将结果保存到硬盘上
>>> # 此外，该方法还可以给文件加上 track line，方便后续传到 UCSC Genome Browser
>>> result = a.each(TSS, upstream=1000, downstream=0).saveas('upstream_regions.bed')  
>>>
>>> # 2. .moveto() 方法是将结果进行重命名，不进行拷贝，速度更快。
>>> # 该方法适合已经将结果写到硬盘上，比如临时文件，但希望给该文件一个更适合的名字进行保存 
>>> c = a.intersect(b).moveto('intersection_of_a_and_b.bed')
>>> 
>>> # 3. 此外，对所有作用于 BedTool 并得到新的 BedTool 对象的方法，
>>> # 均有 output 参数，可直接保存结果，覆盖默认创建的临时中间文件的操作  
>>> c = a.intersect(b, output='intersection_of_a_and_b.bed')
```
  
# Interval Object  
`Interval` 对象是 pybedtools 表示 bed 等文件中一行数据时所使用的。可以通过索引和切片操作获取 BedTool 中的 Interval 对象  
```python  
>>> feature = a[0]
>>> features = a[1:3]  
```  
  
对于每个 Interval 对象，使用 print() 进行打印时输出的就是文件中原来的行，所以对于不同坐标标准的文件会自动变化  
对于所有的 features, 不管最初存储的文件格式是哪一种（可通过 `.file_type` 属性获得）, 转成 Interval 后均有 `chrom`, `start`, `stop`, `name`, `score` 和 `strand` 属性：其中 `start` 和 `stop` 是 integer，均转成了 0-based, 其他的包括 `score` 均为 string，与原文件对应一致；如果对应的值缺失，标记为 `'.'`  
  
除了直接利用属性来获取对应字段之外，也可以利用位置或名字进行索引得到各个字段  
```python  
>>> feature[0]  
'chr1'  
>>> feature['start']
1
```
  
`.fields` 属性存储的是由原行切割后的 string 形式（包括坐标）构成的 list，所以在进行位置索引时实际上是对 fields 属性的 list 进行索引，所使用的位置索引值是其在原文件中的位置  
```python
>>> feature.fields  
['chr1', '1', '100', 'feature1', '0', '+']
```  
进行位置索引时最后得到的是 string 形式，那么对于起止坐标值而言，像 gff 等格式就会导致存在差异了，这点需要注意，可以注意下后面的实例  

通过 `len()` 函数可以获得 Interval 对象的长度  
  
## 实例  
```python
>>> gff = ["chr1",
...        "fake",
...        "mRNA",
...        "51",   # <- start is 1 greater than start for the BED feature below
...        "300",
...        ".",
...        "+",
...        ".",
...        "ID=mRNA1;Parent=gene1;"]
>>> gff = pybedtools.create_interval_from_list(gff)
>>>
>>> bed = ["chr1",
...        "50",
...        "300",
...        "mRNA1",
...        ".",
...        "+"]
>>> bed = pybedtools.create_interval_from_list(bed)
>>>
>>> gff.file_type
'gff'
>>> bed.file_type
'bed'
>>>
>>> bed.start
50
>>> bed[1]          # 注意到返回结果是个字符串，不过由于 bed 文件本身就是 0-based 的，所以值没有变化
'50'
>>>
>>> gff.start
50
>>> gff['start']    # 用名字进行索引得到的还是整数，不受影响
50
>>> gff[3]          # 由于原文件是 gff 格式，1-based，因此就出现差了1的情况      
'51'
```
  
## GFF，GTF attributes
GFF 和 GTF 每行的最后一个字段都包含有大量有用的信息，这些信息将会以字典的形式存储在 `.attrs` 属性中，因此可以通过字典操作去对这些信息进行处理，包括添加，修改，删除等操作  
```python
>>> # original feature
>>> print(gff)
chr1        fake    mRNA    51      300     .       +       .       ID=mRNA1;Parent=gene1;

>>> # original attributes
>>> sorted(gff.attrs.items())
[('ID', 'mRNA1'), ('Parent', 'gene1')]

>>> # add some new attributes
>>> gff.attrs['Awesomeness'] = "99"
>>> gff.attrs['ID'] = 'transcript1'

>>> # Changes in attributes are propagated to the printable feature

>>> # for testing, we make sure keys are sorted. Not needed in practice.
>>> gff.attrs.sort_keys = True
>>> assert gff.attrs.sort_keys
>>> print(gff)
chr1        fake    mRNA    51      300     .       +       .       Awesomeness=99;ID=transcript1;Parent=gene1;
```
  
# 过滤
`.filter()` 方法允许你提供一个函数，该函数能够接受一个 Interval 对象作为其第一个参数，并返回 True 或者 False，从而实现对 BedTool 对象进行类 grep 的操作。filter() 方法会将提供给它的 `args` 和 `kwargs` 均传递给函数  
```python  
>>> # 匿名函数
>>> a = pybedtools.example_bedtool('a.bed')
>>> b = a.filter(lambda x: len(x) > 100)
>>> print(b)
chr1        150     500     feature3        0 
>>>
>>> # 更一般化的函数  
>>> def len_filter(feature, L):
...     "Returns True if feature is longer than L"
...     return len(feature) > L  
>>> a = pybedtools.example_bedtool('a.bed')
>>> print(a.filter(len_filter, L=10))
chr1        1       100     feature1        0       +
chr1        100     200     feature2        0       +
chr1        150     500     feature3        0       -
chr1        900     950     feature4        0       +
```
  
另外，在 `featurefuncs` 模块中包含一些用 `Cython` 预先写好的一些函数，能够实现更快的速度，例如 `greater_than()` 和 `less_than()` 等  
```python
>>> from pybedtools.featurefuncs import greater_than
>>> a.filter(greater_than, 100)  
```
  
# Each
类似之前的 `filter()` 方法，`.each()` 方法也会对每个 feature 应用函数，并返回一个新的，可能经过修饰的 Interval  
```python  
>>> def normalize_count(feature, scalar=0.001):
...     """
...     assume feature's last field is the count
...     """
...     counts = float(feature[-1])
...     normalized = round(counts / (len(feature) * scalar), 2)
...
...     # need to convert back to string to insert into feature
...     feature.score = str(normalized)
...     return feature
>>> normalized = with_counts.each(normalize_count)
>>> print(normalized)
chr1        1       100     feature1        0.0     +       0
chr1        100     200     feature2        10.0    +       1
chr1        150     500     feature3        2.86    -       1
chr1        900     950     feature4        20.0    +       1
```
  
# 基因组  
对于 bed 的操作，有时候需要提供基因组或染色体大小信息以避免得到的结果超过染色体的范围。`genome_registry` 模块中已经准备了一些常用基因组各染色体大小的信息，存储为有序字典，key 为染色体名称，value 则是 0-based 的起止坐标的 tuple  
```python
>>> from pybedtools import genome_registry
>>> genome_registry.dm3['chr2L']
(0, 23011544)
```  
  
如果要在方法中指定基因组，可以是文件名，关键字或者自己提供的字典。有两个参数可以调用：
* `g` 文件名或者字典  
* `genome` 关键字或者字典  
  
```python
>>> a = pybedtools.example_bedtool('a.bed')
>>> b = a.slop(b=100, g='hg19.genome') 
>>> c = a.slop(b=100, genome='hg19')    # 首先在 genome registry 中找，没有的话直接去 UCSC 中下载  
>>> d = a.slop(b=100, g={'chr1':(1, 10000)})
>>> e = a.slop(b=100, genome={'chr1':(1,100000)})
```  

# 操作符重载  
为了方便起见，pybedtools 将 `+` 和 `-` 重载，以方便进行 intersection  
* `+` 等价于 `intersect` 加上 `u` 参数 `a+b <=> a.intersect(b, u=True)`  
* `-` 等价于 `intersect` 加上 `v` 参数 `a-b <=> a.intersect(b, v=True)`  
  
# 比较  
有时候我们会对 feature 进行比较以确定它们是否存在重叠或者相对位置是什么样的，就会用到比较操作符 `<`, `<=`, `==`, `=>`, `>`，返回 True 或者 False。需要注意的是，这种比较会忽略染色体链。  
以下是这些操作符的示意图  
```python
# a == b, a >= b, a <= b
a ---------
b ---------

# a < b, a <= b
a ----
b       -----

# a <= b
a ----
b     -----  (book-ended)

# a >= b
a     -----
b ----      (book-ended)

# a > b, a >= b
a       ------
b ----

# a >= b
a  ------------
b  ---------

# a >= b
a   -----------
b -------------

# a <= b
a -------------
b   -----------

# a <= b
a  ---------
b  ------------

# a <= b
a -----------
b        -----------

# a >= b
a        -----------
b -----------

# undefined!
a    ----
b -----------

# undefined!
a -----------
b    ----

# a <= b
a -----------
b           -

# a >= b
a           -
b -----------

# a == b, a <= b, a >= b
a -
b -  (starts and stops are identical for all features)
```
