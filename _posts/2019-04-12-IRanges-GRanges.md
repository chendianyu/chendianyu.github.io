---
title: IRanges and GRanges
description: IRanges, GRanges介绍
categories:
 - Bioconductor
tags:
 - bioconductor
---

# IRanges
`IRanges` 对象用于描述范围，其有3个基本的属性：  
* `start`  
* `end`  
* `width`  

我们可以通过 `IRanges()` 函数构建 `IRanges` 对象，需要提供上述3个属性中的任意两个:  
```R
library(IRanges)
## 下面3个命令得到的结果是相同的
ir1 <- IRanges(start=1:10, width=10:1)
ir2 <- IRanges(start=1:10, end=11)
ir3 <- IRanges(end=11, width=10:1)

## 获得3个属性的值
## 也可以通过赋值修改这些值
start(ir1)
end(ir1)
width(ir1)
```  
  
对 `IRanges` 对象取子集：  
```R
> ir <- IRanges(c(1, 8, 14, 15, 19, 34, 40),
          width=c(12, 6, 6, 15, 6, 2, 7))
> ir
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15
  [5]        19        24         6
  [6]        34        35         2
  [7]        40        46         7

## 数值索引
> ir[1:4]
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15

## 逻辑索引
> ir[start(ir) <= 15]
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15
```
  
对每个范围值添加名称，一是通过 `IRanges()` 函数中的`names` 参数，而是构建完之后通过 `names()` 函数赋值  
```R
names(ir) <- letters[1:7]
```
  
分析过程中有时我们需要操作一组 `IRanges` 对象，所以提供了 `IRangesList` 类，当中的每个元素均为 `IRanges` 对象。可以通过 `start`，`end` 和 `width` 函数获取 `IntegerList` 对象，记录对应的属性值  
```R
> rl <- IRangesList(ir, rev(ir))
> rl
IRangesList of length 2
[[1]]
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15
  [5]        19        24         6
  [6]        34        35         2
  [7]        40        46         7

[[2]]
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        40        46         7
  [2]        34        35         2
  [3]        19        24         6
  [4]        15        29        15
  [5]        14        19         6
  [6]         8        13         6
  [7]         1        12        12

> start(rl)
IntegerList of length 2
[[1]] 1 8 14 15 19 34 40
[[2]] 40 34 19 15 14 8 1
```
  
## 对IRanges对象的操作
### 合并
```R
> a <- IRanges(start=7, width=4)
> b <- IRanges(start=2, end=5)
> c(a, b)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         7        10         4
  [2]         2         5         4
```
  
### 算术
```R
> x <- IRanges(start=c(40, 80), end=c(67, 114))
> x
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        40        67        28
  [2]        80       114        35

## + 起始均向外拓展
> x + 4
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        36        71        36
  [2]        76       118        43

## - 起始均向内缩短
> x - 10
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        50        57         8
  [2]        90       104        15

## * 改变width值
## 乘上正数相当于“放大”（宽度变窄）
## 负数相当于“缩小”（宽度变长）
> x * 2
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        47        60        14
  [2]        89       105        17

> x * (-2)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        26        81        56
  [2]        62       131        70
```
  
### 截取指定边界内的ranges
```R
> y <- IRanges(start=c(4, 6, 10, 12), width=13)
> y
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        16        13
  [2]         6        18        13
  [3]        10        22        13
  [4]        12        24        13
> restrict(y, 5, 10)
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         5        10         6
  [2]         6        10         5
  [3]        10        10         1
```
  
### 旁侧序列
```R
## 默认为上游  
> x
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        40        67        28
  [2]        80       114        35
> flank(x, width = 7)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        33        39         7
  [2]        73        79         7

## start=FALSE 取下游
> flank(x, width=7, start=FALSE)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        68        74         7
  [2]       115       121         7
```
  
### 起始边界
```R
> ir
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15
  [5]        19        24         6
  [6]        34        35         2
  [7]        40        46         7
> range(ir)
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        46        46
```

### 压缩
```R
## 将重叠或相邻的range连接
> reduce(ir)
IRanges object with 3 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        29        29
  [2]        34        35         2
  [3]        40        46         7
```  
  
### 间隙
```R
## 默认返回区间之间的间隔区段
> gaps(ir)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        30        33         4
  [2]        36        39         4
```
  
### 平移
```R
## 左-右+
> shift(ir,-2)
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        -1        10        12
  [2]         6        11         6
  [3]        12        17         6
  [4]        13        27        15
  [5]        17        22         6
  [6]        32        33         2
  [7]        38        44         7
> shift(ir,2)
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         3        14        12
  [2]        10        15         6
  [3]        16        21         6
  [4]        17        31        15
  [5]        21        26         6
  [6]        36        37         2
  [7]        42        48         7
```
  
### 收窄
```R
> ir
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1        12        12
  [2]         8        13         6
  [3]        14        19         6
  [4]        15        29        15
  [5]        19        24         6
  [6]        34        35         2
  [7]        40        46         7

## start=n 参数表明相较于原起始位点，
## 新的range起始点应从第n个位置开始
> narrow(ir, start=2)
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         2        12        11
  [2]         9        13         5
  [3]        15        19         5
  [4]        16        29        14
  [5]        20        24         5
  [6]        35        35         1
  [7]        41        46         6

## end=n 则是指新的range在从起始点第n个位置处终止
> narrow(ir, end=2)
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         2         2
  [2]         8         9         2
  [3]        14        15         2
  [4]        15        16         2
  [5]        19        20         2
  [6]        34        35         2
  [7]        40        41         2

## width 限定了长度
## 注意向量循环
> narrow(ir, start=1:5, width=2)
IRanges object with 7 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         1         2         2
  [2]         9        10         2
  [3]        16        17         2
  [4]        18        19         2
  [5]        23        24         2
  [6]        34        35         2
  [7]        41        42         2
```
  
### 类集合操作  
* 并集  
```R
> a <- IRanges(start=4, end=13)
> b <- IRanges(start=12, end=17)
> a
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        13        10
> b
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        12        17         6

## range间重叠或相邻就会连接
> union(a,b)
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        17        14

## 没有相邻或重叠的话
> a <- IRanges(start=4, end=10)
> a
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        10         7
> union(a,b)
IRanges object with 2 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        10         7
  [2]        12        17         6
```
  
* 交集  
```R
> intersect(a, b)
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        12        13         2
```
  
* 差集  
```R
## a有b没有
> setdiff(a, b)
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]         4        11         8

## b有a没有
> setdiff(b, a)
IRanges object with 1 range and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        14        17         4
```
  
如果是包含多个range的对象，则实现对每个对象进行reduce操作，然后再去进行类集合操作。此外，针对等长的 `IRanges` 对象，提供了 `psetdiff()`，`pintersect()`，`punion()` 和 `pgap()` 函数进行两两配对操作  
  
## 重叠部分  
寻找重叠部分是许多分析的一个重要部分，许多分析也涉及到如何处理重叠部分，例如计算覆盖率，测序深度等  
首先来看怎么找出两组 `IRanges` 对象之间的重叠部分（结果存储在一个 `Hits` 对象中）：  
```R
> qry <- IRanges(start=c(1, 26, 19, 11, 21, 7), 
                 end=c(16, 30, 19, 15, 24, 8),
                 names=letters[1:6])
> qry
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        16        16
  b        26        30         5
  c        19        19         1
  d        11        15         5
  e        21        24         4
  f         7         8         2
> sbj <- IRanges(start=c(1, 19, 10), end=c(5, 29, 16), names=letters[24:26])
> sbj
IRanges object with 3 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  x         1         5         5
  y        19        29        11
  z        10        16         7

> hts <- findOverlaps(qry, sbj)
## 数字代表其在原IRanges对象中的索引
> hts
Hits object with 6 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         1           1
  [2]         1           3
  [3]         2           2
  [4]         3           2
  [5]         4           3
  [6]         5           2
  -------
  queryLength: 6 / subjectLength: 3

## 提取索引值，是个整数向量
> queryHits(hts)
[1] 1 1 2 3 4 5
> subjectHits(hts)
[1] 1 3 2 2 3 2
## 可以用来做索引提取名字
> names(qry)[queryHits(hts)]
[1] "a" "a" "b" "c" "d" "e"
> names(sbj)[subjectHits(hts)]
[1] "x" "z" "y" "y" "z" "y"
```
  
默认情况下我们认为只要两个范围存在重叠就是重叠区域（type="any"），可以通过 `type` 参数来控制重叠的定义  
  
另一个重要的参数是 `select`，当query的range与多个subject的range重叠时，如何处理。默认为 `select="all"`，返回所有结果，其他可选的还包括 `"first"`，`"last"` 和 `"arbitrary"`  
```R
## NA 表示没有重叠
> findOverlaps(qry, sbj, select="first")
[1]  1  2  2  3  2 NA
> findOverlaps(qry, sbj, select="last")
[1]  3  2  2  3  2 NA
> findOverlaps(qry, sbj, select="arbitrary")
[1]  1  2  2  3  2 NA
```
  
对于重叠部分的结果 `Hits` 对象，一些函数：  
```R
## 强制转换成矩阵
> as.matrix(hts)
     queryHits subjectHits
[1,]         1           1
[2,]         1           3
[3,]         2           2
[4,]         3           2
[5,]         4           3
[6,]         5           2

## 每个query重叠的命中数
## 结果为整数向量
> countQueryHits(hts)
[1] 2 1 1 1 1 0
> setNames(countQueryHits(hts), names(qry))
a b c d e f 
2 1 1 1 1 0

## 类似的，对subject
> countSubjectHits(hts)
[1] 1 3 2
> setNames(countSubjectHits(hts), names(sbj))
x y z 
1 3 2 

## 等价于 setNames(countQueryHits(hts), names(qry))
> countOverlaps(qry, sbj)
a b c d e f 
2 1 1 1 1 0

## 从query中筛选出存在重叠的range
## 等价于 qry[unique(queryHits(hts))]
> subsetByOverlaps(qry, sbj)
IRanges object with 5 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a         1        16        16
  b        26        30         5
  c        19        19         1
  d        11        15         5
  e        21        24         4
```
  
## 找出最近的range并计算距离
```R
> qry <- IRanges(start=6, end=13, name='query')
> sbj <- IRanges(start=c(2, 4, 18, 19), end=c(4, 6, 21, 24), names=1:4)
> qry
IRanges object with 1 range and 0 metadata columns:
            start       end     width
        <integer> <integer> <integer>
  query         6        13         8
> sbj
IRanges object with 4 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  1         2         4         3
  2         4         6         3
  3        18        21         4
  4        19        24         6

## sbj中与qry最近的range的索引
## 允许重叠，后面两个就不允许了
> nearest(qry, sbj)
[1] 2
## qry在前，sbj中与之最近的
> precede(qry, sbj)
[1] 3
## qry在后，sbj与之最近的
> follow(qry, sbj)
[1] 2

## 这些操作符都是向量化的，所以qry中可以有多个range
> qry2 <- IRanges(start=c(6, 7), width=3)
> nearest(qry2, sbj)
[1] 2 2

## 找出最近的range并返回距离
> qry <- IRanges(sample(seq_len(1000), 5), width=10)
> sbj <- IRanges(sample(seq_len(1000), 5), width=10)
> qry
IRanges object with 5 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]       761       770        10
  [2]       972       981        10
  [3]       895       904        10
  [4]       298       307        10
  [5]       623       632        10
> sbj
IRanges object with 5 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]       954       963        10
  [2]       842       851        10
  [3]       816       825        10
  [4]       282       291        10
  [5]       104       113        10
> distanceToNearest(qry, sbj)
Hits object with 5 hits and 1 metadata column:
      queryHits subjectHits |  distance
      <integer>   <integer> | <integer>
  [1]         1           3 |        45
  [2]         2           1 |         8
  [3]         3           2 |        43
  [4]         4           4 |         6
  [5]         5           3 |       183
  -------
  queryLength: 5 / subjectLength: 5

## 两两比较计算距离
> distance(qry, sbj)
[1] 183 120  69   6 509
```
  
## Run Length Encoding
除了基因组范围对应的碱基序列之外，还存在数值型序列数据，例如每个碱基上覆盖度的数据。为了减少这类数据占据的空间，`IRanges` 中引入了 `Run Length Encoding` 对数据进行压缩。`Run` 指一串连续的相同值    
```R
> x <- as.integer(c(4, 4, 4, 3, 3, 2, 1, 1, 1, 1, 1, 0, 0, 0,
                  0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4))
> xrle <- Rle(x)
## 即用3个4，2个3...进行存储
> xrle
integer-Rle of length 28 with 7 runs
Lengths: 3 2 1 5 7 3 7
Values : 4 3 2 1 0 1 4  

## 获取每个run的长度和值
> runLength(xrle)
[1] 3 2 1 5 7 3 7
> runValue(xrle)
[1] 4 3 2 1 0 1 4

## 转换回向量
> as.vector(xrle)
 [1] 4 4 4 3 3 2 1 1 1 1 1 0 0 0 0 0 0 0 1 1 1 4 4 4 4 4 4 4
``` 
  
`Rle` 对象支持大部分R向量操作：  
```R
> xrle + 4L
integer-Rle of length 28 with 7 runs
  Lengths: 3 2 1 5 7 3 7
  Values : 8 7 6 5 4 5 8
> xrle/2
numeric-Rle of length 28 with 7 runs
  Lengths:   3   2   1   5   7   3   7
  Values :   2 1.5   1 0.5   0 0.5   2
> xrle > 3
logical-Rle of length 28 with 3 runs
  Lengths:     3    18     7
  Values :  TRUE FALSE  TRUE
> xrle[xrle > 3]
integer-Rle of length 10 with 1 run
  Lengths: 10
  Values :  4
> sum(xrle)
[1] 56
> summary(xrle)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.00    0.75    1.00    2.00    4.00    4.00
```
  
`Rle` 对象常用于表示深度，其中 `coverage()` 函数接受一组range作为输入，并返回一个 `Rle` 对象记录每个位置的深度  
```R
> set.seed(0)
> rngs <- IRanges(start=sample(seq_len(60), 10), width=7)
> rngs
IRanges object with 10 ranges and 0 metadata columns:
           start       end     width
       <integer> <integer> <integer>
   [1]        54        60         7
   [2]        16        22         7
   [3]        22        28         7
   [4]        33        39         7
   [5]        51        57         7
   [6]        12        18         7
   [7]        49        55         7
   [8]        56        62         7
   [9]        35        41         7
  [10]        57        63         7
> rngs_cov <- coverage(rngs)
## 1-63每个位置的覆盖深度
> rngs_cov
integer-Rle of length 63 with 18 runs
  Lengths: 11  4  3  3  1  6  4  2  5  2  7  2  3  3  1  3  2  1
  Values :  0  1  2  1  2  1  0  1  2  1  0  1  2  3  4  3  2  1
```  
  
## Views
通过 `slice()` 函数，以 `Rle` 对象为输入，返回一个 `view` 对象。该对象中整合了ranges和run-length编码的向量，使得每个range均为序列部分的一个展示  
```R
## 展示序列中覆盖度不低于2X的区域
> min_cov2 <- slice(rngs_cov, lower=2)
> min_cov2
Views on a 63-length Rle subject

views:
    start end width
[1]    16  18     3 [2 2 2]
[2]    22  22     1 [2]
[3]    35  39     5 [2 2 2 2 2]
[4]    51  62    12 [2 2 2 3 3 3 4 3 3 3 2 2]
  
## 如果只对ranges该兴趣，提取这些范围
> ranges(min_cov2)
IRanges object with 4 ranges and 0 metadata columns:
          start       end     width
      <integer> <integer> <integer>
  [1]        16        18         3
  [2]        22        22         1
  [3]        35        39         5
  [4]        51        62        12
``` 
  
存在许多处理 `views` 的函数，如 `viewMeans()`，`viewMaxs()`等  
  
# GenomicRanges
`GRanges` 对象基于 `IRanges` 创建，用于存储基因组范围。其包含了额外的两个信息：`sequence name` （例如位于哪条染色体）和 `strand` (以 `run-length encoded` 方式存储)。此外，`GRanges` 对象还包含 `metadata` (可通过 `mcols()` 函数获取当中的信息)，用于记录每个基因组范围的相关信息，是一个 `DataFrame`  

```R
> library(GenomicRanges)
> gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                ranges=IRanges(start=5:8, width=10),
                strand=c("+", "-", "-", "+"))
> gr
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1      5-14      +
  [2]     chr1      6-15      -
  [3]     chr2      7-16      -
  [4]     chr3      8-17      +
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 添加metadata
> gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                ranges=IRanges(start=5:8, width=10),
                strand=c("+", "-", "-", "+"),
                gc=round(runif(4), 3))
> gr
GRanges object with 4 ranges and 1 metadata column:
      seqnames    ranges strand |        gc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1      5-14      + |     0.147
  [2]     chr1      6-15      - |     0.656
  [3]     chr2      7-16      - |     0.379
  [4]     chr3      8-17      + |      0.11
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

可以看到 `GRanges` 对象的结构：左侧包含序列名称，范围以及链的方向等基因组位置信息，右侧为元数据列  
由于基因组范围数据常对应具体某个版本的基因组，所以一般我们已经知道每条序列或染色体的长度等信息。这些信息在例如计算覆盖度等中会用到，存储于 `seqinfo` 中  

```R
## 构建时加入序列长度seqlengths信息
> seqlens <- c(chr1=152, chr2=432, chr3=903)
> gr <- GRanges(seqname=c("chr1", "chr1", "chr2", "chr3"),
                ranges=IRanges(start=5:8, width=10),
                strand=c("+", "-", "-", "+"),
                gc=round(runif(4), 3),
                seqlengths=seqlens)
> gr
GRanges object with 4 ranges and 1 metadata column:
      seqnames    ranges strand |        gc
         <Rle> <IRanges>  <Rle> | <numeric>
  [1]     chr1      5-14      + |     0.787
  [2]     chr1      6-15      - |     0.278
  [3]     chr2      7-16      - |     0.116
  [4]     chr3      8-17      + |     0.592
  -------
  seqinfo: 3 sequences from an unspecified genome

## 或者在之后添加
> seqlengths(gr) <- seqlens

## 获取这些信息
> seqinfo(gr)
Seqinfo object with 3 sequences from an unspecified genome:
  seqnames seqlengths isCircular genome
  chr1            152         NA   <NA>
  chr2            432         NA   <NA>
  chr3            903         NA   <NA>
```  

## Granges对象的常用操作
```R
## 基本属性值
> start(gr)
> end(gr)
> width(gr)

## 结果均为Rle对象
> seqnames(gr)
> strand(gr)

## 范围信息
> ranges(gr)

## 类向量操作
> length(gr)
[1] 4
> names(gr) <- letters[1:length(gr)]
> gr
GRanges object with 4 ranges and 1 metadata column:
    seqnames    ranges strand |        gc
       <Rle> <IRanges>  <Rle> | <numeric>
  a     chr1      5-14      + |     0.571
  b     chr1      6-15      - |     0.591
  c     chr2      7-16      - |         1
  d     chr3      8-17      + |     0.506
  -------
  seqinfo: 3 sequences from an unspecified genome
> gr[start(gr) > 7]
GRanges object with 1 range and 1 metadata column:
    seqnames    ranges strand |        gc
       <Rle> <IRanges>  <Rle> | <numeric>
  d     chr3      8-17      + |     0.506
  -------
  seqinfo: 3 sequences from an unspecified genome

## 针对存储元数据的DataFrame的操作
> mcols(gr)$gc
[1] 0.571 0.591 1.000 0.506
> gr$gc    ## 简写形式
[1] 0.571 0.591 1.000 0.506
```
  
## GRangesList
组合 `GRanges` 对象构成一个类似原生列表的 `GRangesList` 对象  

```R
> gr1 <- GRanges(c("chr1", "chr2"), IRanges(start=c(32, 95), width=c(24, 123)))
> gr2 <- GRanges(c("chr8", "chr2"), IRanges(start=c(27, 12), width=c(42, 34)))
> grl <- GRangesList(gr1, gr2)
> grl
GRangesList object of length 2:
[[1]] 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1     32-55      *
  [2]     chr2    95-217      *

[[2]] 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames ranges strand
  [1]     chr8  27-68      *
  [2]     chr2  12-45      *

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 列表转换成长向量
> unlist(grl)
GRanges object with 4 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1     32-55      *
  [2]     chr2    95-217      *
  [3]     chr8     27-68      *
  [4]     chr2     12-45      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 合并GRangsList
> doubled_grl <- c(grl, grl)
```
  
`GRangesList` 对象支持许多类似原生列表的操作，例如获得子列表。此外，像 `start()`、`end()` 等函数，会基于 **list-element**（即列表各项内分别作用）  

```R
> start(grl)
IntegerList of length 2
[[1]] 32 95
[[2]] 27 12
  
## 出现编码问题的话尝试
> Sys.setlocale("LC_ALL", "C")
```
  
一种常见的应用场景是我们会基于序列名称将同一条序列上的基因组范围数据合并为一组  

```R
> chrs <- c("chr3", "chr1", "chr2", "chr2", "chr3", "chr1")
> gr <- GRanges(chrs, IRanges(sample(1:100, 6, replace=TRUE),
                              width=sample(3:30, 6, replace=TRUE)))
> gr
GRanges object with 6 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr3     58-65      *
  [2]     chr1      7-24      *
  [3]     chr2      8-17      *
  [4]     chr2     66-92      *
  [5]     chr3    79-103      *
  [6]     chr1     32-46      *
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 分割  
> gr_split <- split(gr, seqnames(gr))
> gr_split
GRangesList object of length 3:
$chr3 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr3     58-65      *
  [2]     chr3    79-103      *

$chr1 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames ranges strand
  [1]     chr1   7-24      *
  [2]     chr1  32-46      *

$chr2 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames ranges strand
  [1]     chr2   8-17      *
  [2]     chr2  66-92      *

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 类似对原生列表应用apply系函数
> lapply(gr_split, function(x) order(width(x)))
$chr3
[1] 1 2

$chr1
[1] 2 1

$chr2
[1] 1 2

> sapply(gr_split, function(x) min(start(x)))
chr3 chr1 chr2 
  58    7    8 
> sapply(gr_split, length)
chr3 chr1 chr2 
   2    2    2

## reduce等函数可以直接按照list-element方式应用
> reduce(gr_split)
GRangesList object of length 3:
$chr3 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr3     58-65      *
  [2]     chr3    79-103      *

$chr1 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames ranges strand
  [1]     chr1   7-24      *
  [2]     chr1  32-46      *

$chr2 
GRanges object with 2 ranges and 0 metadata columns:
      seqnames ranges strand
  [1]     chr2   8-17      *
  [2]     chr2  66-92      *

-------
seqinfo: 3 sequences from an unspecified genome; no seqlengths

## 重组
> unsplit(gr_split, seqnames(gr))
```
  
# GenomicFeatures and rtracklayer  
`GenomicFeatures` 包提供了一系列函数处理 `TranscriptDb` 对象中包含的转录本相关数据，通过下面的例子可以看到实际中转录组数据的相关操作  

```R
## 已经构建好的转录组数据库对象
> library(TxDb.Mmusculus.UCSC.mm10.ensGene)
> txdb <- TxDb.Mmusculus.UCSC.mm10.ensGene

## 提取TxDb-like对象中的基因信息
## 返回结果为GRanges对象
## 类似函数查询 help(trabscripts)
> mm_genes <- genes(txdb)
> head(mm_genes)
GRanges object with 6 ranges and 1 metadata column:
                     seqnames              ranges strand |            gene_id
                        <Rle>           <IRanges>  <Rle> |        <character>
  ENSMUSG00000000001     chr3 108107280-108146146      - | ENSMUSG00000000001
  ENSMUSG00000000003     chrX   77837901-77853623      - | ENSMUSG00000000003
  ENSMUSG00000000028    chr16   18780447-18811987      - | ENSMUSG00000000028
  ENSMUSG00000000031     chr7 142575529-142578143      - | ENSMUSG00000000031
  ENSMUSG00000000037     chrX 161117193-161258213      + | ENSMUSG00000000037
  ENSMUSG00000000049    chr11 108343354-108414396      + | ENSMUSG00000000049
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome

## 也可以按照组合成GRangesList
## 类似函数查询 help(trabscriptsBy)
> mm_exons_by_tx <- exonsBy(txdb, by="tx")
> mm_exons_by_tx
GRangesList object of length 94647:
$1 
GRanges object with 1 range and 3 metadata columns:
      seqnames          ranges strand |   exon_id   exon_name exon_rank
         <Rle>       <IRanges>  <Rle> | <integer> <character> <integer>
  [1]     chr1 3054233-3054733      + |         1        <NA>         1

$2 
GRanges object with 1 range and 3 metadata columns:
      seqnames          ranges strand | exon_id exon_name exon_rank
  [1]     chr1 3102016-3102125      + |       2      <NA>         1

$3 
GRanges object with 2 ranges and 3 metadata columns:
      seqnames          ranges strand | exon_id exon_name exon_rank
  [1]     chr1 3466587-3466687      + |       3      <NA>         1
  [2]     chr1 3513405-3513553      + |       4      <NA>         2

...
<94644 more elements>
-------
seqinfo: 66 sequences (1 circular) from mm10 genome

## 找出与某一区域重叠的特征数据
## 类似函数查询 help(transcriptByOverlaps)
> qtl_region <- GRanges("chr8", IRanges(123260562, 123557264))
> qtl_region_expanded <- qtl_region + 10e3
> transcriptsByOverlaps(txdb, qtl_region_expanded)
GRanges object with 73 ranges and 2 metadata columns:
       seqnames              ranges strand |     tx_id            tx_name
          <Rle>           <IRanges>  <Rle> | <integer>        <character>
   [1]     chr8 119910841-124345722      + |     47374 ENSMUST00000127664
   [2]     chr8 123254195-123269745      + |     47530 ENSMUST00000001092
   [3]     chr8 123254271-123257636      + |     47531 ENSMUST00000150356
   [4]     chr8 123254284-123269743      + |     47532 ENSMUST00000156896
   [5]     chr8 123254686-123265070      + |     47533 ENSMUST00000154450
   ...      ...                 ...    ... .       ...                ...
  [69]     chr8 123559201-123559319      - |     49320 ENSMUST00000178208
  [70]     chr8 123560888-123561006      - |     49321 ENSMUST00000179143
  [71]     chr8 123562595-123562713      - |     49322 ENSMUST00000178297
  [72]     chr8 123564286-123564404      - |     49323 ENSMUST00000179019
  [73]     chr8 123565969-123566087      - |     49324 ENSMUST00000179081
  -------
  seqinfo: 66 sequences (1 circular) from mm10 genome
```
  
除了使用一些 R 包中已经构建好的转录组数据库对象，还可以利用 `GenomicFeatures` 中的函数如 `makeTxDbPackageFromUCSC()`，`makeTxDbPackageFromBiomart()` 等从 UCSC，Biomart等中自行下载数据进行构建，此外，`rtracklayer` 包提供了更加灵活的函数用于处理其他格式的范围数据，如 `GTF/GFF`, `BED` 等
