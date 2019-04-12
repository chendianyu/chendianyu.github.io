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
  
## GenomicRanges
`GRanges` objects contain the two other pieces of information necessary to specify a genomic location: `sequence name` (e.g., which chromosome) and `strand` (stored as `run-length encoded` in Granges). GRanges objects also have `metadata` columns (using `mcols()` function to access), which are the data linked to each genomic range. All metadata attached to a GRanges object are stored in a `DataFrame`  
So the structure of GRanges objects: genomic location specified by sequence name, range, and strand (on the left of the dividing bar), and metadata columns (on the right). Each row of metadata corresponds to a range on the same row  
Because GRanges (and genomic range data in general) is always with respect to a particular genome version, we usually know beforehand the information of each sequence/chromosome, including length, accession number and etc. These data are stored in `seqinfo`  
```R
GRanges with 4 ranges and 1 metadata column:
seqnames ranges strand | gc
<Rle> <IRanges> <Rle> | <numeric>
[1] chr1 [5, 14] + | 0.897
[2] chr1 [6, 15] - | 0.266
[3] chr2 [7, 16] - | 0.372
[4] chr3 [8, 17] + | 0.573
---
seqlengths:
chr1 chr2 chr3
152 432 903
```  
  
## GRangesList
GRanges objects' own version of list, behave like base `list`  
And we will see `RleList` and `IntegerList` corresponding to `GRangesList`  
  
## GenomicFeatures  
`GenomicFeatures` is a Bioconductor package for creating and working with transcript-based annotation, it provides methods for creating and working with `TranscriptDb` objects  
You can create your own `TranscriptDb` by function `makeTxDbPackageFromUCSC()`, `makeTxDbPackageFromBiomart()` and etc. in `GenomicFeatures`. And package `rtracklayer` provides more flexible function to deal with format such as `GTF/GFF`, `BED` and etc.
