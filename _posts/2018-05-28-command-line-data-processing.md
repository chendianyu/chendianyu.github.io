---
title: Data Processing by Linux Command Line
description: Quick and lightweight data processing by Linux command line tools
categories:
 - Linux
tags:
---

> 我们有很多工具可以用于数据操作，例如 **Pandas** 或者 **Tidyverse**，但有时候我们只需要对数据有个大概了解或者进行一次性的简单操作，
这时候写一个Python，Perl或者R脚本未免显得小题大做。事实上，Linux中有不少简单高效的命令可以帮助我们完成这些任务，接下来让我们一起来学习
一下这些命令吧。

<!-- more -->

# head
```shell
head <file>
head -n 3 test.txt`       ## 等价于 `head -n3 test.txt` 等价于 `head -3 test.txt`  
head -n -<num> test.txt   ## 显示除倒数 `<num>` 行外的其他行
```

# tail
```shell 
tail <file>
tail -n +<num> test.txt   ## 显示从第 `<num>` 行开始的所有行（注意与 `head` 在这一点上的区别）
```

# tr
强大的文件清理工具，可快速实现sed的一些基本功能
```shell
tr <set1> [set2]    # 要注意的是set中的字母是单个单个对应
tr 'abc' 'ABC'      # 指将a转换成A，b转换成B，c转换成C，而不是指出现abc转成ABC
```
* `-c`: use the complement of SET1  
* `-d`: delete characters in SET1, do not translate  
* `-s`: 在set1中列举的字符，如果原文本中出现连续重复，则进行压缩只保留单个字符  
  
set1或2可以使用内置的`[:class:]`参数，从而实现更方便的指定：  
* `[:alnmu:]`: all letters and digits  
* `[:alpha:]`: all letters  
* `[:blank:]`: all horizontal whitespace, 即空格和制表符  
* `[:digit:]`: all digits  
* `[:lower:]`: all lower case letters  
* `[:punct:]`: 所有标点符号  
* `[:space:]`: all horizontal or vertical whitespace  
* `[:upper:]`: all upper case letters  
  
# wc  
```shell
wc <file>    # print newline, word, and byte counts for each file
# 对应 `-l`, `-w`, `-m` 参数  
# `-L` 则输出最长行的字符数  
```  

# cloumn
Sometimes it's difficult to see which elements belong to a particular column, 
so we use `column` to help us to visualize
* `-t` : treat data as table  
* `-s` : 指出文件的分隔符（默认为 **TAB**）  

# cut
* `-d` : 指定分隔符（默认为 **TAB**），需与 `-f` 联用指定返回的列，否则报错  
  `-d,` 、 `-d ,` 、 `d ','`均可
* `-f` : 仅与 `-d` 联用，指定返回哪几列（计数从**1**开始）  
  * `1-3` : 1-3列  
  * `1,3` : 1，3列  
  * `-3` : 等价于1-3列  
  * `3-` : 3-最后列  
  * `-` : 所有列

## Warning
`cut` 不能实现列的重排序，所以诸如 `8,1-3` 这样的形式仍会按照 `1-3,8` 返回列

# sort
* `-t`: 指定分隔符（默认为 `blank characters`， 包括 `TAB` 和 `Space`）  
* `-k<start>,<end>`: 指定用于排序的列的范围（即从`<start>`到`<end>`列）  
*eg :* `sort -k1,2 test.txt` `sort -k1,1 test.txt`  
* `-f`: 忽略大小写  
* `-n`: 按照数值顺序排序，而非字母表顺序  
* `-r`: 逆序，即降序（默认为升序）  
* `-V`: 智能化字母表排序，能理解字符串中的数字，特别适合如 `chr1, chr2, ... ,chr22` 等的排序  
**trick**: `-n`, `-r` 以及 `-V` 既可以直接单独使用，表示对排序所用的全部列均有作用，
也可以放在 `-k` 指定的列后，如 `sort -k1,1nrV test.txt` 表示针对第一列

# uniq
* `-c`: 计数  
* `-i`: 忽略大小写  
* `-d`: 仅打印重复行

# grep
```shell
grep <pattern> <files>
# (The quotes around the *pattern* aren't required, but it's safest to use quotes to avoid conflict)  
# `grep` returns any lines that match the pattern, even ones that only **partially match**
```  

* `-v` : invert  
* `-w` : match **entire words** to avoid partially match  
* `-B <num>` : 同时返回匹配行的前 `<num>` 行
  `-B<num>` 等价于 `-B <num>` 等价于 `-B'<num>'` 等价于 `-B '<num>'`  
* `-A <num>` : 同时返回匹配行的后 `<num>` 行  
* `-C <num>` : 同时返回匹配行的前、后 `<num>` 行  
* `-E` : 调用ERE  
* `-c` : 返回匹配行的行数  
* `-o` : 仅返回匹配的部分，而非整行  
* `-h` : 同时处理**多个文件**时，`grep` 默认返回**文件名：匹配行**的形式，使用 `-h` 可以避免输出文件名  

## Regular Expression
`grep` 既能接受普通的字符串作为 `<pattern>`， 也能够使用正则表达式  
包括 ：*POSIX Basic Regular Expression* (BRE) 和 *POSIX Extended Regular Expression* (ERE)

# awk
**Basic Synatx** : `awk '<pattern> { <action> }' <file>`  
`<pattern>` is an expression or regular expression pattern (regular expressions are in **slashes (/)**)  
we can chain multiple pattern-action pairs, separated by **semicolons (;)**  
  
* `-F` : 指定分隔符（默认为若干个空格，包括TAB）  
