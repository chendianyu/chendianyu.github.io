---
title: Perl 学习笔记
description: Perl 学习笔记
categories:
 - Perl
tags:
 - Perl
 - notebook
---

# Introduction
`#` 代表注释  
Perl 不需要进行变量声明  
大部分 Perl 语句用 `;` 结尾，但只是用来分隔不同的语句，而不是断行  
Perl 解释器能够一次性完成编译和运行两个动作：首先载入整个源程序，将其转换成内部使用的 `bytecodes`，这是一种 Perl 在内部用来表示程序语法树的数据结构，然后交给 Perl bytecode 引擎运行  
  
# Scalar Data
标量是 Perl 中最简单的数据类型，包括数字（number）和字符串（string of characters）（Perl 中字符串就是一个独立的标量，但有必要时我们还是能够访问其中的每个字符），两者在大多数情况下可以相互转换  
不管是整数还是浮点数，在 Perl 内部都是用双精度浮点数保存的，所以像除法运算也是按浮点数操作的  
`literal`（直接量）是指某个数字在 Perl 源代码中的写法，并非运算结果，也不是 I/O 操作的结果，就是直接键入程序源代码的写法  
Perl 允许在整数直接量中插入下划线进行分隔，便于阅读，如 11_222_333_444  
