---
title: Handbook of biological statistics notebook
description: 生物统计学习笔记
categories:
 - Statistics
tags:
 - statistics
 - notebook
---

当预计零假设很可能为真时，统计显著的结果很可能是一个假阳性，我们需要选择一个更低的 P 值作为阈值，或者改用贝叶斯统计学。例如，假设实际有1000个处理，当中一个确实有用，造成差异，另999个是无用的；如果选择0.05作为 P 值时，有用的1确实是真阳性，但对于另外999个无用的，会有5%，即约50个也能给出 P<0.05 这样的值（根据 P 值定义，零假设为真但错判成为假），这样的话总计约51个阳性结果中绝大部分是假阳性，显然不太好。所以对于这类情况，我们往往需要降低 P 值。  
  
当只有一个名义变量时，我们可以进行 Exact test of goodness-of-fit；如果观察值数量较大，则用 G–test 或 chi-square test 替代；如果有多个名义变量，则用 multinomial tests of goodness-of-fit。
