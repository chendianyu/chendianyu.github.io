---
title: 单细胞转录组数据轨迹预测方法的比较
description: 
categories:
 - Single Cell
tags:
 - paper
 - single cell
 - tracjectory inference
---

单细胞转录组学数据能够用于推断轨迹，从而对细胞动态变化过程进行无偏倚研究。本文中，作者对45种轨迹预测方法进行了评估，并为研究者选择合适的软件方法提出了实质性的指导意见

<!-- more -->

# Introduction
单细胞数据可以用于研究细胞动态过程，如细胞周期，分化以及细胞活化等。轨迹预测（trajectory inference, TI）方法，又称拟时间排序（pseudotime analysis）方法，能够基于表达模式的相似性，将细胞沿着一条线性，分叉，或者更为复杂的树或图状结构进行排列，从而用于鉴定新的细胞亚群，描绘分化轨迹，以及推断导致分叉的调控程序等。从2014年起，有超过70种 TI 方法发布，但还未有人对这些方法进行一个系统的比较。不同的方法具有不同的特征，包括对输入的要求，输出结果的结构，以及使用的具体算法等，使得使用者难以确定到底使用哪种方法用于自己的研究，所以我们需要从各种方法的表现，扩展性，鲁棒性以及用户使用角度入手，对其进行全面的评估  
  
# 方法筛选以及评价标准
本研究中共囊括了45种 TI 方法。45种方法的主要特征以及总体的评估结果如下图所示  
![method_characterization](/img/2019-04-26-trajectory-inference-method-comparison/method_characterization.png)  
  
部分方法对输入数据有要求，有指定一个起始细胞的，还有则需已知细胞分组。如果某种方法必须，则提供其所需的先验信息  
不同方法的输出结果也不同，但是能够将它们分为7大组，作者为每一组编写了通用的转换器（下图），实现对结果的封装，用于比较  
![wrapper_function](/img/2019-04-26-trajectory-inference-method-comparison/wrapper_function.png)  
  
为了能够直观地比较不同方法的输出结果，作者开发了一套通用概率模型用于表示轨迹（下图）。在该模型中，整个轨迹的拓扑结构通过一个由“里程碑”构成的网络来表示，细胞则置于相连的点构成的空间中  
![common_trajectory_model](/img/2019-04-26-trajectory-inference-method-comparison/common_trajectory_model.png)  
  
不同方法间最大的差异在于是否固定了结果的拓扑结构，如果没有的话，那么它能够检测出哪些拓扑结构。作者共定义了7种不同结构（下图）。绝大部分方法只能预测线性或其他较简单的结构，只有极少数方法能够预测环形或者分离式的拓扑结构  
![trajectory_type](/img/2019-04-26-trajectory-inference-method-comparison/trajectory_type.png)  
  
作者从以下4大核心层面入手评估每种方法：  
1. **预测的准确性（accuracy）**，总计110个真实数据集和229个合成数据集  
2. **方法的扩展性（scalability）**，细胞和特征（如基因）数量增加后的表现  
3. **预测的稳定性（stability）**，数据集抽样后预测的结果情况  
4. **工具的可用性（usability）**，包括软件的实现，说明文档等
  
具体的细则看下一部分结果  
  
# 评价结果
所有方法具体的评价结果如下图所示  
![detailed_evalution_result](/img/2019-04-26-trajectory-inference-method-comparison/detailed_evalution_result.png)  
  
## 准确性  
评估准确性使用的数据集包括两种类型：  
* 真实数据集，涉及不同的 scRNA-seq 技术，物种，动态过程以及拓扑结构等。如果一个真实数据集参考轨迹不是通过表达数据本身得到的，而是基于细胞分选（cell sorting）或细胞混合（cell mixing）得到，则归为 gold standard，其他的均归为 silver standard  
* 合成数据集，使用了多种数据模拟器，其中包括一种利用基因调控热力学模型的基因调控网络模拟器。对于每个模拟，以真实数据集作为参考，并匹配其维度，差异表达基因的数量，droupout 率以及其他统计特性  
  
选择了4项指标（如下图所示），包括：  
* 拓扑结构（topology）（Hamming–Ipsen–Mikhailov, HIM）  
* 细胞分配到分支的质量（F1<sub>branches</sub>）  
* 细胞的位置（cor<sub>dist</sub>）  
* 差异表达特征沿着轨迹分布的准确性（wcor<sub>features</sub>）  
  
![accuracy_metric](/img/2019-04-26-trajectory-inference-method-comparison/accuracy_metric.png)  
  
结果显示**各方法在不同数据集上的表现差异较大**，意味着不存在适用于所有数据集的“万全之策”。不同来源的数据集的整体分值与具有金标准的真实数据的整体分值之间存在中度到较高的相关性，表明**金标准轨迹的准确性以及合成数据集之间的相关性**。**不同指标得到的结果常出现不一致的情况**，例如 Monocle 和 PAGA Tree 在拓扑结构打分上得分最高，而 Slingshot 则是在对细胞进行排序以及将其置于正确的分支方面表现更佳  
