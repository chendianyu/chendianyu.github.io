---
title: 比值比（odds ratio）与相对风险（relative risk）的差异
description: odds ratio 与 relative risk 的定义和区别
categories:
 - Statistic
tags:
 - Statistic
 - odds ratio
 - relative risk ratio
---  
  
首先来看几个概念：  
* 概率（probability），例如1/4的老虎患病，那么当我们从老虎中随机选择一只时，有1/4的概率挑选到一只患病的老虎。所以概率指的是某一阳性时间发生的次数除以总的事件的次数  
* 比值（赔率）（odds），则是阳性事件发生次数除以非阳性事件发生的次数，**P/(1-P)**，以前面的例子为例，比值就是3赔1  
* 相对风险（relative risk），是指一个概率除以另一个概率，例如老虎患病的概率除以熊患病的概率  
* 比值比（odds ratio），则是指一个比值除以另一个比值，例如老虎的比值除以熊的比值  

下面的这张示意图能够帮助我们理清这几个概念：  
![diagram](/img/2018-09-05-odds-ratio-and-risk-ratio/diagram.svg)  
  
看起来其实相对概率会更加直观，比值比反而并不容易理解，而且事实上，**在科学论文中，我们常说的 odds ratio 实际是 relative risk**，所以往往是误用。  
当概率非常趋近0时，比值比跟相对风险很接近，比如在疾病研究中，发病率比较低时，如癌症；离0越远，它们之间的差异越大，混淆两者带来的危害更大，比如一些常见的疾病，如高血压。  
odds ratio 最大的益处在于它们是对称的，例如上图中，当谈到患病时，老虎和熊的比值比是 (1:3):(1:9)=3，而当谈到未患病时，(3:1):(9:1)=1/3；但是如果用相对风险，则分别为2.5和5/6。不过，虽然如此，更多时候还是建议用相对风险，相较对称，使别人更容易理解更重要。  
  
# REF  
1. http://freerangestats.info/blog/2018/08/17/risk-ratios
