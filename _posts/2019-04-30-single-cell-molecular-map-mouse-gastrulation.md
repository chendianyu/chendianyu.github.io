---
title: 小鼠原肠胚形成及早期器官发生的单细胞转录组图谱
description: 原肠胚形成是胚胎多能细胞分化形成谱系特异性前体细胞过程中的一个重要阶段。本研究中作者对小鼠E6.5-E9.5阶段的116,312细胞进行单细胞转录组测序，从而构建原肠胚形成过程的分子图谱
categories:
 - Single Cell
tags:
 - paper
 - single cell
 - lineage tracing
---

# Introduction
小鼠E6.5-E8.5覆盖了原肠胚形成和早期器官发生的关键阶段。该窗口期，上胚层（epiblast）细胞分化成形成所有主要器官的外胚层（ectodermal）、中胚层（mesodermal）和内胚层（endodermal）细胞。为了研究该阶段细胞动态过程，作者在E6.5-E8.5之间每隔6小时进行取样，共收集411个小鼠胚胎进行 scRNA-seq。整个数据集囊括了 TS9，TS10，TS11 和 TS12（Theiler stages，Theiler等人基于一组形态学标准划分的小鼠发育阶段），分别富集于 pre-streak to early streak，mid-streak to latestreak，neural plate 和 headfold to somitogenesis 阶段  
  
# 早期胚胎发生的单细胞图谱  
严格质控后最终得到116,312个细胞，每个细胞检测到的基因中位数为3436。聚类得到**37个主要的细胞群体**（下图左），**不同类型的细胞的出现与采样的时间点存在关联**。可以看到，多能外胚层细胞的比例随着时间减少；中胚层和最终内胚层细胞最早在E6.75时期即出现；从E7.5开始，随着器官发生过程开启外胚层细胞出现，并伴随各胚层细胞类型的多样化（下图右）  
![cluster](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/cluster.png)  
  
各分支之间的转录相似性与之前的研究结果相符。外胚层与神经外胚层和原条相似，而原条则与中胚层和内胚层相关，对应三大胚层的分化。器官发生阶段（E8.25–E8.5），神经和中胚层通过神经中胚层祖细胞群相连，后者曾被报道能够生成躯干中胚层和脊索神经组织  
  
# 内胚层发育
之前的谱系追踪结果显示胚胎外和胚胎内的内胚层细胞能够插入形成单个组织，显示了胚胎细胞的可塑性，本研究中对原肠胚形成阶段的胚胎的胚外结构进行了取样，所以可在分子层面探索原条衍生的最终内胚层细胞向着内脏内胚层收敛的过程。为此，对内脏内胚层、原条前部、终末内胚层和肠道细胞类型（共5015个细胞）进行分析，结果显示肠道内胚层确实来自内脏和终末内胚层（下图左）（分别通过 *Ttr* 和 *Mixl1* 的表达鉴定）。细胞采集时间点的分析也支持这两种细胞系的转录收敛性（下图右）  
![endoderm_cell_subset](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/endoderm_cell_subset.png)  
  
为了确定成熟中肠道的转录多样性，仅选取E8.25和E8.5的细胞进行分析，得到7个对应肠管不同细胞群的分支（下图左），包括咽部内胚层（表达 *Nkx2-5*），前肠（foregut, *Pyy*），中肠（midgut, *Nepn*）和后肠（hindgut, *Cdx2*）。其中前肠分为两个分支1和2，很可能分别对应肝（相关基因 *Hhex*, *Sfrp5*, *Ttr*）和肺（*Ripply3*, *Irx1*）前体细胞（下图右）；后肠同样可以分为两个分支，其中后肠1分支显著高表达 X 染色体基因 *Trap1a* 和 *Rhox5*  
![maturing_gut](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/maturing_gut.png)  
  
考虑到肠管的空间复杂性，利用 DPT 软件对这些分支进行拟空间排序，重现了其从前端到后端的分布情况（下图）  
![pseudospace](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/pseudospace.png)  
  
