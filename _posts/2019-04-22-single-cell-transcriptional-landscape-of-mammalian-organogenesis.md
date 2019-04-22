---
title: 小鼠器官发生的单细胞转录组图谱
description: 
categories:
 - Single Cell
tags:
 - paper
 - single cell
 - monocle
---

哺乳动物器官发生是一个奇妙的过程：在较短的时间内，三个胚层的细胞转化形成包含绝大部分内部和外部器官结构的胚胎。本研究中，作者通过 sci-RNA-seq3，对来自61个胚胎的约2百万个E9.5-13.5的小鼠细胞进行转录组测序，获得小鼠器官发生图谱（MOCA），并利用 Monocle3，鉴定了上百种细胞类型，56个轨迹以及对应的标志基因

<!-- more -->

# Introduction
小鼠生长很快，从受精到出生仅需21天。从胚胎日起第4天（E4）发生囊胚植入，E6.5-E7.5原肠胚和胚层形成。在早期体节阶段，原肠胚过渡到早期器官发生，形成神经板和心管（E8.0-E8.5）。之后（E9.5-E13.5），胚胎从数万个细胞扩大到超过1000万个细胞，并发育出几乎所有主要的器官系统，因此这4天时间被广泛研究。事实上，绝大部分与重大发育缺陷相关的基因可以在该窗口期进行研究  
  
# 单细胞转录组测序
为提升通量，作者对 sci-RNA-seq 实验技术进行了摸索和改进。sci-RNA-seq3相较之前的方法主要的提升有：  
1. 不用酶处理直接从新鲜组织提取细胞核，然后进行固定和保存  
2. 在第三级的索引，不使用 Tn5 tagmentaion，而是使用发夹结扎技术  
3. 对酶反应进行优化  
4. 用稀释法代替流式分选，并通过声波降解和过滤来减少细胞聚集  

作者共收集了61个E9.5，E10.5，E11.5，E12.5或E13.5阶段的 C57BL/6 小鼠胚胎，液氮中快速冰冻，然后从胚胎中分离出细胞核，分别置于4个96孔板中。作为对照，在其中2个孔中加入人类 HEK-293T 和小鼠 NIH/3T3 的细胞核混合物。单细胞转录组文库通过单个 Illumina NovaSeq run 进行测序，最终得到110亿条reads（下图）  
![sci_rna_seq3_workflow](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/sci_rna_seq3_workflow.png)  
  
该实验中，共恢复出2058652个小鼠胚胎细胞和13359个 HEK-293T or NIH/3T3 细胞（UMI计数≥200）。对照组孔中，来自人或小鼠的转录组表现出极高的物种一致性（即单个细胞中只含有人或小鼠的转录本），仅有3%（420/13359）出现冲突，与之前实验的表现相近（下图）。其中存在的一个缺陷是仅有约7%的细胞被测到，剩下大部分在过滤除去聚集的细胞核时丢失了  
![performance_in_control](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/performance_in_control.png)  
  
每个胚胎被测到的细胞中位数为35272个（下图左），虽然测序深度较低（平均每个细胞5000条原始序列，46%的重复率），但每个细胞最终获得的 UMI 数量的中位数为671（519个基因）（下图右）  
![cell_gene_umi](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/cell_gene_umi.png)
