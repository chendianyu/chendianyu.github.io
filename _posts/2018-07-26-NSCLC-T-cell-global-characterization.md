---
title: scRNA-seq 揭示非小细胞肺癌中 T 细胞整体特征
description: 
categories:
 - Lung Cancer
tags:
 - paper
 - scRNA-seq
 - Lung Cancer
 - T cell
---

肿瘤免疫疗法在非小细胞肺癌（non-samll-cell lung cancer, NSCLC）中表现出持续性的效果，但其效率部分依赖于肿瘤浸润淋巴细胞（tumor infiltrating lymphocytes, TIL）的数量和组成。为了描绘 TIL 组成，谱系和功能状态的基准图谱，作者对从14位未经治疗的 NSCLC 患者中提取的12346个 T 细胞进行深度单细胞转录组测序。

<!-- more -->

NSCLC 是肿瘤相关致死率中的主要原因，占了所有肺癌病例的约85%。NSCLC 肿瘤通常具有极多的基因组变异，高免疫负荷，因此普遍对免疫检查点抑制疗法反应更好。为了揭示 NSCLC 中肿瘤浸润性 T 细胞的复杂性，作者从14位未接受治疗的病人的肿瘤，相邻正常组织以及外周血中分离出 T 细胞。这14位病人中有11位为腺癌（adenocarcinomas），3位为鳞癌（squamous cell carcinomas）。确认 NSCLC 中存在 T 细胞浸润后，作者利用流式细胞术，基于 CD3/CD4/CD8/CD25 细胞表面标记，将各种 T 细胞亚型进行分选，然后进行深度单细胞转录组测序。最终，作者共计获得12346个细胞，平均测序深度为每个细胞104万对单一匹配的 read，足够对检测地表达的细胞因子和转录因子。  
对表达数据进行 t-SNE 降维分析后，发现 T 细胞基于其组织来源和亚型实现聚类。为了进一步探索当中 T 细胞的异质性，作者基于 t-SNE+densityClust 进行无监督式聚类，确认存在7个 CD8 分支和9个 CD4 分支（如下图）。  
![T_cell_cluster](/img/2018-07-26-NSCLC-T-cell-global-characterization/T_cell_cluster.png)
  
根据特征基因表达情况，对比已知各类型细胞的功能性标志，可以确定各类细胞分支的身份（如下图）。  
![cell_marker](/img/2018-07-26-NSCLC-T-cell-global-characterization/cell_marker.png)  
  
如下图所示，不同的细胞分支在各组织中的分布存在差异。例如，CD8-C1-LEF1 和 CD4-C1-CCR7 均代表初始 T 细胞，在血液中大量存在；CD8-C6-LAYN（CD8<sup>+</sup>耗竭 T 细胞）和 CD4-C9-CTLA4（CD4<sup>+</sup> Tregs）则主要在肿瘤中聚集。两个效应 T 细胞分支，CD8-C3-CX3CR1 和 CD4-C3-GNLY 则主要血液和（或）正常组织中聚集。这两个分支高表达趋化因子受体 CX3CR1 以及与细胞毒性相关的基因如 PRF1，GZMA 和 GZMB，低表达 T 细胞耗竭标志基因如 PDCD1，CTLA4 和 HAVCR2。针对另外3个 NSCLC 患者进行流式细胞术分析发现血液中 CX3CR1<sup>+</sup>CD8<sup>+</sup> T 细胞要比肿瘤中所占的比例更高，进一步验证了之前的发现。此外，在 TCGA 数据中，肺癌内 CX3CR1 的表达较对应正常组织或者 GTEx 中的无关正常肺要少。  
![tissue_distribution](/img/2018-07-26-NSCLC-T-cell-global-characterization/tissue_distribution.png)  
  
为了探索各个 T 细胞间的克隆关系，作者获取了16个分支共8038个细胞的全长 TCR α 和 β 链，发现其中5015个细胞具有独特的 TCRs，3023个细胞的 TCRs 则存在重复。这些克隆扩增在血液和实体组织中均存在，克隆大小从2到75不等。根据相同克隆是在单个还是多个组织中存在，可以将这些克隆细胞分为组织内和组织间两大类。如下图所示，可以看到大部分细胞属于组织内克隆，但也存在一定比例的组织间克隆。  
![clonal_distribution](/img/2018-07-26-NSCLC-T-cell-global-characterization/clonal_distribution.png)  
（对角线存在富集，为组织内克隆；但组织间区域也有一些红点，表明存在组织间克隆）  
  
不同分支细胞克隆扩展程度不同（下图）。  
![clonal_each_cluster](/img/2018-07-26-NSCLC-T-cell-global-characterization/clonal_each_cluster.png)  
  
上图中我们可以看到，初始 T 细胞分支 CD8-C1-LEF1 和 CD4-C1-CCR7 中克隆扩增几乎可以忽略不计，而效应 T 细胞分支 CD8-C3-CX3CR1 和 CD4-C3-GNLY 中，不仅克隆细胞占比最高，而且相当一部分为组织间克隆（下图），表明这些 T 细胞是由相同的祖先迁移至血液和实体组织区域。  
![inter_clonal](/img/2018-07-26-NSCLC-T-cell-global-characterization/inter_clonal.png)  
  
此外，参与“黏着斑”，“细胞黏附分子（CAMs）”，“白细胞跨内皮迁移”，“肌动蛋白细胞骨架调节”等通路的基因在 CD8-C3-CX3CR1 和 CD4-C3-GNLY 分支中高表达（下图），表明这些细胞可能存在迁移特征。另外，CD4-C3-GNLY 中高表达 TBX21 基因，说明可能存在 Type 1 T helper（Th1）细胞的富集。  
![migratory_feature](/img/2018-07-26-NSCLC-T-cell-global-characterization/migratory_feature.png)  
  
作者发现9%（158/1669）的 CD8<sup>+</sup> T 细胞克隆存在于至少两个分支中（下图，同一颜色表示不同分支中相同的克隆），表明不同分支的 CD8<sup>+</sup> T 细胞并非完全独立，而是可能存在状态转换。  
![identical_T_clones](/img/2018-07-26-NSCLC-T-cell-global-characterization/identical_T_clones.png)  
  
为了探索这种状态转换，作者利用 Monocle 构建了6个 CD8 分支（除去了 CD8-C7-SLC4A10 的粘膜相关不变 T 细胞，mucosal associated invariant T，MAIT，因为它们的 TCR 太过不同）的发育轨迹图，如下图所示，可以看到整个轨迹图存在一个分支结构，CD8-C6-LAYN 耗竭分支与 CD8-C1-LEF1 初始分支和 CD8-C3-CX3CR1 效应分支位于相反的位置，而 CD8-C2-CD28，CD8-C4-GZMK 和 CD8-C5-ZNF683 分支位于中间，表明它们处于中间功能状态。  
![developmental_trajectory](/img/2018-07-26-NSCLC-T-cell-global-characterization/developmental_trajectory.png)  
  
基于表达数据的发育轨迹图也支持了之前的 TCR 分析，通路上相邻分支相较距离较远的分支有更多共同的 TCR（下图）。  
![shared_clonotypes](/img/2018-07-26-NSCLC-T-cell-global-characterization/shared_clonotypes.png)  
  
为了更好的理解发育轨迹，作者根据以前研究中确定的特征基因定义了初始性分值和细胞毒性分值，以及基于肿瘤浸润性耗竭 CD8<sup>+</sup> T 细胞中高表达的90个基因的平均表达量定义了 T 细胞耗竭分值。这90个基因与从肝细胞癌，黑色素瘤以及感染了淋巴细胞性脉络丛脑膜炎病毒的小鼠模型中鉴定出来的耗竭相关基因存在明显的重叠。基于人类肿瘤研究得到的23个基因中有14个在小鼠研究中并为发现，可能代表人类特有的耗竭标志基因。对比 Monocle 中成分与这些功能性分值之间的关系，发现成分1与 T 细胞耗竭高度相关，成分2与细胞毒性正相关，与初始性负相关（下图）。这些结果表明 NSCLC 中 CD8<sup>+</sup> T 细胞的状态是由固有的 T 细胞发育程序和肿瘤诱导的 T 细胞耗竭两个不同过程决定的。  
![functional_scores](/img/2018-07-26-NSCLC-T-cell-global-characterization/functional_scores.png)  
  
