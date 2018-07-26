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
![CD8_developmental_trajectory](/img/2018-07-26-NSCLC-T-cell-global-characterization/CD8_developmental_trajectory.png)  
  
基于表达数据的发育轨迹图也支持了之前的 TCR 分析，通路上相邻分支相较距离较远的分支有更多共同的 TCR（下图）。  
![shared_clonotypes](/img/2018-07-26-NSCLC-T-cell-global-characterization/shared_clonotypes.png)  
  
为了更好的理解发育轨迹，作者根据以前研究中确定的特征基因定义了初始性分值和细胞毒性分值，以及基于肿瘤浸润性耗竭 CD8<sup>+</sup> T 细胞中高表达的90个基因的平均表达量定义了 T 细胞耗竭分值。这90个基因与从肝细胞癌，黑色素瘤以及感染了淋巴细胞性脉络丛脑膜炎病毒的小鼠模型中鉴定出来的耗竭相关基因存在明显的重叠。基于人类肿瘤研究得到的23个基因中有14个在小鼠研究中并为发现，可能代表人类特有的耗竭标志基因。对比 Monocle 中成分与这些功能性分值之间的关系，发现成分1与 T 细胞耗竭高度相关，成分2与细胞毒性正相关，与初始性负相关（下图）。这些结果表明 NSCLC 中 CD8<sup>+</sup> T 细胞的状态是由固有的 T 细胞发育程序和肿瘤诱导的 T 细胞耗竭两个不同过程共同决定的。  
![CD8_functional_scores](/img/2018-07-26-NSCLC-T-cell-global-characterization/CD8_functional_scores.png)  
  
最近组织驻留记忆 T 细胞（tissue-resident memory T, Trm）特征被发现与肺癌中细胞毒性 T 细胞应答的程度有关。虽然 Trm 标志基因如 CD69，ITGAE (CD103)，CXCR6 和 ITGA1 在整个肿瘤驻留 T 细胞中高表达，但在本研究中的3个相关的 CD8 分支：CD8-C4-GZMK，CD8-C5-ZNF683 和 CD8-C6-LAYN 中，其表达特征存在差异。CD8-C4-GZMK 表达低水平的 ITGAE，但高表达 PDCD1；CD8-C5-ZNF683则恰好相反；CD8-C6-LAYN 同时高表达 ITGAE 和 T cell 耗竭标志基因，包括 PDCD1，CTLA4，HAVCR2，LAG3 和 TIGIT。在其他 NSCLC 肿瘤中基于 ITGAE，PD1 和 CTLA4 抗体的流式细胞术分析表明3个分支在蛋白水平上存在差异。以上结果表明 **Trm 细胞不能被看成是一个同质的群体，尤其是在肿瘤中**。  
  
CD8-C4-GZMK 和 CD8-C5-ZNF683 位于发育轨迹中心，耗竭分值相较 CD8-C6-LAYN 更低，表明可能是处于“耗竭前”状态。这些细胞在肺癌中也存在，不过 CD8-C5-ZNF683 细胞在 NSCLC 中更加富集（下图），意味着它们可能在 NSCLC 免疫中发挥作用。  
![HCC](/img/2018-07-26-NSCLC-T-cell-global-characterization/HCC.png)  
  
针对 TCGA 中肺腺癌（lung adenocarcinoma, LUAD）队列的数据分析发现高表达 CD8-C5-ZNF683 特征基因，低表达 CD8-C6-LAYN 特征基因的病人整体存活率要明显优于 CD8-C5-ZNF683/CD8-C6-LAYN 比率低的病人（下图）。类似的，对比高 CD8-C4-GZMK/CD8-C6-LAYN 和低比率的病人，也存在存活率的差异。因此，**“耗竭前” T 细胞比上耗竭 T 细胞的比例与 LUAD 预后相关**。另外，该模式在肺鳞细胞癌中并未发现，表明不同肺癌类型之间可能存在差异。  
![survival](/img/2018-07-26-NSCLC-T-cell-global-characterization/survival.png)  
  
CD4<sup>+</sup> T 细胞发育轨迹分析也呈现出分支结构（下图），可以看到 CD4-C2-ANXA1，CD4-C4-CD69，CD4-C5-EOMES 和 CD4-C6-GZMA 作为桥梁串联起 naïve CD4-C1-CCR7，effector CD4-C3-GNLY 和 exhausted CD4-C7-CXCL13 clusters。  
![CD4_developmental_trajectory](/img/2018-07-26-NSCLC-T-cell-global-characterization/CD4_developmental_trajectory.png)  
  
Monocle 成分1与 T 细胞毒性正相关，与初始性负相关，成分2则与耗竭性相关（下图），与 CD8<sup>+</sup> T 细胞类似。  
![CD4_functional_scores](/img/2018-07-26-NSCLC-T-cell-global-characterization/CD4_functional_scores.png)  
  
两个 Treg 分支中，CD4-C8-FOXP3 主要存在于血液中, 而 CD4-C9-CTLA4 则主要存在于肿瘤中。Treg 的克隆扩增是最具分支特异性的，绝大部分只在肿瘤中存在，表明**抑制性肿瘤驻留 Tregs 的局部性克隆性扩展**。CD4-C9-CTLA4 相较 CD4-C8-FOXP3 具有更多的克隆扩展细胞，且有404个基因表达量至少高2倍。这些基因与之前针对黑色素瘤，乳腺癌，结直肠癌和肺癌的研究中所发现的肿瘤驻留 Tregs 中的基因存在高度重叠性。  
针对 CD4-C9-CTLA4 分支进一步分析发现，根据与 Tregs 功能相关的基因的表达，该分支内存在相当程度的异质性。尤其是 TNFRSF9（4-1BB），一种已知的抗原特异性 Tregs 活化标志基因，呈现出显著的双峰表达分布（下图），可能代表**最近刚激活的和早就存在的两个群体**。  
![TNFRSF9](/img/2018-07-26-NSCLC-T-cell-global-characterization/TNFRSF9.png)  
  
相较 TNFRSF9<sup>-</sup> 群体，TNFRSF9<sup>+</sup> 高表达260个基因，与之前提到的404个基因重叠。这些基因中包括免疫抑制功能相关的基因如 REL 和 LAYN，进一步证明 **TNFRSF9<sup>+</sup> 细胞不仅是经过抗原刺激的，而且是肿瘤功能性 Tregs 的主要组成部分**。对 TCGA LUAD 数据集的生存分析发现260个基因的特征集与更差的预后相关（下图），而404个基因的特征集则无明显关联。因此，对于活化的肿瘤 Tregs 的临床意义，我们需要进一步的深入研究。  
![TNFRSF9_survival](/img/2018-07-26-NSCLC-T-cell-global-characterization/TNFRSF9_survival.png)  
  
将来自所有肿瘤富集 T 细胞簇的特征基因作用于 TCGA LUAD 样本，作者发现病人们具有异质的肿瘤浸润淋巴细胞特征，但总体而言可以分为两大组（如下图，9种 T 细胞分支对应的特征基因）。组1中耗竭前 CD8<sup>+</sup> T cells(CD8-C4-GZMK)，未活化 Tregs(TNFRSF9<sup>-</sup> cells of CD4-C9-CTLA4) 和活化 CD4<sup>+</sup> cells (CD4-C5-EOMES)存在富集；组2中耗竭 T 细胞(CD8-C6-LAYN 和 CD4-C7-CXCL13) 以及活化 Tregs(TNFRSF9<sup>+</sup> cells of CD4-C9-CTLA4)存在富集。对比发现组1病人相较组2具有更好的预后。因此，**T 细胞组成情况可能是一个用于人类 NSCLC 分型的重要标志**。  
![LUAD_classify](/img/2018-07-26-NSCLC-T-cell-global-characterization/LUAD_classify.png)  
  
对 TIL 实施精细地分类为当下免疫治疗药物的选择性使用提供了机会。对目前临床试验中作为肿瘤免疫治疗靶标的基因进行分析发现，属于同一类别的靶标（例如效应分子中激活或抗 Treg）在数据集的细胞中呈现出不同的表达模式（下图）。例如 PDCD1，在耗竭 CD8<sup>+</sup> T 细胞（CD8-C6-LAYN），两类 CD4<sup>+</sup> T 细胞（CD4-C5-EOMES 和 CD4-C7-CXCL13）高表达，而 CTLA4 则在抑制性肿瘤 Tregs 和耗竭 CD8<sup>+</sup> T 细胞中高表达。以上结果与最近研究发现抗-CTLA4 和抗-PD1 分别靶向不同 TIL 群体实现肿瘤退化的结果一致。LAG3 基本上只在 CD8<sup>+</sup> T 细胞中表达，而 CD27 在所有肿瘤浸润 T 细胞中均表达。因此，**不同 T 细胞群体可能可以被不同的免疫疗法调节，而这类细胞子集的组成可能能够实现高效的病人分层**。  
![target_genes_expression](/img/2018-07-26-NSCLC-T-cell-global-characterization/target_genes_expression.png)  
  
最后总结，本文献中作者对 NSCLC 患者进行高深度单细胞转录组测序分析，描绘出其 T 细胞图谱，找到了2个可能处于耗竭前状态的 CD8<sup>+</sup> T 细胞分支，能够在 LUAD 中预示更佳的预后表现。严重耗竭的 T 细胞由于其表观遗传学变化导致其对免疫检查点抑制产生抗性，因此耗竭前细胞可能可以替代成为免疫疗法的靶标。此外，由于组织间存在高度迁徙性效应 T 细胞，以及肿瘤内存在的耗竭前 T 细胞，可能使得 NSCLC 对免疫治疗呈现出积极的响应。耗竭 CD8<sup>+</sup> T 细胞，耗竭前细胞，以及活化的肿瘤 Tregs 的标志基因，都可能可以作为 LUAD 患者的临床标志物。由于不同的 T 细胞组成可能呈现不同的临床表征，因此整个 T 细胞分支的构建为病人分类提供了新的策略。
