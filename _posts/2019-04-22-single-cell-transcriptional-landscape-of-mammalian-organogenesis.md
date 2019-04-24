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

哺乳动物器官发生是一个奇妙的过程：在较短的时间内，三个胚层的细胞转化形成包含绝大部分内部和外部器官结构的胚胎。本研究中，作者通过 sci-RNA-seq3，对来自61个E9.5-13.5阶段的小鼠胚胎的约2百万个细胞进行转录组测序，获得小鼠器官发生图谱（MOCA），并利用 Monocle3，鉴定了上百种细胞类型，56个轨迹以及对应的标志基因

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
  
该实验中，共恢复出2,058,652个小鼠胚胎细胞和13,359个 HEK-293T or NIH/3T3 细胞（UMI计数≥200）。对照组孔中，来自人或小鼠的转录组表现出极高的物种一致性（即单个细胞中只含有人或小鼠的转录本），仅有3%（420/13359）出现冲突，与之前实验的表现相近（下图）。其中存在的一个缺陷是仅有约7%的细胞被测到，剩下大部分在过滤除去聚集的细胞核的过程中丢失了  
![performance_in_control](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/performance_in_control.png)  
  
每个胚胎被测到的细胞中位数为35,272个（下图左），虽然测序深度较低（平均每个细胞5000条原始序列，46%的重复率），但每个细胞最终获得的 UMI 数量的中位数为671（519个基因）（下图右）。在部分孔中将测序深度提升至3.7倍后导致复杂程度也加倍（每个细胞 UMI 中位数为1142，87%的重复率）  
![cell_gene_umi](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/cell_gene_umi.png)  
  
由于是对细胞核中的 RNA 进行测序，因此59%的 UMI 映射至内含子上，25%映射至外显子上，反映了早期的转录情况，虽然存在时间上的偏移，但也足够用于预测细胞中转录组的状态。处于晚期状态的胚胎 UMI 数量减少，可能反映了核中 mRNA 含量的减少（下图）。利用 Scrublet 检测到4.3%的细胞可能为 doublet  
![decreased_umi](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/decreased_umi.png)
  
基于对不同时间点胚胎中细胞数量的粗略估计，可以估算出不同时间点胚胎的“细胞覆盖率”分别为0.8X（E9.5，每个胚胎200,000个细胞，总共获得152,000个该阶段的细胞），0.3X（E10.5，378000/1100000），0.2X（E11.5，616000/2600000），0.08X（E12.5，475000/6000000），0.03X（E13.5，437000/13000000）  
基于比对到 *Xist* 和 Y 染色体上的reads数量，能够很容易地将胚胎分为雄性（31）和雌性（30）。将各胚胎所有细胞的转录组整合获得整体的转录组，然后 t-SNE 可视化可以看到对应5个发育阶段的分支（下图左）。此外，对这些胚胎进行拟时序排列（下图右），发现 E9.5–E10.5 和 E11.5–E12.5 存在较大的空隙，表明这两个窗口期出现了特别显著的变化。基于这些“拟混池”数据，共鉴定出12,236个不同发育阶段差异表达的基因  
![pseudobulk_result](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/pseudobulk_result.png)
  
# 细胞类型及亚型的鉴定
对2,058,652个细胞进行 Louvain 聚类以及 t-SNE 可视化，发现同一发育阶段的胚胎细胞分布类似，不同阶段的则不是。总计获得40个不同的分支，手动注释细胞类型。其中2个分支均对应 definitive erythroid lineage，因此将它们合并；另外有一个分支 doublet 率达52%，因此移除该分支，所以最终得到38个主要的分支（下图）  
![all_cell_tsne](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/all_cell_tsne.png)
  
高度特异性的标志基因能够帮助注释细胞类型。下图为每一类型挑选出一个标志基因得到的可视化结果，点的大小对应该分支细胞中有多少检测到了该基因，颜色对应平均表达水平。可以看到大部分细胞类型特意表达标志基因，即便是小分支，也能够比较容易地注释。部分标志基因在多种细胞类型中均有出现，但是在某一分支中会高度表达。对应胚胎间质和结缔组织的分支的注释比较困难，因为已知的标志基因较少  
![marker_gene](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/marker_gene.png)
  
主要分支之间，68%的基因（17789/26183）存在差异表达，其中2863个基因为细胞类型特异性标志基因（定义为在表达最多和次多的细胞类型中的表达量相差2倍以上），平均每个每种类型细胞约75个（下图）。这些基因中的大部分之前未被认为是细胞类型标志基因，比如在脊索（分支30）中，*Tox2*，*Stxbp6*，*Schip1*，*Spon1* 等之前并未被认为是标志基因，但在本研究中是分支30的标志。  
![cluster_specific_marker](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/cluster_specific_marker.png)
  
器官发生过程中，各类型细胞占比发生了明显的变化。大部分主要细胞类型呈指数增长，但有些细胞类型则是短暂出现，然后在E13.5消失（下图左和中）。例如在E9.5存在起源于卵黄囊（分支26，标志基因 *Hbb-bh1*）的 primitive erythroid lineage，但是来自胎儿肝脏（分支22，标志基因 *Hbb-bs*）的 definitive erythroid lineage 逐渐取代前者，最终成为E13.5时期唯一的红细胞系  
![cell_change](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/cell_change.png)
  
38个主要细胞类型包含细胞数量的中位数为47,073，同一大类细胞之间也存在明显的异质性，因此对这些主要类型再一次进行 Louvain 聚类，移除那些1-2个胚胎却占据绝大部分的分支，并合并高度相似的分支，最终得到655个亚型（下图）。检测细胞类型及亚型的敏感性取决于本实验大量的细胞数目，降低取样后检测到的分支数减少  
![subcluster](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/subcluster.jpg)
  
13%的亚型被认为可能是人为错误（这些分支中doublets超过10%）而剔除。剩下572个亚型，特异性标志基因（同一主细胞类型内，表达量最高和次高的亚型之间差异超过2倍）中位数为20。绝大部分亚型能够基于标志基因集与其他571个亚型区分开来，其中63%通过2个标志基因即可实现区分，95%通过4个标志基因可实现区分。下图显示用于区分亚型所需的标志基因数量  
![number_of_markers_to_distinguish](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/number_of_markers_to_distinguish.png)
  
为了与其他细胞图谱进行区分，作者将该数据集命名为小鼠器官发生细胞图谱（mouse organogenesis cell atlases, MOCA）。将MOCA的细胞亚型与MCA中130个胎儿（E14.5）细胞类型比较，将96个MCA细胞类型与58个MOCA亚型对应。MOCA中无法与MCA对应的那些亚型多来自早期阶段（如神经管）或比较罕见（如晶状体），而MCA中无法与MOCA对应上的细胞主要是组织特异性免疫或上皮细胞，可能是E13.5之后才出现的。两个数据集之间可以相互印证：MCA中的解剖学信息可用于MOCA中亚型的定位，比如MOCA中一种内分泌上皮细胞亚分支对应MCA中位于胎儿胃部的腺泡和内分泌细胞；MOCA则可用于确定MCA中细胞的胚胎起源，例如MCA中定位于胎儿肾脏的处于细胞周期的细胞对应MOCA中过渡态中胚层的一个亚型，可能是肾脏的祖细胞。类似地，MOCA中的68个亚型与BCA中的48细胞类型对应  
  
# 外胚层顶鞘的鉴定
作者对上皮细胞（分支6）（下图）所有的亚型进行了注释
![epithelial_subtypes](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/epithelial_subtypes.png)  
  
为了深入研究某一亚型，作者将目光聚焦到外胚层顶鞘（apical ectodermal ridge, AER）（亚型6.23）上，一类参与手指或足趾发育的高度特异性上皮细胞亚型。除了已知的标志基因外，本研究中还发现该亚型细胞可通过 *Fndc3a*，*Adamts3*，*Slc16a10*，*Snap91* 和 *Pou6f2* 的表达情况与其他亚型区分开来。WISH实验结果确认了 *Fgf8*（已知的标志基因），*Fndc3a*，*Adamts3* 和 *Snap91* 在腋芽最远端特异性表达，而该区域在E10.5或E11.5阶段对应AER（下图）  
![AER_marker](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/AER_marker.png)  
  
之后作者分析了AER增殖以及基因表达的动态过程。虽然在所有时间点以及几乎所有胚胎中都有测到，但是估算结果显示AER细胞数量在E10.5-E11.5之间达到峰值，与之前的报道相符。对AER细胞进行拟时间排序，得到一条简单的从早期到晚期的轨迹以及510个差异表达的基因。可以看到，像 *Fgf8*，*Fgf939* 和 *Rspo234*，活化过程后于 *Fndc3a*，*Mik67* 和 *Igf2* 等推动细胞增殖的基因的表达则是明显下降。通路层面的分析结果也显示该窗口期增殖相关程序下调  
![AER_marker_dynamics](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/AER_marker_dynamics.png)  
  
# 重构发育轨迹
之后作者探索器官发生过程中各类型细胞的发育轨迹。大部分现有的算法假设过程是连续的，但本研究是从E9.5起始的，不可避免地存在部分祖细胞或状态缺失的情况。此外，大部分算法也不允许出现细胞命运的收敛，但我们知道一些细胞类型可以通过不同的途径衍生而来。因此，作者开发了 Monocle3。Monocle3 首先通过 UMAP 将细胞映射至一个编码转录状态的低维空间，然后通过 Louvain community detection 算法将相似的细胞组合，然后将相近的组整合成“超组”，最后确定细胞发育过程中经过的路径或轨迹，分支的位置以及每个超组内的收敛性  
Monocle3 将1,524,792个高质量细胞（UMI大于400）分成12组轨迹，其中2组对应感知神经元和2组对应血细胞的分别合并，最终剩下10组（下图）其中最复杂的为神经管-脊索轨迹（轨迹3，包括了脊索，神经管，神经和胶质细胞类型祖细胞和发育中细胞）和间质轨迹（轨迹2，包括所有间质和肌肉细胞类型）。此外还有3条神经嵴轨迹，分别对应感知神经元，Schwann 前体细胞和生黑色素细胞。造血轨迹（轨迹5）包括巨核细胞，红细胞和白细胞。剩下的4条轨迹（内皮，上皮，肝和晶状体）均对应一种主要的细胞类型  
![ten_tracjectories](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/ten_tracjectories.png)  
  
38个主类型基本上都仅归属于当中某一个组轨迹中（下图，颜色对应各类型细胞分到各组轨迹中占据的比例）  
![cell_type_2_tracjectories](/img/2019-04-22-single-cell-transcriptional-landscape-of-mammalian-organogenesis/cell_type_2_tracjectories.png)  
  
与t-SNE不同，UMAP能将相关的细胞类型放置在邻近的位置。比如在较晚发育节点发现的细胞类型如抑制性神经元，与早期中枢神经系统前体细胞（辐射胶质细胞）通过神经元前体细胞构成的“桥”相连
