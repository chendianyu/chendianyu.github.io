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

# Highlights
* 绘制了小鼠原肠胚形成以及早期器官发生的单细胞转录组图谱，得到37个主要细胞类群及其分布情况  
* 探索了胚外和胚内内胚层的发育过程，鉴定出内脏内胚层-后肠1和定形内胚层-后肠2轨迹，以及如何共同构成肠道  
* 对红细胞的 primitive 和 yolk-sac-definitive wave 进一步探索，推测出一条原始红细胞发育轨迹，不同区域的内皮细胞具有转录组异质性。此外在 primitive wave 中鉴定出2个罕见的细胞群，注释为巨核和骨髓细胞  
* 构建嵌合型小鼠胚胎用于基因突变的研究，明确 TAL1 对血液细胞系的重要作用，其缺失导致 primitive 和 second haematopoietic wave 均中断，细胞被封锁在 second haematopoietic wave 上类似于EC3亚分支的转录组状态，且由于无法继续向着血液表型进发，这些细胞转而开始激活其他的中胚层程序  
  
# 基本概念  
定形内胚层在原肠胚形成过程中出现并取代胚外内脏内胚层，参与肠管形态发生过程，并构成相关的内脏器官。定形内胚层插入到原始内胚层，使得原始内胚层细胞扩散。之前研究认为原始内胚层会被定形内胚层完全取代，但后来发现即便是到E8.75阶段，在胚胎肠道中仍存在部分原始内胚层来源的细胞，尤其是后肠部分，占到30%  
  
# Introduction
小鼠E6.5-E8.5覆盖了原肠胚形成和早期器官发生的关键阶段。该窗口期，上胚层（epiblast）细胞分化成形成所有主要器官的外胚层（ectodermal）、中胚层（mesodermal）和内胚层（endodermal）细胞。为了研究该阶段细胞动态过程，作者在E6.5-E8.5之间每隔6小时进行取样，共收集411个小鼠胚胎进行 scRNA-seq。整个数据集囊括了 TS9，TS10，TS11 和 TS12（Theiler stages，Theiler等人基于一组形态学标准划分的小鼠发育阶段），分别富集于 pre-streak to early streak，mid-streak to latestreak，neural plate 和 headfold to somitogenesis 阶段  
  
# 早期胚胎发生的单细胞图谱  
严格质控后最终得到116,312个细胞，每个细胞检测到的基因中位数为3436。聚类得到**37个主要的细胞群体**（下图左），**不同类型的细胞的出现与采样的时间点存在关联**。可以看到，多能外胚层细胞的比例随着时间减少；中胚层和定形内胚层细胞最早在E6.75时期即出现；从E7.5开始，随着器官发生过程开启，外胚层细胞出现，并伴随各胚层细胞类型的多样化（下图右）  
![cluster](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/cluster.png)  
  
各分支之间的转录相似性与之前的研究结果相符。外胚层与神经外胚层和原条相似，而原条则与中胚层和内胚层相关，对应三大胚层的分化。器官发生阶段（E8.25–E8.5），神经和中胚层通过神经中胚层祖细胞群相连，后者曾被报道能够生成躯干中胚层和脊索神经组织  
  
# 内胚层发育
之前的谱系追踪研究显示，胚胎外和胚胎内的内胚层细胞能够交叉共同形成单个组织，显示了胚胎细胞的可塑性，本研究中对原肠胚形成阶段的胚胎的胚外结构也进行了取样，所以可用于在分子层面探索原条衍生的**定形内胚层细胞与内脏内胚层细胞汇聚的过程**。为此，对内脏内胚层、原条前部、定形内胚层以及肠道细胞（共5015个）进行分析，结果显示肠道内胚层确实来自内脏和定形内胚层（下图左）（分别通过 *Ttr* 和 *Mixl1* 的表达鉴定）。检查细胞的采集时间点也支持这两种细胞系在发育过程中出现汇聚（下图右）  
![endoderm_cell_subset](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/endoderm_cell_subset.png)  
  
为了探索成熟中肠道的转录多样性，仅选取E8.25和E8.5的细胞进行分析，得到**7个对应肠管不同细胞群的分支**（下图左），包括咽部内胚层（pharyngeal endoderm，表达 *Nkx2-5*），前肠（foregut, *Pyy*），中肠（midgut, *Nepn*）和后肠（hindgut, *Cdx2*）。其中前肠分为两个分支1和2，很可能分别对应肝（相关基因 *Hhex*, *Sfrp5*, *Ttr*）和肺（*Ripply3*, *Irx1*）前体细胞（下图右）；后肠同样可以分为两个分支，其中后肠1分支显著高表达 X 染色体基因 *Trap1a* 和 *Rhox5*  
![maturing_gut](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/maturing_gut.png)  
  
考虑到肠管的空间复杂性，利用 DPT 对这些分支进行拟空间排序，重现了其从前端到后端的分布情况（下图）  
![pseudospace](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/pseudospace.png)  
  
利用 **transport maps** 方法来推断细胞是如何沿着连续采集时间点转变的，从而确定E8.25和E8.5阶段的细胞是从E7.0的内脏内胚层还是定形内胚层衍生而来。为了解释具有永久性胚外命运的细胞，将胚外（extra-embryonic, ExE）内胚层细胞也加入到分析中。结果显示**胚外内胚层和后肠1分支中细胞主要由内脏内胚层衍生而来，其余分支则主要为定形内胚层起源**（下图）。后肠1分支特异基因 *Trap1a* 和 *Rhox5* 在胚外内胚层和胚外外胚层均表达，与该分支起源于胚外的事实一致。后肠1和2分支虽共享核心后肠特征，但前者还保留有胚外来源的转录组印记。利用 Ttr–YFP（*Ttr*，上文用于确定来自内脏内胚层的标志基因） 进行谱系追踪，结果同样显示 Ttr–YFP-traced 细胞在E8.5胚胎的后段富集  
![predicted_origin](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/predicted_origin.png)  
  
接下来，作者先推断了哪些细胞属于**内脏内胚层-后肠1**的发育轨迹，哪些又属于**定形内胚层-后肠2**的发育轨迹（下图，transport maps 方法得到，此外还用 PAGA 方法也得到类似的结果）  
![trajectory_for_hindgut_formation](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/trajectory_for_hindgut_formation.png)  
  
作者利用 DPT 将细胞排序，并将基因按照其表达沿着轨迹的动态变化过程进行聚类。基因沿着两条轨迹变化的图式都分成两块区域，其中第一区域对应在E7.5完成内胚层插入之前的基因表达情况，第二区域则对应后续时间点的基因表达（原文献 Extended Data Fig. 5c）。结果显示，内脏内胚层-后肠1轨迹中，内脏内胚层基因在第一区域上调，之后随着细胞转向肠道命运，急剧减少（下图，左边粉红色块对应内胚层插入，右边蓝色块对应肠道成熟区域）。这一结果表明有**一组内脏内胚层细胞，在定形内胚层插入之前，正进行内脏成熟过程**  
![VE_gene_dynamics](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/VE_gene_dynamics.png)  
  
在两条轨迹中，有一组基因在插入过程中均上调表达，包括参与上皮重塑的基因如 *Pcna*，*Epcam* 和 *Vim*，与该阶段预计出现的上皮排列相契合。后续肠道成熟以及形态发生过程，共同上调的基因则富集于转录因子，当中的66%为呈现出顺序激活的同源结构域蛋白，表明在后肠特化过程中的时间共线性。此外还找到了内脏内胚层-后肠1早期特异性表达的基因，包括 *Hes1*，*Pou5f1* 和 *Sox4*，可作为后续研究的候选者  
  
# 血内皮细胞系的起源
卵黄囊中红细胞可分别通过两种方式形成。一是 （primitive）从E7.5开始，产生带核的红细胞，出生后不久便消失，二是（yolk-sac-definitive）在E8.25以血管内皮出现红系骨髓祖细胞（erythro-myeloid progenitors, EMP）起始，后续迁移至胎儿肝脏形成终末红细胞。虽然已经知道了原始和卵黄囊终末造血过程在表型以及分子层面的一些关键性差异，但是对其各自在体内的祖细胞仍知之甚少  
为深入研究这些过程，将红系细胞，血内皮细胞，血细胞祖细胞，内皮细胞以及混合中胚层细胞（共15875个细胞）分离出来重新进行聚类，结果如下图所示（BP, blood progenitor; EC, endothelial cell; Ery, erythrocyte; Haem, haemato-endothelial progenitor; Mes, mesodermal cell; Mk, megakaryocyte, 巨核细胞; My, myeloid cell），其中左图内方框是对骨髓，巨核细胞以及血液内皮细胞部分的放大，右边两图则是表示亚分支之间关系的简化图  
可以看到推测出一条向着原始红细胞延伸的轨迹，穿过血内皮祖细胞1和2（Haem 1 和 Haem 2），血液祖细胞1和2（BP1 和 BP2），红系1-4（Ery 1–Ery 4）  
![blood_lineage_cluster](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/blood_lineage_cluster.png)  
  
上述轨迹未包含内皮区域（亚分支EC1–EC8）。该区域富集E8.25–E8.5的细胞，呈现复杂的结构，高表达 *Kdr*（故也称为 *Kdr<sup>hi</sup>* 区域），后者编码 FLK1 蛋白。另外，部分内皮细胞表达造血标志基因，如 *Spi1* (also known as *PU.1*) 和 *Itga2b*，可能提示**血液在 yolk-sac-definitive wave 中由内皮产生**  
已知内皮细胞在卵黄囊、尿囊以及胚体中独立生成，且推测尿囊内皮具有特异性转录特征。为了确认 *Kdr<sup>hi</sup>* 区域的异质性是否与不同的解剖学位置有关，从一批新的 E8.25 胚胎中分离出卵黄囊，尿囊以及胚体，流式挑选出 *FLK1+* 群体，得到288个内皮细胞进行 scRNA-seq。将之前EC1-EC8亚分支的细胞分配到最可能的胚胎位置（下图），表明**不同的解剖学结构来源能够解释部分内皮中发现的转录异质性**  
![endothelial_cell_location](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/endothelial_cell_location.png)  
  
之前有研究表明，除了红细胞，primitive wave 还能够生成巨噬和巨核祖细胞，但对这些细胞的分子性质了解甚少。本研究中作者鉴定出2个罕见分支（频率约0.1%），注释为巨核细胞和骨髓细胞  
作者发现骨髓细胞中表达 *Ptprc*（编码CD45），*Kit*，*Csf1r* 和 *Fcgr3*（编码CD16）（下图），均被报道是E8.5时期 EMP 样细胞的标志基因，后者能够产生小神经胶质巨噬细胞，从而与早期骨髓祖细胞能够产生大脑小神经胶质的报道对上；不过没有检测到更加成熟小神经胶质的标志基因 *Cx3cr1*，*Adgre1*（编码F4/80），只有 *Tmem119* 少量表达。通过分离E8.5胚胎区域并基于标志 CD16/32 和 CSF1R 进行流式分析，发现罕见的 CD16/32<sup>+</sup>CSF1R<sup>+</sup> 在所有区域中均存在，表明**E8.5时期该类细胞早已开始从卵黄囊中向外迁移了**  
![myeloid_marker](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/myeloid_marker.png)  
  
# 研究基因突变的平台  
之前的研究强调了转录因子 TAL1（又称SCL）在造血过程中的重要性，*Tal1<sup>−/−</sup>* 小鼠胚胎在E9.5左右由于严重贫血而死亡。由于小鼠繁殖以及确定胚胎基因型十分耗时，且突变的影响常会被明显的发育畸形或胚胎致死掩盖，因此通过敲除小鼠探究活体内这些主要调控基因的机制比较困难。因此，作者通过将 *Tal1<sup>−/−</sup>* tdTomato<sup>+</sup> 小鼠的胚胎干细胞注射到野生型囊胚中，形成嵌合型小鼠胚胎。产生的嵌合体中，野生型细胞仍能生成血细胞，使得我们能够在健康胚胎中研究 TAL1 缺失带来的影响  
为了确定 *Tal1* 突变细胞与特定细胞系的异常是否相关，从E8.5胚胎中将 tdTomato<sup>−</sup>（野生型）和 tdTomato<sup>+</sup>（Tal1<sup>−/−</sup>） 细胞分离，进行单细胞测序，然后将细胞转录组映射至之前获得的野生型图谱上来对细胞进行注释（下图）。tdTomato<sup>+</sup> 细胞中没有血细胞（下图中白色三角形指示的区域），与 *Tal1* 在造血过程中的重要作用对应。另外，作为对照，将野生型 tdTomato<sup>+</sup> *Tal1<sup>+/+</sup>* 胚胎干细胞注射到野生型胚胎中，这些细胞对造血作用的贡献与 tdTomato<sup>−</sup> 原胚胎的细胞相近（意味着确实是 *Tal1* 缺失带来的影响）  
![chimaera_cell_cluster](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/chimaera_cell_cluster.png)  
  
对比野生型和 *Tal1<sup>−/−</sup>* 的细胞映射至前文血液相关细胞群上，发现 TAL1 缺失破坏了原始红系细胞的出现，意味着 **primitive wave 被破坏**（下图，图中黑色箭头表示 *Tal1<sup>−/−</sup>* 细胞中血液发育在该处被封锁）。本研究中新鉴定的巨核细胞和骨髓细胞的出现同样被阻断；另外，虽然有一部分 *Tal1<sup>−/−</sup>* 细胞映射到血内皮亚分支EC6和EC7上，但它们缺少血液发育相关基因的表达，如 *tga2b* 或已知的 TAL1 靶基因 *Cbfa2t3*，与野生型的原细胞相反，意味着 **second haematopoietic wave 也中断了**  
![blood_related_cell](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/blood_related_cell.png)  
  
为了进一步探索 second haematopoietic wave 中断的情况，量化 *Tal1<sup>−/−</sup>* 和野生型细胞对内皮细胞（EC1–EC8）和血内皮细胞（Haem 1 and Haem 2）的贡献比例（下图）  
![abundance_comparison](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/abundance_comparison.png)  
  
上图中我们可以看到EC3亚分支中，E8.5期的 *Tal1<sup>−/−</sup>* 细胞远多于野生型，不过这有可能是映射到该亚分支的 *Tal1<sup>−/−</sup>* 具有与EC3相似但不相同的转录组。为了验证这一假设，将这些 *Tal1<sup>−/−</sup>* 细胞先与整个参考图谱中最相似的细胞比较，再与同样映射至EC3亚分支的野生型原细胞比较，发现了一些差异表达的基因，如下图所示。CM 指的是参考图谱中的心肌细胞群（cardiomyocytes），WT指同样映射到EC3亚分支的野生型细胞，EC3指参考图谱中最相似的细胞。可以看到 *Tal1<sup>−/−</sup>* 细胞中特异性上调表达 *Pcolce*，*Tdo2* 和 *Plagl1* 等基因。在整个参考图谱中检测这些基因，发现它们在间质以及其他中胚层分支，如尿囊，旁轴，咽部和中间中胚层中高表达  
此外，有小部分 *Tal1<sup>−/−</sup>* 细胞上调表达心脏相关基因如 *Nkx2-5*，*Mef2c* 和 *Tnnt2*，与之前报道 Tal1−/− 卵黄囊细胞呈现出心肌细胞样表型相对应。不过这些细胞没有表达所有的心肌细胞转录程序，仍在继续表达内皮基因如 *Esam* 和 *Sox17*，虽然相较野生型存在下调  
以上结果表明，***Tal1* 缺失将细胞封锁在 second haematopoietic wave 上类似于EC3亚分支的转录组状态上；由于无法继续向着血液表型进发，这些细胞转而开始激活其他的中胚层程序**  
![dge](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/dge.png)  
  
# REF
1. Pijuan-Sala B, *et al*. A single-cell molecular map of mouse gastrulation and early organogenesis. Nature, 2019, 566(7745):490-495. doi: 10.1038/s41586-019-0933-9. Epub 2019 Feb 20.  
