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

# 基本概念  
定形内胚层在原肠胚形成过程中出现并取代胚外内脏内胚层，参与肠管形态发生过程，并构成相关的内脏器官。定形内胚层插入到原始内胚层，使得原始内胚层细胞扩散。之前研究认为原始内胚层会被定形内胚层完全取代，但后续发现即便是到E8.75阶段，在胚胎肠道中仍存在部分原始内胚层来源的细胞，尤其是后肠部分，占到30%  
  
# Introduction
小鼠E6.5-E8.5覆盖了原肠胚形成和早期器官发生的关键阶段。该窗口期，上胚层（epiblast）细胞分化成形成所有主要器官的外胚层（ectodermal）、中胚层（mesodermal）和内胚层（endodermal）细胞。为了研究该阶段细胞动态过程，作者在E6.5-E8.5之间每隔6小时进行取样，共收集411个小鼠胚胎进行 scRNA-seq。整个数据集囊括了 TS9，TS10，TS11 和 TS12（Theiler stages，Theiler等人基于一组形态学标准划分的小鼠发育阶段），分别富集于 pre-streak to early streak，mid-streak to latestreak，neural plate 和 headfold to somitogenesis 阶段  
  
# 早期胚胎发生的单细胞图谱  
严格质控后最终得到116,312个细胞，每个细胞检测到的基因中位数为3436。聚类得到**37个主要的细胞群体**（下图左），**不同类型的细胞的出现与采样的时间点存在关联**。可以看到，多能外胚层细胞的比例随着时间减少；中胚层和最终内胚层细胞最早在E6.75时期即出现；从E7.5开始，随着器官发生过程开启，外胚层细胞出现，并伴随各胚层细胞类型的多样化（下图右）  
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
  
作者利用 DPT 将细胞排序，并将基因按照其表达沿着轨迹的动态变化过程进行聚类。基因沿着两条轨迹变化的图式都分成两块区域，其中第一区域对应在E7.5完成内胚层插入之前的基因表达情况，第二区域则对应后续时间点的基因表达（原文献 Extended Data Fig. 5c）。结果显示，内脏内胚层-后肠1轨迹中，内脏内胚层基因在第一区域上调，之后随着细胞转向肠道命运，急剧减少（下图，左边粉红色块对应内胚层插入，右边蓝色块对应肠道成熟区域）。这一结果表明有**一组内脏内胚层细胞，在终末内胚层插入之前，正进行内脏成熟过程**  
![VE_gene_dynamics](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/VE_gene_dynamics.png)  
  
在两条轨迹中，有一组基因在插入过程中均上调表达，包括参与上皮重塑的基因如 *Pcna*，*Epcam* 和 *Vim*，与该阶段预计出现的上皮排列相契合。后续肠道成熟以及形态发生过程，共同上调的基因则富集于转录因子，当中的66%为呈现出顺序激活的同源结构域蛋白，表明在后肠特化过程中的时间共线性。此外还找到了内脏内胚层-后肠1早期特异性表达的基因，包括 *Hes1*，*Pou5f1* 和 *Sox4*，可作为后续研究的候选者  
  
# 血内皮细胞系的起源
卵黄囊中红细胞可分别通过两种方式形成。一是（primitive）从E7.5开始，产生带核的红细胞，出生后不久便消失，二是（yolk-sac-definitive）在E8.25以血管内皮出现红系骨髓祖细胞（erythro-myeloid progenitors, EMP）起始，后续迁移至胎儿肝脏形成终末红细胞。虽然已经知道了原始和卵黄囊终末造血过程在表型以及分子层面的一些关键性差异，但是对其各自在体内的祖细胞仍知之甚少。为深入研究这些过程，将红细胞，红系骨髓祖细胞，血细胞祖细胞，内皮细胞以及混合中胚层细胞（共15875个细胞）分离出来重新进行聚类，结果如下图所示（BP, blood progenitor; EC, endothelial cell; Ery, erythrocyte; Haem, haemato-endothelial progenitor; Mes, mesodermal cell; Mk, megakaryocyte, 巨核细胞; My, myeloid cell），其中左图内方框是对骨髓，巨核细胞以及血液内皮细胞部分的放大，右边两图则是表示亚分支之间关系的简化图。可以看到  
![blood_lineage_cluster](/img/2019-04-30-single-cell-molecular-map-mouse-gastrulation/blood_lineage_cluster.png)
