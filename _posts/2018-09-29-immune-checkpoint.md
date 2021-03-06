---
title: Immune checkpoint therapy
description: 免疫检查点治疗相关汇总
categories:
 - cancer
tags:
 - cancer
 - immune checkpoint
 - immunotherapy
---

# Introduction
免疫检查点治疗，通过靶向 T 细胞调控通路来增强抗肿瘤免疫应答，为我们提供了对抗肿瘤的新武器。相较传统的肿瘤治疗方法，免疫检查点抑制药物并非直接针对肿瘤细胞，而是靶向参与 T 细胞调控的分子，通过移除封堵了有效抗肿瘤 T 细胞应答的抑制信号通路，激活免疫系统，清除肿瘤细胞。目前，检查点抑制剂包括 CTLA-4 抑制剂 和 PD-1/PD-L1 抑制剂两大类，本篇文档将对肿瘤免疫检查点治疗方法相关内容做一个汇总记录  
  
# CTLA-4
T 细胞活化需要两个信号（下图）：  
* **T 细胞抗原受体（TCR）与抗原-MHC复合物结合**  
* **T 细胞表面的 CD28 与抗原提呈细胞（APC）表面的 B7 分子（CD80 和 CD86）结合**  
![costimulatory](/img/2018-09-29-immune-checkpoints/costimulatory.png)  
  
其中，B7 分子仅在造血细胞，尤其是树突状细胞（DC）中表达。除了少数淋巴瘤特例外，肿瘤细胞并不会表达 B7 分子，因此会被免疫系统忽视；不过这可以通过炎症反应克服，例如杀死肿瘤细胞，使得 APC 能够结合并提呈抗原，并提供 B7 分子，从而激活 T 细胞。当肿瘤抗原参与到 B7 共激活中时，激活的肿瘤特异性 T 细胞获得效应功能，迁移至肿瘤位点，攻击肿瘤细胞  
`毒性 T 淋巴细胞相关蛋白4（cytotoxic T lymphocyte–associated protein 4, CTLA-4）`分别被 James Allison 和 Jeffrey Bluestone 两个研究组发现可能对 T 细胞应答具有抑制作用。*ctla-4* 基因与 CD28 具有高度同源性，而 CTLA-4 与 CD28 一样，**能够结合 B7 分子，且亲和力更强**。在休眠的 T 细胞中，CTLA-4 作为一个细胞内蛋白存在；但在 TCR 以及 CD28 共刺激信号参与下，CTLA-4 转移至细胞表面，并与 CD28 竞争与 CD80，CD86 等共刺激分子的结合，从而抑制 T 细胞的扩增和活化。Allison 提出，假如能够利用抗体将该分子刹车封锁片刻，可能能使 T 细胞增殖并活化  
经过临床前试验发现 CTLA-4-blocking 抗体确实能在动物体内产生作用。之后，两个全人源 CTLA-4-blocking 抗体（ipilimumab and tremelimumab）进入临床试验，发现确实能够产生持续性的肿瘤退化，虽然并非所有病人都有如此效果。除此之外，该疗法往往还伴随一系列由于组织特异性炎症，包括小肠结肠炎，肝炎和皮炎等导致的毒副作用，不过这些炎症反应可以通过皮质类固醇或者其他免疫抑制剂来控制，且对抗肿瘤活性没什么太大的影响。至于另一些相对较罕见的甲状腺，垂体以及肾上腺等炎症副作用，就需要终身激素替代  
ipilimumab 治疗后，放射影像显示的临床反应与常规的治疗方法有一些区别：部分患者肿瘤可能会先继续生长，或者先出现新的肿瘤，然后消失，之后原来的肿瘤才会缩小，这就给基于过去常用的客观反应率或者无进展生存期等指标进行安全监管批准带来了挑战，因此需要用更长期的结果--总体存活率来作为注册试验研究的主要终点。FDA 在2011年批准了 ipilimumab，Tremelimumab 仍处于临床试验观察中，另外有一个新的 CTLA-4-blocking 抗体进入临床试验  
CTLA-4 抑制剂相对较低的反应率以及相对发生率较高的毒副作用使得找到预测性和药效标志物成为研究的优先方向。研究表明更高的肿瘤突变负荷（tumour mutation burden, TMD）更有可能发生应答，此外外周血中绝对淋巴细胞数量和诱导性共刺激分子 ICOS 也与最终的治疗效果有关。另外，虽然大量小鼠研究发现具有正确 Fc 结构域的 CTLA-4-blocking 抗体能够消耗退化肿瘤中的调节 T 细胞（Tregs），但在人体中，相关的数据还是不足，因此仍需进一步的研究  
  
# PD-1
程序性细胞死亡1（programmed cell death 1, PD-1）受体与其表达于肿瘤细胞表面的配体 PD-L1 结合后成为抗肿瘤 T 细胞效应作用的负向调控剂。PD-1 最初因能够诱导活化的 T 细胞杂种瘤凋亡而得名，但后续研究发现它实际上是一个免疫检查点，其抑制功能受酪氨酸磷酸酶 SHP-2 调控，后者能使 TCR 下游信号分子去磷酸化。PD-1 有两个配体，PD-L1（CD274 或 B7-H1），广泛表达于暴露于促炎细胞因子的体细胞；PD-L2（CD273 或 B7-DC），仅表达于抗原提呈细胞。肿瘤微环境中炎症诱导的 PD-L1 表达导致 PD-1 调控的 T 细胞耗竭，从而抑制抗肿瘤毒性 T 细胞应答  
随着肿瘤从原位病灶转移，抗肿瘤 T 细胞持续识别同源的肿瘤抗原。TCR 触发导致产生促炎细胞因子，包括最强有力的 PD-L1 表达刺激分子--干扰素γ（IFN-γ）。T 细胞长期暴露于同源抗原导致靶细胞表达 PD-L1，而 T 细胞内持续性的 PD-1 信号诱导 T 细胞耗竭的表观程序。其他还有一些 PD-1 信号通路的功能尚不清楚，PD-L1 能够结合 T 细胞上表达的共刺激分子 CD80（B71），传递抑制信号；排斥守护分子 b（repulsive guidance molecule b, RGMb）能与 PD-L2 结合，可能与肺部耐受性有关  
缺乏 PD-1 的小鼠的表型相较缺乏 CTLA-4 的小鼠更有限，几乎没有自身免疫疾病，意味着 PD-1 的作用有限，因此阻断 PD-1 信号通路对抗肿瘤 T 细胞的作用更加特异，相较阻断 CTLA-4，治疗效果更好，毒性更小  
  
# REF
1. Sharma P, *et al*. The future of immune checkpoint therapy. Science. 2015, 348:56-61. doi: 10.1126/science.aaa8172.
