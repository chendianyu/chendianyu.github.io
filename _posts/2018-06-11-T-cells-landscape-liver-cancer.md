---
title: 单细胞测序揭示肝癌中浸润 T 细胞图谱
categories:
 - Cancer
tags:
 - scRNA-seq
 - cancer
 - liver cancer
---  
  
针对肿瘤浸润淋巴细胞进行系统性的研究是发展免疫治疗以及预测其在肿瘤临床治疗效果的关键。本文献中，作者从6个肝癌患者的外周血，肿瘤及邻近正常组织在分离出5063个 T 细胞，进行高深度 scRNA-seq。利用这些细胞的转录组图谱和组装的 T 细胞受体序列，作者鉴定出11个 T 细胞亚组，描绘了其发育轨迹，并确定了每个亚组的特征基因。这些转录组数据为我们理解癌症中的免疫图谱提供了非常重要和丰富的信息。  

<!-- more -->

# Highlights

# Introduction  
过去十年，肿瘤免疫治疗显著的改变了癌症治疗的面貌。但是，虽然带来了显著的临床应答，免疫疗法在不同病人以及不同类型癌症之间的效果并不一致。因此，我们需要找到能够有效预测治疗效果的生物标志物。尽管出现了诸如突变负荷 (mutation load)，肿瘤浸润淋巴细胞 (tumor-infiltrating lymphocytes, TILs) 以及药物特异性靶标的表达量等指标，其预测作用仍不够可靠。  
新型癌症免疫疗法的开发和有效生物标记的确定需要我们对肿瘤中存在的 T 细胞有深入的理解。针对结直肠癌，肺癌以及乳腺癌等中浸润调控T细胞 (Tregs) 的转录组分析揭示了这些细胞的高度抑制性，它们能将免疫抑制剂解放出来的部分活性 CD8+ cytotoxic T cells 很快又压制回无活性状态，或者使之功能失调，称为 `exhaustion`。后者的特征是效应功能低下，持续性表达抑制性受体以及独一无二的转录状态。因此，对**导致肿瘤中 CD8+ T cell exhaustion 增强和 Tregs 累积**的机制和通路的深入理解将有助于提供更好的策略来协调免疫系统清除癌症。  
对浸润淋巴细胞的单细胞分析能使得我们对这些细胞有更细致的认识。除了能够揭示 T 细胞耗竭特征及其与 T 细胞活化间的联系外，scRNA-seq 还能确定每个细胞中的 TCR 序列。这些序列对于识别病毒抗原或者肿瘤细胞 MHC 提呈的肿瘤特异性新生抗原至关重要。尽管 TCRs 数量众多，且存在随机重组，但**通过寻找一致的 TCR 序列能够帮助我们理解 T 细胞克隆扩展模式和 T 细胞谱系**。  
肝癌 (hepatocellular carcinoma, HCC) 是世界范围内肿瘤相关死亡的主要原因之一，但截至目前，治疗手段有限，且临床上缺乏免疫疗法成功的案例。尽管肝癌中有大量 TILs，但据推断这些细胞并不能杀死肿瘤细胞。因此 HCC 是一个用于描述 TILs 功能紊乱的理想模型。  
  
# 肿瘤特征和单个 T 细胞转录组数据的获得  
作者首先对6位未经过治疗的 HCC 患者的肿瘤样本进行混池外显子组和转录组测序来获得他们的基本特征，并利用 Opal multicolor IHC staining 确定了不同亚型 T 细胞的存在。其中，作者发现肿瘤中心 CD8+ T cell 的数量远低于外皮层，而 FOXP3+ regulatory T 细胞则没有这一现象，表明 **CD8+ T cell 浸润效率低**。  
之后，作者对 CD3<sup>+</sup>CD4<sup>+</sup> 和 CD3<sup>+</sup>CD8<sup>+</sup> T 细胞进行分选，并将 CD4<sup>+</sup> T 细胞进一步分为 CD25<sup>-</sup> 和 CD25<sup>high</sup>。肿瘤中 CD4<sup>+</sup>CD25<sup>high</sup> 群体占所有 CD4<sup>+</sup> 细胞的比例要远高于正常组织和外周血，表明**肿瘤微环境中可能存在 Tregs 的富集**。  
为方便起见，作者将分别来自外周血，相邻正常组织以及肿瘤组织的 CD8<sup>+</sup> cytotoxic T cells 记为 PTC，NTC 和 TTC；CD4<sup>+</sup>CD25<sup>-</sup> T cells 记为 PTH，NTH 和 TTH；CD4<sup>+</sup>CD25<sup>high</sup> T cells 记为 PTR，NTR 和 TTR。  
分选出来的5063个细胞的转录组数据结果如下图所示：  
![scRNA_seq](/img/2018-06-11-T-cells-landscape-liver-cancer/scRNA_seq.png)  
  
从图中我们可以看到：  
* 总体来看，细胞是按照来源和细胞类型进行聚类的，表明**不同组细胞具有不同的特征和内在结构**  
* 一些情况下，不同组的细胞会混在一起，如部分 TTC 和 NTC，表明它们有相似的功能属性  
* PTC 细胞被分成两个明显的分支，与血液中同时存在初始和活化 CD8<sup>+</sup> T 细胞一致  
* TTH 表现出高度多样性，其中小部分与 TTR 类似  

以上发现的 T 细胞复杂的内部组成表明**通过单细胞技术仔细剖析肿瘤相关 T 细胞很有必要**。  
  
# T 细胞聚类和亚型分析  
利用 SC3 包中的 `spectral clustering method` 对 T 细胞进行无监督聚类，最后得到11个类，其中5个为 CD8+ 细胞，另6个为 CD4+ 细胞，每个类别均有其特异的特征基因。结果如下图所示：  
![subtypes](/img/2018-06-11-T-cells-landscape-liver-cancer/subtypes.png)  
![CD8_signature_genes](/img/2018-06-11-T-cells-landscape-liver-cancer/CD8_signature_genes.png)  
  
**不同的 CD8+ 细胞簇在不同病人各组织中分布的模式具有可比性**。如下图所示，表达 CCR7 的初始 T 细胞（C1）和表达 CX3CR1 效应记忆 T 细胞（C2）主要在外周血中相对较为常见；而 MAIT 细胞（C3）在邻近的正常肝组织所占比例最高，其在 HCC 中所占比例显著降低。对 TCGA 队列分析发现 HCC 肿瘤中 MAIT 细胞标志基因 SLC4A10 的表达显著低于对应正常样本，表明 **MAIT 细胞的减少在肝癌中很常见**。此外，TCGA 队列中，**HCC 低 SLC4A10 表达量与不良预后相关**。综上，这些结果表明 **HCC 肿瘤微环境中一大特征便是 CD8+ T 细胞的差异化分布**。  
![CD8_distribution](/img/2018-06-11-T-cells-landscape-liver-cancer/CD8_distribution.png)  
  
类似地，作者找出了6种不同的 CD4+ T 细胞簇，不同类型细胞在组织中的分布情况不同。其中 C8_CD4_CTLA4 是一种 FOXP3<sup>+</sup> Tregs，高表达 FOXP3 及其他多种 Treg 相关基因如 TNFRSF9，TIGIT 和 CTLA4，其主要组成是从肿瘤中分选出的 CD4<sup>+</sup>CD25<sup>high</sup> 细胞 （TTR）。这些细胞更倾向于在肿瘤中富集，如下图所示。  
![C8_CD4-CTLA4](/img/2018-06-11-T-cells-landscape-liver-cancer/C8_CD4-CTLA4.png)  
  
# 与 HCC 浸润 Tregs 和 耗竭 CD8<sup>+</sup> T细胞相关基因的鉴定和确认  
从上面的结果我们知道 CTLA4<sup>high</sup> Tregs **(C8_CD4_CTLA4)** 在肿瘤中富集（上图），另外根据下图，我们可以看到，高表达耗竭标志基因 CTLA4，PDCD1 和 HAVCR2 的 **C4_CD8-LAYN**（耗竭 CD8<sup>+</sup> T 细胞） 亚型同样在肿瘤中富集。 由于这两种 T 细胞亚组表达的共抑制受体如 PDCD1 和 TIGIT 等是癌症免疫治疗的靶标，因此作者主要针对这些类型的细胞进行后续的分析。  
![CD8_cancer_enrichment](/img/2018-06-11-T-cells-landscape-liver-cancer/CD8_cancer_enrichment.png)  
  
作者共鉴定出401个肿瘤浸润 Tregs 特异性表达的基因，包括 FOXP3，CTLA4，TNFRSF18，TNFRSF4 和 CCR8 等。这些基因与之前其他针对黑色素瘤，乳腺癌，结直肠癌以及肺癌的研究所发现的基因存在明显的重合。33个之前发现的常见 Treg 特征基因，有31个在本次研究中也被找到；剩下两个表达量提升相对较低。因此，尽管病人样本数量有限，**单细胞测序仍能够对这些 TILs 进行详细的刻画**。此外，有146个基因是本数据集中特有的。  
之后，作者对肿瘤浸润耗竭 CD8<sup>+</sup> T 细胞也做了分析，对比耗竭和未耗竭 TTC 细胞，得到了82个 耗竭特异性基因。不少已知的耗竭标志基因如 HAVCR2，PDCD1，ENTPD1，CTLA4，TIGIT，TNFRSF9 和 CD27 等都排名靠前；除此之外也包含了相对描述较少的基因如 MYO7A，WARS 和 CXCL13，及新的标志基因如 LAYN，PHLDA1 和 SNAP47。利用这些基因评估耗竭状态，发现**处于晚期的病人相对其他人呈现出更高的耗竭水平**。另外值得注意的一点是，有22个基因在肿瘤浸润 Tregs 中也高表达。  
综上，本研究**既验证了之前研究所发现的与肿瘤浸润耗竭 CD8+ T 细胞和 Tregs 相关的基因，还进一步发现了更多的标志基因**。  
  
# 血液分离出来的 CD8<sup>+</sup> T 细胞和 Tregs 一旦激活，诱导 LAYN 表达  
LAYN，编码透明质酸受体抗体，最近被报道在从肺癌和结直肠癌分离出来的 Tregs 中高表达，但其功能尚未被充分研究。根据本研究的数据，**LAYN 在肿瘤 Tregs (C8_CD4-CTLA4) 和耗竭 TTC 细胞 (C4_CD8-LAYN) 中特异性高表达**。此外，通过对 TCGA HCC 存活率数据与标准化后的 LAYN 表达量的相关分析发现，**LAYN 的高表达量与无病存活率的降低相关**。  
由于 LAYN 编码细胞表面蛋白，因此作者利用流式细胞仪对从人类外周血单核细胞 （PBMCs）分离得到的 CD8<sup>+</sup> 和 Tregs 中该蛋白的表达量进行了验证分析，发现在静止期，LAYN 不在 CD14<sup>+</sup> myeloid cells，B cells，CD4<sup>+</sup> T cells，Tregs 或 CD8<sup>+</sup> T cells 中表达；但当 T 细胞被抗 CD3 和 抗 CD28 抗体激活两天后，可以很容易地在超过30%的 Tregs 和 CD8<sup>+</sup> T 细胞中检测到 LAYN，虽然其在 CD4<sup>+</sup> T 细胞中的表达量仅微微上调。这些数据与之前发现的 LAYN mRNA 在肿瘤浸润 CD8<sup>+</sup> 和 Tregs 表达上调相符。  
虽然 FOXP3 是人类 T 细胞中调控 Treg 分化和功能的主要转录因子，但其在 CD4<sup>+</sup> T 细胞中也能瞬时表达上调；另一个转录因子 Helios，最近被报道能够稳定 FOXP3 的表达以及 Tregs 的抑制作用。因此，作者检测了 FOXP3<sup>+</sup> 和 FOXP3<sup>+</sup>Helios<sup>+</sup> Tregs 中 LAYN 的表达情况，发现其更倾向于在双阳性细胞中表达上调，表明**LAYN 的表达与更具抑制性和稳定的 Tregs 相关**。  
由于之前从未在 HCC 中报道过 LAYN 与 肿瘤浸润耗竭 CD8<sup>+</sup> T 细胞之间的联系，因此作者又做了体外试验，发现 CD8 T 细胞被激活1天后，就诱导表达了 LAYN 蛋白，且能持续存在之第6天。在 PD-1<sup>+</sup> 和 PD-1<sup>-</sup> CD8<sup>+</sup> T 细胞中，均能检测到 LAYN。值得注意的是 LAYN 仅在 LAG-3 阴性 CD8<sup>+</sup> T 细胞中表达。之前研究表明 LAG3<sup>+</sup>CD8<sup>+</sup> T 细胞代表耗竭 CD8+ T 细胞的一种独特的亚型，因此 LAYN 和 LAG-3 的互斥表明 **LAYN 可能是另一种 CD8<sup>+</sup> T 细胞亚型的标志基因**。  
为了揭示 LAYN 在 CD8+ T 细胞中可能起到的作用，作者利用逆转录病毒调控 LAYN 在上述细胞中高表达以模拟肿瘤浸润 CD8+ T 细胞中的情境，发现 LAYN 过表达的细胞相较对照组生成的 IFN-γ 明显减少，从而支持其**发挥抑制作用**的角色。  
上述结果表明 **LAYN 蛋白能在 Tregs 和 CD8 T 细胞中表达，发挥抑制作用，与之前转录组结果相符**  
  
# 利用 TCR 揭示 HCC 中 CD8 T 细胞和 Tregs 的富集  
TCR 经常被用作确定 T 细胞祖源的唯一标识符，因此作者用 TraCeR 方法，将从5位 HBV 阳性的 HCC 患者细胞获取的单细胞转录组数据进行拼接，得到全长 TCR α 和 β 序列。虽然绝大部分细胞表达唯一的 TCR α 和 β 等位基因，但仍有一部分 T 细胞中的 α 和/或 β 序列并非独一无二。为了准确定义 T 细胞克隆性，作者定义至少有一对相同的 α-β 链的 T 细胞克隆自同一祖先，扩增得到的克隆则被定义为在给定细胞群体中，至少有三个细胞共享相同的 TCR α 和 β 序列。基于这些定义，作者发现 T 细胞间不同 TCR α 链复现的频率和 TCR β 链存在显著相关性，表明这些细胞起源于同一祖先。此外作者还进一步研究了 TCR 重排，以确定 V-D-J 利用偏好性，也确实找到了预期中的 TCR 片段使用偏好性。  
获得每个细胞 TCR 序列为我们研究不同 T 细胞之间关系提供了机会。虽然绝大部分细胞包含独一无二的 TCRs，但仍能在不同类型细胞中发现不同程度的重复利用现象，尤其是肿瘤浸润 CD8<sup>+</sup> T 细胞和 Tregs。**相较于血液和正常肝组织中的 CD8<sup>+</sup> T 细胞，肿瘤组织中的细胞存在克隆性 TCRs 所占比例明显更高**。类似的，**肿瘤浸润 Tregs 同样呈现出克隆富集**，尤其是 C8_CD4-CTLA4 这一分支。  
肿瘤组织中 T 细胞相较血液和正常组织的克隆富集看起来可能是由于**肿瘤微环境下局部性的 T 细胞增殖和活化导致**，这和之前的研究结果一致。但我们仍无法排除这些克隆 T 细胞是在外周淋巴器官中扩增后迁移至肿瘤这一可能性。由于这些病人都是 HBV 阳性患者，克隆 T 细胞也可能是在慢性感染下，T 细胞活化以应对病毒抗原得结果，不过这样的话相邻 HBV<sup>+</sup> 正常组织也应该有类似的克隆富集现象，但事实上并没有，所以排除了这一可能。  
综上，**TCR 信息与不同状态 T 细胞之间存在关联。在 HCC 中，CD8+ T 细胞和 Tregs 可能存在克隆扩增**。  
  
# 利用拟时间状态转换及克隆性 TCRs 寻找亚组之间的关联
获得大量 T 细胞的转录本数据和 TCR 信息，使得我们能够深入理解这些细胞的功能状态和相互之间的关系。作者利用 `Monocle 2` 算法将 CD8<sup>+</sup> T 细胞或 CD4<sup>+</sup> T helper 细胞拟时间排序研究其发展轨迹。此外还用了 EMBEDRR，SCORPIUS 和 TSCAN 等算法进行比较。  
对于 CD8<sup>+</sup> T 细胞，作者移除了 MAIT（C3）细胞，因为它们的 TCR 特征不同。根据表达谱相似性，不同的细胞簇形成了一个拟时间化的相对发展过程。如下图所示，以 C1_CD8-LEF1 cells (naive CD8<sup>+</sup> T cells) 开始, 之后是 C2_CD8-CX3CR1 (effector memory CD8<sup>+</sup> T cells)， C5_CD8-GZMK, 至 C4_CD8-LAYN cells (exhausted CD8<sup>+</sup> T cells) 为止。所以 **T 细胞状态是从活化到耗竭为止**。  
![CD8_pseudotime](/img/2018-06-11-T-cells-landscape-liver-cancer/CD8_pseudotime.png)  
  
从上图中我们看到，C5_CD8-GZMK 细胞呈现出处于效应和耗竭 T 细胞的中间状态。通过对其进行克隆性分析，发现那些来自同一祖先的相同的 TCRs 分布在不同的细胞簇中，如下图所示。在共计61个至少两种细胞簇共享的 TCR 对中，有30种是 C2_CD8-CX3CR1 和 C5_CD8-GZMK 共享的，有20种则是 C5_CD8-GZMK 和 C4_CD8-LAYN 共享的，而其他细胞簇之间共享的 TCR 对不超过5个，进一步证明 **过渡期的 C5_CD8-GZMK 与效应细胞和耗竭细胞之间的发育联系**。  
![CD8_shared_TCR](/img/2018-06-11-T-cells-landscape-liver-cancer/CD8_shared_TCR.png)  
  
类似地，作者对 CD4<sup>+</sup> T helper 也做了分析，结果如下图所示。C6_CD4-CCR7 naive T cells 和 C9_CD4-GZMA T helper cells 在主干处聚集，而 C10_CD4-CXCL13 (exhausted CD4<sup>+</sup> T cells) 和 C11_CD4-GNLY (cytotoxic CD4<sup>+</sup> T cells 则分居拟时间轨迹的不同方向，说明这两者之间的功能分化。共享 TCR 序列分析发现绝大部分共享 TCRs 发生在 C9_CD4-GZMA 和 C10_CD4-CXCL13 以及 C9_CD4-GZMA 和 C11_CD4-GNLY 之间；而 C10_CD4-CXCL13 和 C11_CD4-GNLY 之间几乎没有。  
![CD4_T_helper_pseudotime](/img/2018-06-11-T-cells-landscape-liver-cancer/CD4_T_helper_pseudotime.png)  

基于拟时间轨迹和 TCR 的分析，作者推测**耗竭 CD4 和 CD8 T 细胞分别与具有 GZMA 和 GZMK 标志基因的中间态细胞关联更紧密，而不是效应细胞群体**。至于将这些中间态细胞作为免疫治疗靶标，以使它们向效应细胞而不是耗竭细胞方向转化这一策略是否有效，仍需要更多的研究。  

# 耗竭 CD8+ T 细胞和 Tregs 在 HCC 微环境中的克隆富集  
过去公认肿瘤微环境中的 T 细胞倾向于成为耗竭或者受 Treg 抑制，从而避免这些细胞诱发 T 细胞介导的肿瘤细胞杀伤。本研究实验数据发现，肿瘤微环境中效应 CD8+ T 细胞数量极少，而耗竭状态的细胞更多。此外，相较于早期患者，晚期患者中耗竭 CD8+ T 细胞呈现增长的趋势。此外，如下图所示，克隆性 CD8+ T 细胞更可能表现出耗竭的表型，尤其是在晚期 HCC 肿瘤中；高克隆性群体（每个克隆细胞数多于4个）相较非克隆性群体更可能表现为耗竭状态。  
![colonal_TTC](/img/2018-06-11-T-cells-landscape-liver-cancer/colonal_TTC.png)  
  
除了耗竭 CD8+ T 细胞，之前有研究在结直肠和前列腺癌中检测出 CD8+FOXP3+ 调节 T 细胞。在本研究中，也发现有少量耗竭 CD8+ T 细胞表达 Treg 标志基因 FOXP3。根据 TCR 分析，这些 FOXP3+ 细胞中有部分与 FOXP3- 耗竭 CD8+ T 细胞具有相同的 TCR，表明**这些 FOXP3+ 细胞在发育上与典型的 FOXP3- 耗竭 T 细胞存在关联**。此外，如下图所示，这些细胞的表达谱与 Treg 和 cytolytic 相似，能够表达 Treg 标志基因如 FOXP3，CTLA4，TNFRSF18 和 TNFRSF9 以及 溶细胞相关基因 PRF1，GZMA 和 NKG7，表明这些细胞**同时具有抑制和杀伤细胞的能力**。综上表明，**HCC 微环境能够推动浸润 CD8+ T 细胞向耗竭状态转换，且偶尔能使之获得抑制功能**。  
HCC 微环境中高比例的克隆性扩增的 Treg 可能代表另一种抑制效应 CD8+ T 细胞介导的杀伤作用的机制。作者对肿瘤浸润性 CD4+ T 细胞也做了拟时间和 TCR 分析（下图），发现虽然肿瘤浸润性 Tregs (C8_CD4-CTLA4) 与 耗竭 CD4 T 细胞 (C10_CD4-CXCL13) 之间在发育轨迹上很接近，但仅共享极有限的 TCRs。  
![CD4_pseudotime](/img/2018-06-11-T-cells-landscape-liver-cancer/CD4_pseudotime.png)  
  
绝大部分（71/87）肿瘤特异性 Treg 克隆是唯一的，表明其与其他肿瘤浸润 CD4 T helper 细胞之间具有相互独立的发育路径。这与肿瘤相关耗竭 CD8+ T 细胞不同，后者与肿瘤中 CD8 T 细胞的其他簇共享大量常见的 TCRs。Tregs 和耗竭 CD4 T 细胞簇均只有少量常见 T 细胞克隆表明耗竭 CD4 T 细胞可能向可诱导 Tregs 转换。
  
# REF  
1. Zheng C, *et al*. Landscape of Infiltrating T Cells in Liver Cancer Revealed by Single-Cell Sequencing. Cell, 2017, 169:1342-1356. doi: 10.1016/j.cell.2017.05.035.
