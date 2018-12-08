---
title: scRNA-seq 技术介绍
description: 单细胞转录组测序技术汇总
categories:
 - scRNA-seq
tags:
 - scRNA-seq
---

# Smart-seq2 （2013）
**Smart-seq2** 的流程如下图所示：  
![Smart-seq2](/img/2018-08-04-scRNA-seq-technology/Smart-seq2.jpg)  
  
1. 细胞裂解。使用相对较温和的方式（低渗溶液）裂解细胞。裂解液中包含了 free dNTPs 和具有 oligo-dT 的寡核苷酸序列（包含30nt的 oligo-dT 和25nt的通用5'锚定序列），后者能够起始具有 polyA 尾的 RNA 的 RT 反应。free dNTPs 能够提高 RT-PCR 的产量。  
2. RT。当 RT 进行至 RNA 5'端时（合成是从5'-3'，对 RNA 模版而言是从其3'-5'），会在合成的 cDNA 3' 末端添加2-5个 untemplated nucleotides，而模版转换反应（template-switching reaction）正依赖于这几个碱基。TSO（template-switching oligos）3'末端会携带与这几个碱基互补的序列，以及与 oligo-dT 引物相同的锚定序列（使得后面的 PCR 只需要使用单种引物），与RNA 5'端连接，然后反转录酶实现模版转换，在 cDNA 3'端合成与 TSO 互补的序列。  
3. 预扩增。当合成第一条链后，进行有限次数的循环，得到的材料足够进行后续分析即可。  
4. Tagmentation。对扩增的 cDNA 快速高效地构建测序文库，DNA 片段化和接头连接在同一步完成，使用 Illumina 双接头策略。  
5. 扩增并测序。对带接头的 DNA 进行扩增并测序。  
  
与 **SMART-seq** 的差异：  
1. SMART-seq 是在42度下进行RT，但由于部分 RNA 具有二级结构，出现位阻现象，导致链的延伸提前终止。SMART-seq2 中加入了甜菜碱，是一种甲基供体，能够增加蛋白质热稳定性。SMART-seq2 先在50度下反应2 min，破坏 RNA 二级结构，然后回到42度，进行2 min RT，循环多次，提高cDNA产量  
2. 甜菜碱的加入需要提升氯化镁浓度。镁离子能与甜菜碱中的羧酸阴离子结合，作为缓冲液帮助细胞应对渗透压的变化。有报道称过量氯化镁可能对 PCR 保真性有负作用，但在 SMART-seq2 的 RT-PCR 以及 PCR 阶段并未发现  
3. SMART-seq 使用的 TSO 3'端是3个 riboguanosines，并使用 Moloney murine leukemia (M-MLV) 逆转录酶，后者能够进行模版转换，合成 TSO 互补链；SMART-seq2 中前两个碱基还是 riboguanosines，第3个是修饰过的鸟嘌呤得到的 locked nucleic acids (LNA)。LNA 单体具有更强的热稳定性，以及在退火阶段高效结合到cDNA非模版3'端。RT 酶则换成了 Superscript II (Invitrogen)，它是 M-MLV 编辑得到的，Rnase H 酶活性降低，热稳定性增强  
4. SMART-seq2 预扩增时用 KAPA HiFi HotStart ReadyMix 替换了 Advantage 2 polymerase mix，取消了先用磁珠富集 cDNA 这一步，简化了 protocol  
   
# CEL-seq （2012）
**CEL-seq** 通过体外扩增（In vitro transcription, IVT）实现单细胞测序。  
![CEL-seq](/img/2018-08-04-scRNA-seq-technology/CEL-seq.png)  
  
1. RT。引物中包含 poly(T)，barcode，5' Illumina 测序接头和 T7 启动子，结合到具有 ploy(A) 尾的 mRNA 上，合成 cDNA。  
2. 合成第二条链后，形成的就是双链的 DNA。将所有的 cDNA 样本混合到一起，达到足以进行 IVT 的模版量。T7 promoter 用于启动IVT，以第二条合成的 DNA 为模版转录 RNA（所以此时的 RNA 序列与原 mRNA 序列是反向互补的，即 antisense RNA，aRNA），此时5'端最外为 Illumina 5'adaptor。经过多轮 IVT，实现模版的线性扩增  
3. 扩增后的 RNA 直接进入 RNA 模版准备步骤，被打断适合测序的大小，连接上 Illumina 3' adaptor  
4. RNA 反转录成 DNA。两端都具有接头，且含有 barcode 的片段被用于双端测序。read1 记录 barcode 信息，read2 用于确定 RNA  
  
# CEL-seq2 （2016）
![CEL-seq2](/img/2018-08-04-scRNA-seq-technology/CEL-seq2.png)  
  
相较 CEL-seq，做出了一些改进：  
1. RT 引物中加入了一段6nt的 UMI，将8nt 的 barcode 缩短至6nt，并缩短 T7 引物和 Illumina 5' 接头，使引物总长从92nt缩短至82nt，提升了 RT 效率  
2. RT 时使用 SuperScript II，并在第二条链合成时使用 SuperScript II Double-Stranded cDNA Synthesis Kit  
3. 改进了去除 dsDNA 和 aRNA 的过程，提高产量  
4. IVT 得到的 RNA 在反转录时，直接插入接头，不需要进行连接  
  
# Drop-seq（2015）
首先是构造大量不同的 barcoded beads，示意图如下：  
![Drop-seq-barcode](/img/2018-08-04-scRNA-seq-technology/Drop-seq_barcode.png)  
  
直接在每个 beads 上从5'-3'合成寡核苷酸引物，3'用于启动逆转录。每条寡核苷酸链包含4个部分：  
1. PCR handle，一段共有的序列，作为 PCR 和测序的起始位置；  
2. 细胞标签，用于标记细胞身份，同一 bead 中均相同；  
3. UMI，用于排除 PCR 偏好性带来的干扰  
4. oligo-dT 序列，捕获带有 poly(A) 尾的 mRNA，起始 RT  
  
barcode 的合成采用 split-and-pool 策略，将上百万个磁珠分装成4等份，分别加入A，G，C，T，然后全部混合，再随机封装，共计循环12次，得到足够多的不同barcode；之后整个微珠池进行8轮的简并寡核苷酸（degenerate oligonucleotide）合成，得到 UMI  
  
然后让我们来看一下后续的整个分析过程示意图：  
![Drop-seq](/img/2018-08-04-scRNA-seq-technology/Drop-seq.png)  
  
1. 设备中两条通道分别导入含有微珠的细胞裂解液和细胞悬液，穿过油通道，形成微滴  
2. 在微滴中裂解细胞，mRNA 结合到磁珠上，形成 STAMPs (single-cell transcriptomes attached to microparticles)  
3. 加入试剂破坏水油表面，使得微滴破坏，所有的 STAMPs 混合到一起，进行 RT（使用的是模版转换来获取全长 cDNA），生成共价结合，稳定的 STAMPs（cDNA 已经与磁珠上的引物连接在一起了）  
4. 之后就对 cDNA 进行 PCR 扩增，然后进行测序  
  
# inDrop-seq（2015）
首先来看 barcoded hydrogel microspheres（BHM，带条形码的水凝胶微球） 的合成过程，示意图如下：  
![inDrop-seq-BHM](/img/2018-08-04-scRNA-seq-technology/inDrop-seq_BHM.png)  
![inDrop-seq-barcode](/img/2018-08-04-scRNA-seq-technology/inDrop-seq_barcode.png)  
  
其中看上图B中的结构，acrylic phosphoroamidite moiety 是常用的 DNA 引物，后面是一段可光解的间区和 T7 RNA 聚合酶启动子序列，以及测序引物。之后，经过两部的引物扩展（图D，\*代表反向互补序列）：先进行引物延伸，加上了 barcode1 和 W1 site，然后转换成单链 DNA，再进行一次杂交和延伸，加上 barcode2，UMI 和 poly(T)  
  
再来看整一个测序的过程：  
![inDrop-seq](/img/2018-08-04-scRNA-seq-technology/inDrop-seq.png)  
  
1. 在微流设备中有4个输入通道，分别是 BHMs，细胞，RT/lysis 试剂以及油，封装成微滴  
2. 利用 UV 释放 BHMs 上的引物，与裂解细胞得到的 mRNA 结合，进行 RT  
3. 破坏微滴后将所有 cDNA 整合，按照 CEL-seq protocol 进行测序  
  
# CITE-seq （2012）  
**Cellular indexing of transcriptomes and epitopes by sequencing (CITE-seq)** 利用寡核苷酸标记的抗体，对细胞表面蛋白和和转录组同时进行测定，其流程如下图所示：  
![CITE-seq](/img/2018-08-04-scRNA-seq-technology/CITE-seq.jpg)  
  
1. 如图a所示构造 DNA-barcoded antibodies。通过链霉亲和素-生物素互作，将寡核苷酸5'端与抗体相连，而这种通过二硫键的结合会在还原条件下断裂，释放寡核苷酸（图b所示）。寡核苷酸中包括用于 PCR 扩增的 handle，用于抗体身份识别的 barcode，以及能与 oligo-dT 引物结合的一段 polyA 序列。  
2. 将 antibody-oligo 复合物与单细胞悬液共孵育，条件设置参照流式细胞染色，之后将未结合的抗体移除，然后进行 scRNA-seq。  
3. 以 Drop-seq 为例（beads 构成如图c所示），在微流体设备中将单细胞封装成小液滴，然后细胞裂解，使得带有 oligo-dT 的 beads 与细胞 mRNA 及抗体-寡核苷酸复合物结合，然后对同一细胞的 mRNA 和寡核苷酸进行 RT 和 PCR（图d），并标记上相同的 cell barcode。  
4. 扩增得到的 cDNA 和 antibody-derived tags (ADTs) 可以通过大小进行区分，从而分别进行 Illumina 测序文库的制备。由于两种文库可以独立生成，因此可以调整其在单个 lane 中的比例以确保各自所需的测序深度。
  
# REAP-seq  
**RNA expression and protein sequencing assay (REAP-seq)** 实现了在单个细胞中同时测定蛋白质和 mRNA。细胞通过连接了 DNA barcode 的抗体进行标记，其中抗体 DNA 标签的设计示意图如下所示：  
![REAP_AbB](/img/2018-08-04-scRNA-seq-technology/REAP_AbB.jpg)  
  
可以看到抗体与长为65-66bp的寡核苷酸偶联。寡核苷酸包含3个部分：（1）33 bp Nextera Read 1 sequence，作为后续扩增和测序的引物；（2）8bp长的独一无二的抗体标签（antibody barcode, Ab BC）；（3）24-25 bp poly(dA) sequence。  
REAP-seq 测序流程如下所示：  
![REAP_seq](/img/2018-08-04-scRNA-seq-technology/REAP-seq.jpg)  
  
1. 将 antibody DNA label 与细胞共孵育，之后洗脱未连接的 AbBs。  
2. 如图a所示，将标记了抗体-DNA 标签的细胞用 10x Genomics single-cell (sc)RNA-seq platform 进行处理。通过 DNA 聚合酶实现同时扩增 AbBs 和 mRNA，并利用核酸外切酶Ⅰ将过量的未结合的单链寡核苷酸消化，防止不同细胞间的 barcode 与 AbBs 结合。  
3. 破坏液滴后通过大小对 AbBs 和 mRNA 进行分选，然后制备文库和测序。  
  
`REAP-seq` 与 `CITE-seq` 的差异在于前者在抗体和胺化了的 DNA barcode 之间构建了稳定的共价键，能够减小空间位阻。减小空间位阻对于蛋白试验的扩展性至关重要，将来也能拓展至细胞内标记。  
