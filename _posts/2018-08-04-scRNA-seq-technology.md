---
title: scRNA-seq 技术介绍
description: 单细胞转录组测序技术汇总
categories:
 - scRNA-seq
tags:
 - scRNA-seq
---

# Smart-seq2
**Smart-seq2** 的流程如下图所示：  
![Smart-seq2](/img/2018-08-04-scRNA-seq-technology/Smart-seq2.jpg)  
  
1. 裂解。使用相对较温和的方式（低渗溶液）裂解细胞。裂解液中包含了 free dNTPs 和具有 oligo-dT 的寡核苷酸序列（包含30nt的 oligo-dT 和25nt的通用5'锚定序列），后者能够起始具有 polyA 尾的 RNA 的 RT 反应。free dNTPs 能够提高 RT-PCR 的产量。  
2. RT。当 RT 进行至 RNA 5'端时（合成是从5'-3'，对 RNA 模版而言是从其3'-5'），会在合成的cDNA 3' 末端添加2-5个untemplated nucleotides，而模版转换反应（template-switching reaction）正依赖于这几个碱基。TSO（template-switching oligos）3'末端会携带与这几个碱基互补的序列，与RNA 5'端连接，然后反转录酶实现模版转换，在cDNA 3'端合成与TSO互补的序列。  
3. 预扩增。当合成第一条链后，进行有限次数的循环，得到的材料足够进行后续分析即可。  
4. Tagmentation。对扩增的 cDNA 快速高效地构建测序文库，DNA 片段化和接头连接在同一步完成，使用 Illumina 双接头策略。  
5. 扩增并测序。对带接头的 DNA 进行扩增并测序。  
  
与SMART-seq的差异：  
1. SMART-seq是在42度下进行RT，但由于部分RNA具有二级结构，出现位阻现象，导致链的延伸提前终止。SMART-seq2中加入了甜菜碱，是一种甲基供体，能够增加蛋白质热稳定性。SMART-seq2先在50度下反应2 min，破坏RNA二级结构，然后回到42度，进行2 min RT，循环十次，提高cDNA产量  
2. 甜菜碱的加入需要提升氯化镁浓度。镁离子能与甜菜碱中的羧酸阴离子结合，作为缓冲液帮助细胞应对渗透压的变化。有报道称过量氯化镁可能对PCR保真性有副作用，但在SMART-seq2的RT-PCR以及PCR阶段并未发现  
3. SMART-seq使用的TSO 3'端是3个riboguanosines
  
# CEL-seq
**CEL-seq** 通过体外扩增实现单细胞测序。  
![CEL-seq](/img/2018-08-04-scRNA-seq-technology/CEL-seq.png)  

# CITE-seq  
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
