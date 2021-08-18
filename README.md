## Anti-viral CD8 T cell developmental trajectory inference with single cell RNA-seq

> Huitian Diao (Yolanda)<br/>
> Matthew Pipkin

- *Manuscript & Figure Update*: 
  - Gsuite:
    - [Figures](https://drive.google.com/drive/folders/1a3ZAHna3BRwikHiOBOHLk-nxJNz5t076?usp=sharing)
    - [Manuscript](https://drive.google.com/drive/folders/1H5hCcD-f-DKY9Gsr309r0NCCHH4vL4h3?usp=sharing)
    - [Slides](https://drive.google.com/drive/folders/1-FMAa7kR0APNTA1pJ0ZcNdUFxjSL5J0K?usp=sharing)
    - [Sup tables](https://drive.google.com/drive/folders/1ffSiZvIphrq_wSHkYEGpZHQfreT8Mb_V?usp=sharing)
  - [Dropbox](https://www.dropbox.com/sh/lrswxf2msgenqcj/AADE3R-FuQcxOk59wkrtzQ5Ja?dl=0)
- *Original sequencing data*: 
  - [pipkinngs/Exp391_Acute-Chronic_SC]
  - []
  - []
- *Linked flow data*: [FCS files](https://drive.google.com/open?id=1-dELlhTREXr1Opehsvfqiok8N6dYf9OD)
- *Experiment design*:

- *Related Repos*
    - [Exp391 Acute-Chronic scRNAseq](https://github.com/Yolanda-HT/Exp391_Acute-Chronic_SC)
    - [Exp392 shCRF scRNAseq](https://github.com/Yolanda-HT/Exp392_shCRF_SC)
    - [Exp334 CD25KO scRNAseq](https://github.com/Yolanda-HT/Exp334CD25KOSc)
    - [Exp337 CD25KO Nascent RNAseq](https://github.com/Yolanda-HT/Exp337CD25KONascent)
    - [CRF screen](https://github.com/ScrippsPipkinLab/CRF_Screen)
    - [Exp349 shSmarca4 Bulk RNAseq](https://github.com/Yolanda-HT/Exp349_shBrg1_RNAseq)
    - [Exp276 Chd7KO Bulk RNAseq](https://github.com/Yolanda-HT/Exp276_Chd7KO_RNAseq)
    - [TFclass Database Collection](https://github.com/Yolanda-HT/TFclassDataCollection)
    - [HSAP](https://github.com/Yolanda-HT/HSAP)
    - [Scribe-py](https://github.com/Yolanda-HT/Scribe-py)

## Results
1. [Exp391: Acute v.s. Chronic infection](0_Acute-Chronic.md)
2. [Exp334: Il2raKO]()
3. [Exp392: shCRF]()

## Abstract
Individual naive CD8 T cells activated in lymphoid organs differentiate into functionally diverse and anatomically distributed T cell phylogenies in response to intracellular microbes. During infections that resolve rapidly, including live viral vaccines<sup>1</sup>, distinct effector (T<sub>EFF</sub>) and memory (T<sub>MEM</sub>) cell populations develop that ensure long term immunity<sup>2</sup>. During chronic infections, responding cells progressively become dysfunctional and “exhaust”<sup>3</sup>. A diverse taxonomy of T<sub>EFF</sub>, T<sub>MEM</sub> and exhausted (T<sub>EX</sub>) CD8 T cell populations is known, but the initial developmental basis of this phenotypic variation remains unclear<sup>4-10</sup>. Here, we defined single-cell trajectories and identified chromatin regulators that establish antiviral CD8 T cell heterogeneity using unsupervised analyses of single-cell RNA dynamics<sup>11-13</sup> and an <i>in vivo</i> RNAi screen<sup>14</sup>. Activated naive cells differentiate linearly into uncommitted effector-memory progenitor (EMP) cells, which initially branch into an analogous manifold during either acute or chronic infection. Disparate RNA velocities in single EMP cells initiate divergence into stem, circulating, and tissue-resident memory lineages that generate diverse T<sub>MEM</sub> and T<sub>EX</sub> precursor states in specific developmental orders. Interleukin-2 receptor (IL-2R) signals are essential for formation and transcriptional heterogeneity of EMP cells, and promote trajectories toward T<sub>EFF</sub> rather than T<sub>EX</sub> states. Nucleosome remodelers <i>Smarca4</i> and <i>Chd7</i> differentially promote transcription that delineates divergent T<sub>MEM</sub> lineages before cooperatively driving terminal T<sub>EFF</sub> cell differentiation. Thus, the lineage architecture is established by specific chromatin regulators that stabilize diverging transcription in uncommitted progenitors.

![Best](https://www.researchgate.net/profile/Michael_Dustin/publication/279840489/figure/fig1/AS:409936333950985@1474747847703/Gene-expression-profiles-associated-with-the-activation-and-memory-formation-of-CD8-T.png)
