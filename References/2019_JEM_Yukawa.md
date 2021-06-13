# AP-1 activity induced by co-stimulation is required for chromatin opening during T cell activation
> Masashi Yukawa, Artem Barski, et. al. <br>
> Cincinnati Children's Hospital Medical Center

## Abstract
- Background
	- Naive T cells quickly open chromatin with 5 hours of activation
	- Newly opened regions strongly enriched for AP-1 motif
	- ChIP-seq: >70% new open region bound by AP-1
- Methods:
	- Broad inhibition of AP-1
		- Prevented chromatin opening at AP-1 sites and reduced exression of nearby genes
	- Induction of anergy (absenece of co-stimulation during activation)
		- Reduced induction of AP-1 and failure of proper chromatin remodeling
- Conclusion: AP-1 links T cell activation and chromatin remodeling

## Introduction
- T cell activation:
	- TCR, costimulation
		- Downstream: NFATs, AP-1 (heterodimer of Fos and Jun), NF-kB activated via Ca2+ - calcineurin, MAPK and PI3k/PKC pathways
	- Cytokines: differentiation signals
		- JAK-STAT -> lineage specific transcription factors (TFs) -> lineage specific gene expression
	- Model: Il2 locus
		- Il2 promoter: several AP-1 and NFAT binding sites (conserved between human & mice)
		- Binding sites are adjacent (AP-1 and NFAT form heteromer) -> synergize to induce Il2 expression. Mutation of these binding sites prevents Il2 expression
		- NF-kB and other TFs participate in Il2 regulation during T cell activation via binding sites near promoter
	- Different mechanisms:
		- IL2 expression: dependent on new protein synthesis
		- IL0, IFNG, TNF: not dependent on new protein synthesis

## Results:
*NOR: Naive open region* <br>
*EOR: Effector open region* <br>
*COR: Common open region*

**Human CD4 T cell**


### Characterizing open chromatin regions (Fig1, Fig2)
- Expression pattern (RNAseq) matches with chromain open region
- Comparing with NOR and COR
	- EOR has higher % in gene body
	- EOR has higher H3K4me3 and H3K27ac


### The TFs AP-1 and NFAT1 bind effector open chromatin regions (Fig3)
- Motif enrichment: EOR enriched with ATF3, JUNB, FRA2, BATF, AP-1, FRA2 etc.
- Comparing with NOR, COR: EOR has higher read density of AP-1, NFAT factors (lower cMC)

### Formation of SEs involves chromatin opening (Fig4)
*SE: super enhancer*

### AP-1 activity is required for open chromatin formation during T cell activation (Fig5)
- A Fos electroporation of T cells -> RNA-seq

### Co-stimulation is required for open chromatin formation (Fig6)
- Naive T cell activation: anti-CD3 + anti-CD28 / anti-CD3 only
	- EOR most affected when anti-CD28 is absent
	- JACK2 locus: when anti-CD28 is absent:
		- Reduced ATAC-seq open 
		- Altered H3K27ac
		- cFos, JUNB binding reduced
		- RNA-seq: expression reduced

### AP-1 sites overlap risk loci for immunological diseases (Fig7)







