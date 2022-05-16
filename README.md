# READ ME
Documentation for the paper "Examining the epigenetic mechanisms of childhood adversity and sensitive periods:
a gene set-based approach"
- Author: 	Yiwen Zhu
- Last updated:	May 16, 2022

In this folder, you can find the following information regarding our manuscript. 
## Main analysis scripts
- `01_main-analysis-PCA.R`: main analysis script 1, principal component analysis
	(1) Annotate 450K CpGs to promoters of sensitive period genes 
	(2) Perform PCA to summarize variation within each gene
    (3) Create Figure 2 (% variation explained by the first two PCs)
- `02_main-analysis-LDA.R`: main analysis script 2, linear discriminant analysis
	(1) Perform the main analysis of paper: gene-set-level linear discriminant analysis
	(2) Generate bootstrapping p-values

## Sensitivity & supplementary analysis scripts
- `03_descriptive-analysis.R`: 
	- (1) Compare analytic sample to ALSPAC entire sample (now Table 1)
  - (2) Generate prevalence - Figure S1
  - (3) Calculate correlations across adversity types (Figure S2, Figure S3)
- `04_sensitivity-analysis-LDA-genebody.R`:
	- (1) Perform the sensitivity analysis: gene-set-level linear discriminant analysis, using *gene body* CpGs only
	
  *The code to generate genebody PCs is identical to the code for promoter regions except that the GeneGroup filtering after FEM annotations is set to "5=gene body"*
- `05_sensitivity-analysis-timing.R`: 
	- (1) Perform the sensitivity analysis merging across *time periods* to examine the effect of ever vs. never exposed
- `06_sensitivity-analysis-adversity-type.R`: 
	- (1) Perform the sensitivity analysis merging across *adversity types* to examine the effect of being exposed to any adversity
- `07_sensitivity-analysis-individual-CpGs.R`:
	- (1) Perform the sensitivity analysis: CpG level results look up in Lussier et al. (2022, Epigenetics)
- `08_supplementary_blood-brain-corr.R`:
	- (1) Create an additional figure to show blood-brain correlations. Data queried from IMAGE-CpG (Braun et al., 2019)
