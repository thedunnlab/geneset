

# script info -------------------------------------------------------------
## Program name : sensitivity-analysis-individual-CpGs
## Name of programmer(s): Yiwen Zhu
## Dates : created on May 18, 2021
## Purpose of the program: 
##  To perform the sensitivity analysis: CpG level results look up in Lussier et al. (2022, Epigenetics)
##  FDR correctiin for # of CpGs in the gene set analyses (530)
## Input files: 
##  list of gene set CpGs: "data/2021-03-26_senspgenes-betas-promoter-resid-FEM.Rdata"
##  CpG level results (Lussier et al.) in the following directory: [EWAS directory] 
##  rerun Fscore results: 
## Output files:
##  full list: "2021-06-22_CpG-level-lookup_FEM.csv"

# set up ------------------------------------------------------------------
library(readr)
library(readxl)
library(dplyr)
load("data/2021-03-26_senspgenes-betas-promoter-resid-FEM.Rdata")

## look up CpG level results 
advers <- c("abuse", "Fscore",  "oneadult", 
            "r_faminst", "nbhqual","parc","mompsy")

## the full list of genes, Entrez IDs, and corresponding CpGs are available in 
## Supplementary Table 1
senspgenes_group <- read_excel("data/sensitive_period_genes.xlsx") %>% 
  dplyr::select(gene = NCBI.Gene.ID, group = `Functional group for analysis`)

amg_cpgs_slcma_promoter<- data.frame(matrix(NA, nrow = 0, ncol = 7))

for (adver in advers){
  if (adver == "Fscore"){
    load("data/SLCMA-newDNAm-rerun-sI-noFWL-Fscore_fixed_20210622.Rdata")
  } else {
    load(paste0("data/SLCMA-data-diff-06-newDNAm-educFWL-sI/06-newDNAm-educFWL-sI-", adver, "-20200415.Rdata"))
  }
  
  res <- res[res$probe %in% colnames(betas), ]
  res$adversity <- adver
  res$FDRqval <- p.adjust(res$`P-value`, method = "BH")
  amg_cpgs_slcma_promoter <- rbind(amg_cpgs_slcma_promoter, res)
}

# examining different FDR cutoffs
amg_cpgs_slcma_promoter %>% filter(`P-value` < 0.05/530)
amg_cpgs_slcma_promoter %>% filter(FDRqval < 0.05)
amg_cpgs_slcma_promoter %>% filter(FDRqval < 0.10)

## merge in gene name 
load("data/annotations_promoter_FEM.Rdata")

amg_cpgs_slcma_promoter <- left_join(amg_cpgs_slcma_promoter, senspgenes.promoter,
                                     by = c("probe" = "probeID"))
amg_cpgs_slcma_promoter[, 1:4] <- lapply(amg_cpgs_slcma_promoter[, 1:4], unlist)
amg_cpgs_slcma_promoter <- amg_cpgs_slcma_promoter %>% 
  dplyr::select(probe, adversity, `Variable_Name`, R2, `P-value`, GeneGroup,
                FDRqval, 
                NCBI.Gene.ID, `Functional group for analysis`)

# amg_cpgs_slcma_promoter$`P-value` <- round(amg_cpgs_slcma_promoter$`P-value`, 6)
amg_cpgs_slcma_promoter$FDRqval <- round(amg_cpgs_slcma_promoter$FDRqval, 3)
amg_cpgs_slcma_promoter$R2 <- round(amg_cpgs_slcma_promoter$R2, 3)
amg_cpgs_slcma_promoter %>% filter(`P-value` < 0.05)
amg_cpgs_slcma_promoter %>% filter(FDRqval < 0.10)
amg_cpgs_slcma_promoter %>% filter(FDRqval < 0.20)

write.csv(amg_cpgs_slcma_promoter, "results/2021-06-22_CpG-level-lookup_FEM.csv")
