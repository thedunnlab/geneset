# script info -------------------------------------------------------------
## Program name : main-analysis-PCA
## Name of programmer(s): Yiwen Zhu
## Dates : created on March 23, 2021; cleaned up on May 31, 2021
## Purpose of the program: 
##  To annotate 450K CpGs to promoters of sensitive period genes
##  Perform PCA to summarize variation within each gene
##  And create Figure 2 (% variation explained by the first two PCs)
## Input files: 
##  list of sensitive period genes: "data/Animal Model Genes-2019-07-11.xlsx"
##  winsorized DNAm data at age 7: [file path for DNAm data]

## Output files:
##  list of promoter CpGs: "data/annotations_promoter_FEM.Rdata"
##  promoter only residuals: "data/2021-03-26_senspgenes-betas-promoter-resid-FEM.Rdata"
##  unstandardized PC scores: "results/2021-03-26_PCA_promoter_var_unstd-FEM.Rdata"
##  Figure 2: variance explained by PCs

# set up ------------------------------------------------------------------


## packages 
library(dplyr)
library(readxl)
library(ggplot2)
library(tibble)
library(tidyr)
library(ggpubr)
library(FEM)

## load list of sensitive period genes
## updated list and grouping 
## contains the following columns: gene, EntrezID, NCBI gene name, group/gene set
senspgenes <- read_excel("data/sensitive_period_genes.xlsx")
## load DNAm data (winsorized, after QC)
## each row represents a participant, and each column is DNAm for a CpG site
load("[file path for DNAm data]")


# subset to CpG sites in sensitive period genes ---------------------------

## for annotation: using the bioconductor package FEM
## Reference: A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control. Jiao Y, Widschwendter M, Teschendorff AE. Bioinformatics. 2014;30(16):2360-2366
## [https://rdrr.io/bioc/FEM/man/probeEPICfemanno.html] 

data(probe450kfemanno) # extract annotations from the FEM package for Illumina 450k 
probe450kfemanno <- as.data.frame(probe450kfemanno)
## subset probe450kfemanno to the genes we are interested in, by EntrezID
senspgenes.annot <- probe450kfemanno[probe450kfemanno$eid %in% senspgenes$EntrezID,]

## define promoter areas:TSS1500, TSS200, and 5’-UTR
## Reference: Lokk K, et al. (2014) DNA methylome profiling of human tissues identifies global and tissue-specific methylation patterns. Genome Biol 15(4):r54.
## [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4053947/]
## FEM documentation: Regions are annotated as follows: 
## 1=TSS1500, 2=TSS200, 3=5'UTR, 4=1stExon, 5=gene body, 6=3'UTR. 
## Probes with ambiguous mappings are assigned an NA.
senspgenes.promoter <- senspgenes.annot %>% 
  filter(GeneGroup %in% c(1,2,3)) %>% 
  unique(.) # n=567

## merge alias back 
senspgenes.promoter$eid <- as.character(senspgenes.promoter$eid)
senspgenes$EntrezID <- as.character(senspgenes$EntrezID)
senspgenes.promoter <- left_join(senspgenes.promoter, senspgenes, by = c("eid" = "EntrezID"))
# restrict to CpGs available in our dataset 
senspgenes.promoter <- senspgenes.promoter %>% filter(probeID %in% colnames(aries))
length(unique(senspgenes.promoter$probeID)) # n=530 CpGs

## save CpGs in the promoter, which we will use for the main analysis
# save(senspgenes.promoter, file = "data/annotations_promoter_FEM.Rdata") 
# get a list of CpG probe names 
promoter_cpgs <- as.character(unique(senspgenes.promoter$probeID)) 

## check if all sensitive period genes are represented
sum(!senspgenes$EntrezID %in% senspgenes.promoter$eid) ## 4 not on the list 
senspgenes$NCBI.Gene.ID[!senspgenes$EntrezID %in% senspgenes.promoter$eid]
# "ADAMTS4" "MMP8" "DLG4"

## non-cpg columns: phenotype variables (the first 102 columns)
non_cpg_cols <- colnames(aries)[1:102]
aries <- aries[, colnames(aries) %in% c(non_cpg_cols, 
                                        as.character(senspgenes.promoter$probeID))]
## double checking that we did include 530 CpGs in the dataset
ncol(aries) - 102 
# 530 (number of columns after subtracting the number of phenotype variables)

# regress betas on technical vars -----------------------------------------
cell_mat <- model.matrix(~ Bcell + CD4T + CD8T + Gran + Mono + NK + sample_type, 
                         data = aries)[, -1]
betas <- lm(as.matrix(aries[, 103:ncol(aries)]) ~ cell_mat)$resid
rm(aries)
## save the CpG level residuals
# save(betas, file = "data/2021-03-26_senspgenes-betas-promoter-resid-FEM.Rdata")

# summarize each gene into 2 PCs ------------------------------------------

# load("data/2021-03-26_senspgenes-betas-promoter-resid-FEM.Rdata")
# for the main analysis, restricting to promoters only
senspgenes <- senspgenes.promoter
genes <- unique(senspgenes$NCBI.Gene.ID) # 58 genes
senspgenes$probe <- as.character(senspgenes$probeID)
# a "dictionary" for groups - to look up gene set labels
senspgenes_group <- senspgenes %>% 
  dplyr::select(gene = NCBI.Gene.ID, group = `Functional group for analysis`) %>% 
  unique(.)


## PCA functions
# runGenePCA: running PCA for each gene, extracting the first two PCs 
runGenePCA <- function(standardize = FALSE){
  genes_PCs <- do.call("cbind", 
                       lapply(genes, function(gene){
                         print(gene)
                         # filter to CpGs annotated to each gene
                         gene_cpgs <- senspgenes %>% 
                           filter(NCBI.Gene.ID == gene) %>% 
                           .$probe %>% 
                           unique(.)
                         gene_data <- betas[, gene_cpgs]
                         ## PCA (default is unstandardized PCA)
                         pca_gene <- prcomp(gene_data, scale = standardize)
                         ## extract the first PC score
                         pc.score1 <- pca_gene$x[,"PC1"]
                         # for most genes, 2 PCs are extracted
                         # dealing with the exception where only 1 CpG is available
                         if (length(gene_cpgs) >= 2){
                           pc.score2 <- pca_gene$x[,"PC2"]
                           res <- cbind(pc.score1, pc.score2)}
                         else {
                           res <- cbind(pc.score1)}
                         colnames(res) <- paste0(gene, "_", colnames(res))
                         return(res)
                       }))
  return(genes_PCs)
}

## summPCvar: function to calculate variance explained, for plotting
summPCvar <- function(standardize = FALSE){
  genes_PC_summ <- do.call("rbind", lapply(genes, function(gene){
    print(gene)
    gene_cpgs <- senspgenes %>% 
      filter(NCBI.Gene.ID == gene) %>% 
      .$probe %>% 
      unique(.)
    gene_data <- betas[, gene_cpgs]
    pca_gene <- prcomp(gene_data, scale = standardize)
    summ <- summary(pca_gene)
    # extract the % of variance explained by the first PC
    PC1 <- summ$importance[2,1]
    if (length(gene_cpgs) >= 2){
    # extract the % of variance explained by the second PC
      PC2 <- summ$importance[2,2]
    } else {
      PC2 <- 0 ## no data; because there was only 1 CpG
    }
    r <- cbind(gene = gene, PC1 = PC1, PC2 = PC2, numCpGs = length(gene_cpgs))
    return(r)
  })) %>% as_tibble()
  
  genes_PC_summ <- genes_PC_summ %>% 
    gather(key = "PC", value = "variance_explained", 2:3) 
  ## add geneset labels
  genes_PC_summ <- left_join(genes_PC_summ, senspgenes_group, by = "gene")
  genes_PC_summ <- genes_PC_summ %>% 
    mutate(label = paste0(gene, " (", numCpGs, ")"))
  genes_PC_summ$variance_explained <- 
    as.numeric(as.character(genes_PC_summ$variance_explained))
  
  return(genes_PC_summ)
}


# plot variance explained ------------------------------------------------------

genes_PC_summ_unstd <- summPCvar(standardize = FALSE)
genes_PC_summ_unstd$group <- factor(genes_PC_summ_unstd$group, 
                                    levels = c("opening",
                                               "expression",
                                               "closing"))

# png("plots/2021-03-26_barplot-variancePCA-allgenes-promoter-unstd-FEM.png", height = 7,
#     width = 19, unit = "in", res = 150)
(p <- ggplot(aes(x = label, y = variance_explained*100), data = genes_PC_summ_unstd) +
    geom_bar(aes(fill = PC), stat = "identity") + 
    # geom_text(aes(label = round(variance_explained, 2)),
    # position = position_stack(vjust = 0.6), size = 3) + 
    # geom_text(data = genes_PC_summ, aes(x = label, y = variance_explained, 
    # label = round(variance_explained,2)), 
    # colour="black", show.legend = F) + 
    coord_flip() + 
    geom_hline(yintercept = 50, color = "red", linetype = "dashed") + 
    facet_wrap(~group, scales = "free_y") +
    labs(x="Gene (number of CpGs)", y="% Variance explained") +
    # ggtitle("Variance explained by the first two PCs for each gene") +
    #scale_fill_manual(values=fill) + 
    theme_pubclean() + 
    labs(fill = "") + 
    theme(text = element_text(family = "Gill Sans"), 
          axis.title = element_text(size = 18), 
          strip.text.x = element_text(size = 18), 
          plot.title = element_text(size = 18), 
          legend.title=element_text(size = 16), 
          legend.text=element_text(size = 16), 
          legend.position = "right") + 
    scale_fill_brewer(palette = 3))
# dev.off()

# save the PCs ------------------------------------------------------------
genes_PCs_unstd <- runGenePCA(standardize = FALSE)
# save(genes_PCs_unstd, file = "results/2021-03-26_pca_promoter_unstd-FEM.Rdata")
# format of saved results: each row is a sample/partcipant, 
# each column is a PC score, at most two PCs per gene 
# column names: "Gene Name_pc.score1" "Gene Name_pc.score2" 

