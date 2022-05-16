
# script info -------------------------------------------------------------
## Program name : sensitivity-analysis-timing
## Name of programmer(s): Yiwen Zhu
## Dates : created on June 1, 2021
## Purpose of the program: 
##  To perform the sensitivity analysis: gene-set-level linear discriminant analysis, 
##  Merging across *time periods* to examine the effect of ever vs. never exposed
## 
## Input files: 
##  un-standardized gene-level principal components: "results/2021-03-26_pca_promoter_unstd-FEM.Rdata"
##  ARIES phenotype data: [file path phenotype data]
##  financial hardship variables: "data/20210526_Fscore-fixed.Rdata"
##  list of sensitive period genes: "data/sensitive-period-genes.xlsx"
## Output files:
##  results for supplementary table: "results/20210601_supp-analysis-ever-TableS4-LDA.csv"

# set up ------------------------------------------------------------------

## packages 
library(dplyr)
library(ggplot2)
library(candisc)
library(boot)
library(readxl)

## load data 
## pc scores
load("results/2021-03-26_pca_promoter_unstd-FEM.Rdata")
## ARIES phenotype data
load("[file path for DNAm data]")
df <- aries[, 1:102] ## phenotype data only
rm(aries)

## list of sensitive period genes
senspgenes <- read_excel("data/sensitive_period_genes.xlsx")

## PCA was done using the same data file
## so IDs are in the same order; can just cbind 
df <- cbind(df, genes_PCs_unstd)

## fixed financial hardship variables
load("data/20210526_Fscore-fixed.Rdata")
df <- left_join(df, fscore_fixed, by = "ID")

# code exposure variables -----------------------------------------------

covars <- c("WHITE", "Female", "mom_birthage", 
            "ppregnum", "birth.weight",
            "sustained.smoke", "ed_momgest")

advers <- c("abuse", 
            "Fscore_fixed",  
            "oneadult", 
            "r_faminst", 
            "nbhqual",
            "parc",
            "mompsy")

# creating variables for being ever exposed to each adversity before age 7
for (adver in advers){
  print(adver)
  adver.names <- grep(paste0("^", adver), colnames(df), value = T)
  df[, adver.names] <- lapply(df[, adver.names], function(x) as.numeric(as.character(x)))
  
  
  ## drop = FALSE retains matrix format even if only one column is selected
  ## If ‘TRUE’ (default) the result is coerced to the lowest possible dimension
  ## need to create two separate variables because for exposed, we require any exposure
  ## for unexposed, we require being unexposed throughout the three periods
  any_expos <- ifelse(rowSums(df[, adver.names, drop = FALSE], na.rm = T) > 0, 1, 0)
  all_unexpos <- ifelse(rowSums(df[, adver.names, drop = FALSE], na.rm = F) == 0, 1, 0)
  df[, paste0(adver, "_bin_any")] <- ifelse(any_expos == 1, 1, 
                                            ifelse(all_unexpos == 1, 0, NA))
  
}

lapply(df[, grep("bin_", colnames(df), value = T)], table)

# LDA ---------------------------------------------------------------------

## remove genes that do not have any promoter CpG in this data
included_genes <- unique(gsub("_pc.score.*", "", colnames(genes_PCs_unstd)))
senspgenes <- senspgenes %>% 
  filter(NCBI.Gene.ID %in% included_genes)

LDA_out <- lapply(advers, function(adver){
  adver.names <- grep(paste0("^", adver, "_bin_"), colnames(df), value = T)
  
  ## regressing the PCs on the covariates
  df1 <- df[, c("ID", adver.names, covars, 
                grep("_pc", colnames(df), value = T))]
  df1 <- df1[complete.cases(df1), ]
  colnames(df1) <- gsub("-", "_", colnames(df1))## because the dash was problematic
  
  ## regress each pc score on potential confounders
  ## note that technical variables including cell counts aren't included because 
  ## we already adjusted for them before PCA
  df1[, grep("_pc", colnames(df1), value = T)] <- lapply(
    grep("_pc", colnames(df1), value = T), function(x){
      lm(as.formula(paste0(x, "~", 
                           paste(covars, collapse = "+"))), 
         data = df1)$residuals
    })
  
  results_3sets <- lapply(c("opening", "closing", "expression"), function(set_name){
    
    gene_names <- senspgenes %>% 
      filter(`Functional group for analysis` == set_name) %>% 
      .$NCBI.Gene.ID %>% 
      unique(.) %>% 
      gsub("-", "_", .)
    pc_names <- c(paste0(gene_names, "_pc.score1"), paste0(gene_names, "_pc.score2"))
    # PILRB only had one cpg; need to take out the 2nd pc name
    pc_names <- pc_names[pc_names %in% colnames(genes_PCs_unstd)]
    
    adver.r <- do.call("rbind", lapply(adver.names, function(hypo){
      
      f <- as.formula(paste0("cbind(", paste(pc_names, 
                                             collapse = ","), ") ~", hypo))
      mdl <- lm(f, data = df1)
      lda.out <- candisc(mdl, term = hypo)
      lda.canrsq <- as.numeric(as.character(lda.out$canrsq))
      
      wilks.out <- Wilks(lda.out)
      wilks.pval <- as.numeric(as.character(wilks.out$`Pr(> F)`))
      wilks.lambda <- as.numeric(as.character(wilks.out$`LR test stat`))
      r <- cbind(hypo, lda.canrsq, wilks.pval, wilks.lambda)
      
      return(r)
    }))
    
    selected <- as.data.frame(adver.r) %>% 
      mutate(lda.canrsq = as.numeric(as.character(lda.canrsq)),
             wilks.pval = as.numeric(as.character(wilks.pval)),
             wilks.lambda = as.numeric(as.character(wilks.lambda))) %>% 
      filter(lda.canrsq == max(.$lda.canrsq)) %>% 
      mutate(geneset = set_name) %>% 
      mutate(adversity = adver)
    return(selected)
  }) %>% bind_rows() 
  return(results_3sets)
} %>% bind_rows()
)

LDA_out <- LDA_out %>% bind_rows()

write.csv(LDA_out, file = "results/20210601_supp-analysis-ever-TableS4-LDA.csv")
# no need to pursue bootstrapping, which was meant to correct for testing multiple hypotheses