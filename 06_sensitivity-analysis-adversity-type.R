
# script info -------------------------------------------------------------
## Program name : sensitivity-analysis-adversity-type
## Name of programmer(s): Yiwen Zhu
## Dates : created on June 1, 2021
## Purpose of the program: 
##  To perform the sensitivity analysis: gene-set-level linear discriminant analysis, 
##  Merging across *adversity* to examine the effect of adversity type
## Input files: 
##  un-standardized gene-level principal components: "results/2021-03-26_pca_promoter_unstd-FEM.Rdata"
##  ARIES phenotype data: [file path phenotype data]
##  financial hardship variables: "data/20210526_Fscore-fixed.Rdata"
##  list of sensitive period genes: "data/sensitive-period-genes.xlsx"
## Output files:
##  results for supplementary table: "results/20210601_supp-analysis-combined-adversity-TableS5-LDA.csv"

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
load("[file path phenotype data]")
df <- aries[, 1:102] ## phenotype data only
rm(aries)

## list of sensitive period genes
senspgenes <- read_excel("data/sensitive_period_genes.xlsx")

## PCA was done using the same data file
## so IDs are in the same order; can just cbind 
df <- cbind(df, genes_PCs_unstd)

## financial hardship variables
load("data/20210526_Fscore-fixed.Rdata")
df <- left_join(df, fscore_fixed, by = "ID")

# code exposure variables -----------------------------------------------

covars <- c("WHITE", "Female", "mom_birthage", 
            "ppregnum", "birth.weight",
            "sustained.smoke", "ed_momgest")

advers <- c("abuse", 
            ## fixed version
            "Fscore_fixed",  
            "oneadult", 
            "r_faminst", 
            "nbhqual",
            "parc",
            "mompsy")

## collapse age groups 
## very early: 0-2y (<36m)
## early: 3-5y (36<= age < 72)
## middle: 6-7y (72 <= age <= 84)


for (adver in advers){
  print(adver)
  adver.names <- grep(paste0("^", adver), colnames(df), value = T)
  timepoints <- gsub(".*_", "", adver.names)
  timeyear <- ifelse(grepl("y", timepoints), as.numeric(gsub("y.*", "", timepoints)), 
                     as.numeric(gsub("m.*","", timepoints))/12)
  vechild <- which(timeyear < 3)  
  echild <- which(timeyear >= 3 & timeyear < 6)
  mchild <- which(timeyear >= 6)
  df[, adver.names] <- lapply(df[, adver.names], function(x) as.numeric(as.character(x)))
  
  
  ## drop = FALSE retains matrix format even if only one column is selected
  ## If ‘TRUE’ (default) the result is coerced to the lowest possible dimension
  ## very early 
  any_expos <- ifelse(rowSums(df[, adver.names[vechild], drop = FALSE], na.rm = T) > 0, 1, 0)
  all_unexpos <- ifelse(rowSums(df[, adver.names[vechild], drop = FALSE], na.rm = F) == 0, 1, 0)
  df[, paste0(adver, "_bin_vechild")] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))
  
  ## early 
  any_expos <- ifelse(rowSums(df[, adver.names[echild], drop = FALSE], na.rm = T) > 0, 1, 0)
  all_unexpos <- ifelse(rowSums(df[, adver.names[echild], drop = FALSE], na.rm = F) == 0, 1, 0)
  df[, paste0(adver, "_bin_echild")] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))
  
  ## middle
  any_expos <- ifelse(rowSums(df[, adver.names[mchild], drop = FALSE], na.rm = T) > 0, 1, 0)
  all_unexpos <- ifelse(rowSums(df[, adver.names[mchild], drop = FALSE], na.rm = F) == 0, 1, 0)
  df[, paste0(adver, "_bin_mchild")] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))
  
}

lapply(df[, grep("bin_", colnames(df), value = T)], table)


# create merged exposure variables ----------------------------------------
## any, very early
any_expos <- ifelse(rowSums(df[, grep("_bin_vechild", colnames(df)), drop = FALSE], na.rm = T) > 0, 1, 0)
all_unexpos <- ifelse(rowSums(df[, grep("_bin_vechild", colnames(df)), drop = FALSE], na.rm = F) == 0, 1, 0)
df[, "ACEs_bin_vechild"] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))

## any, early
any_expos <- ifelse(rowSums(df[, grep("_bin_echild", colnames(df)), drop = FALSE], na.rm = T) > 0, 1, 0)
all_unexpos <- ifelse(rowSums(df[, grep("_bin_echild", colnames(df)), drop = FALSE], na.rm = F) == 0, 1, 0)
df[, "ACEs_bin_echild"] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))

## any, middle
any_expos <- ifelse(rowSums(df[, grep("_bin_mchild", colnames(df)), drop = FALSE], na.rm = T) > 0, 1, 0)
all_unexpos <- ifelse(rowSums(df[, grep("_bin_mchild", colnames(df)), drop = FALSE], na.rm = F) == 0, 1, 0)
df[, "ACEs_bin_mchild"] <- ifelse(any_expos == 1, 1, ifelse(all_unexpos == 1, 0, NA))

lapply(df[, grep("ACEs_", colnames(df), value = T)], table)

# LDA ---------------------------------------------------------------------

## remove genes that do not have any promoter CpG in this data
included_genes <- unique(gsub("_pc.score.*", "", colnames(genes_PCs_unstd)))
senspgenes <- senspgenes %>% 
  filter(NCBI.Gene.ID %in% included_genes)

adver.names <- grep(paste0("ACEs_bin_"), colnames(df), value = T)
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
  pc_names <- pc_names[pc_names %in% colnames(genes_PCs_unstd)]
  
  adver.r <- do.call("rbind", lapply(adver.names, function(hypo){
    
    f <- as.formula(paste0("cbind(", paste(pc_names, 
                                           collapse = ","), ") ~", hypo))
    mdl <- lm(f, data = df1)
    lda.out <- candisc(mdl, term = hypo)
    lda.canrsq <- as.numeric(as.character(lda.out$canrsq))
    png(paste0("plots/LDA_geneset_pc_ever/", hypo, "_", set_name, ".png"), 
        width = 15, height = 5, unit = "in", res = 150)
    print(plot(lda.out, which=c(1,2), scale=8, var.col="#777777", var.lwd=1,
               col=c("red","green","blue")))
    dev.off()
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

LDA_out <- results_3sets[-6]
write.csv(LDA_out, file = "results/20210601_supp-analysis-combined-adversity-TableS5-LDA.csv")
