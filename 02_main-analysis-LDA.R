
# script info -------------------------------------------------------------
## Program name : main-analysis-LDA
## Name of programmer(s): Yiwen Zhu
## Dates : created on May 26, 2021
## Purpose of the program: 
##  To perform the main analysis of paper: gene-set-level linear discriminant analysis, 
##  Using unstandardized PCs of promoter CpGs (FEM annotations)
## Input files: 
##  un-standardized gene-level PCs: "results/2021-03-26_pca_promoter_unstd-FEM.Rdata"
##  ARIES phenotype data: [file path for phenotype data]
##  Financial hardship variables: "data/20210526_Fscore-fixed.Rdata"
##  list of sensitive period genes: "data/sensitive_period_genes.xlsx"
## Output files:
##  results for Table 2: "results/20210528_main-analysis-Table1-LDA-boot.csv"

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
## each row corresponds to a participant
## columns include all exposure variables and covariates
load("[File path for phenotype data]")
df <- aries[, 1:102] ## extract phenotype variables
rm(aries)

## list of sensitive period genes
senspgenes <- read_excel("data/sensitive_period_genes.xlsx")

## merging the PCs (same ordering of rows)
df <- cbind(df, genes_PCs_unstd)

## load financial hardship variables
load("data/20210526_Fscore-fixed.Rdata")
df <- left_join(df, fscore_fixed, by = "ID")


# code binned exposure variables ------------------------------------------

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

## create sensitive period age groups 
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
  
  ## we code individuals who had ANY exposure during that time period to be exposed
  ## but only individuals with complete data all 0s to be unexposed
  
  ## drop = FALSE retains matrix format even if only one column is selected
  ## If ‘TRUE’ (default) the result is coerced to the lowest possible dimension
  ## this is needed to accommodate settings when only one time point is available
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

## check frequencies 
lapply(df[, grep("_bin_", colnames(df), value = T)], table)
# save(df, file = "data/2021-05-31_analytic-sample.Rdata")

# Linear discriminant analysis--------------------------------------------------

## remove genes that do not have any promoter CpG in this data
included_genes <- unique(gsub("_pc.score.*", "", colnames(genes_PCs_unstd)))
senspgenes <- senspgenes %>% 
  filter(NCBI.Gene.ID %in% included_genes)

LDA_out <- lapply(advers, function(adver){
  ## variables that include "_bin_" are recoded exposure variables 
  ## corresponding to the three sensitive period hypotheses
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
  
  # compile LDA results from all three gene sets 
  results_3sets <- lapply(c("opening", "closing", "expression"), function(set_name){
    ## extract gene names (reformat)
    gene_names <- senspgenes %>% 
      filter(`Functional group for analysis` == set_name) %>% 
      .$NCBI.Gene.ID %>% 
      unique(.) %>% 
      gsub("-", "_", .)
    pc_names <- c(paste0(gene_names, "_pc.score1"), paste0(gene_names, "_pc.score2"))
    # PILRB only had one cpg; need to take out the 2nd pc name
     pc_names <- pc_names[pc_names %in% colnames(genes_PCs_unstd)]

    adver.r <- do.call("rbind", lapply(adver.names, function(hypo){
      # regression model: regressing PCs in each gene set 
      # on a sensitive period hypothesis
      f <- as.formula(paste0("cbind(", paste(pc_names, 
                                             collapse = ","), ") ~", hypo))
      mdl <- lm(f, data = df1)
      # linear discriminant analysis, with the sensitive period hypothesis
      # being the binary indicator variable in the model
      lda.out <- candisc(mdl, term = hypo)
      # extract canonical R-squared 
      lda.canrsq <- as.numeric(as.character(lda.out$canrsq))
      
      # extract other statistics
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
      # select the hypothesis with the maximal canonical R-squared
      filter(lda.canrsq == max(.$lda.canrsq)) %>% 
      mutate(geneset = set_name) %>% 
      mutate(adversity = adver)

    return(selected)
  }) %>% bind_rows() 
  return(results_3sets)
} %>% bind_rows()
)

LDA_out <- LDA_out %>% bind_rows()

# expected output format: 
# Columns: hypothesis selected, LDA canonical R-squared, P-value, 
# Wilks' lambda, gene set name, adversity name
write.csv(LDA_out, file = "results/20210526_main-analysis-Table1-LDA.csv")


# bootstrapping p-value ---------------------------------------------------

# function that calculates largest canonical R-squared under null hypothesis
# i.e., because data are scrambled, 
# there should be no association between exposure and DNAm of each gene set
calcNull <- function(Data, indices, n, npred) {
  # randomly sample exposures 
  X_hypos <- Data[sample(indices[1:n]),1:npred]
  y <- Data[indices[1:n],(npred+1):ncol(Data)]
  y <- as.matrix(y)
  
  ## get the max canonical R-squared after running LDA with all hypos 
  max <- max(sapply(X_hypos, function(x){
    mdl <- lm(y~x)
    lda.out <- candisc(mdl)
    lda.canrsq <- as.numeric(as.character(lda.out$canrsq))
  }))
  return(max)
}

sets <- c("opening", "closing", "expression")

# initialize
lda.boot <- lapply(1:7, function(x) vector("list", length = length(sets)))

# number of bootstrap samples
R <- 5000

for (i in 1:length(sets)) {
  set_name <- sets[i]
  print(set_name)
  for (j in 1:length(advers)) {
    adver <- advers[j]
    adver.names <- grep(paste0("^", adver, "_bin_"), colnames(df), value = T)
    print(adver)
    
    ## get residuals
    df1 <- df[, c("ID", adver.names, covars, 
                  grep("_pc", colnames(df), value = T))]
    df1 <- df1[complete.cases(df1), ]
    colnames(df1) <- gsub("-", "_", colnames(df1))
    
    df1[, grep("_pc", colnames(df1), value = T)] <- lapply(
      grep("_pc", colnames(df1), value = T), function(x){
        lm(as.formula(paste0(x, "~", 
                             paste(covars, collapse = "+"))), 
           data = df1)$residuals
      })
    
    ## construct X hypothesis matrix
    X_hypos <- df1 %>% 
      dplyr::select(one_of(adver.names))
    
    ## construct y matrix (residuals of PCs)
    gene_names <- senspgenes %>% 
      filter(`Functional group for analysis` == set_name) %>% 
      .$NCBI.Gene.ID %>% 
      unique(.) %>% 
      gsub("-", "_", .)
    pc_names <- c(paste0(gene_names, "_pc.score1"), paste0(gene_names, "_pc.score2"))
    # PILRB only had one cpg
    pc_names <- pc_names[pc_names %in% colnames(genes_PCs_unstd)]    
    y <- as.matrix(df1[, pc_names])
    n <- nrow(y)
    npred <- length(adver.names)
    
    boot.data <- cbind(X_hypos, y)
    
    # Run bootstrap
    set.seed(i*j)
    bsresults <- boot(data = boot.data, statistic = calcNull, R = R, n = n, npred = npred)
    
    # Find observed value of 
    t_obs <- LDA_out %>% filter(geneset == set_name, adversity == adver) %>% .$lda.canrsq
    
    # # plot histogram for reassurance
    # hist(bsresults$t, breaks=100, xlab="",
    #      main="Bootstrap null distribution of largest canonical R-squared",
    #      sub=paste(c("Observed value was",t_obs), collapse=" "))
    # abline(v=t_obs, lwd=2, col = "red")
    # 
    # Calculate p-value 
    ## (proportion of null tests with value more extreme than the observed)
    boot.pval <- sum(abs(bsresults$t)>=abs(t_obs))/R
    # boot.ci.95 <- boot.ci(bsresults, type = "perc")
    lda.boot[[j]][[i]] <- c(adver, set_name, boot.pval)
  }
}

lda.boot <- do.call("rbind", lapply(lda.boot, function(x) as.data.frame(do.call("rbind", x))))
colnames(lda.boot) <- c("adversity", "gene_set", "boot.pval")
lda.boot$boot.pval <- as.numeric(as.character(lda.boot$boot.pval))
colnames(lda.boot)[2] <- "geneset"

LDA_out <- left_join(LDA_out, lda.boot, by = c("adversity", "geneset"))

write.csv(LDA_out, file = "results/20210528_main-analysis-Table1-LDA-boot.csv", 
          row.names = F, quote = F)

