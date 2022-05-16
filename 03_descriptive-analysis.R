# script info -------------------------------------------------------------
## Program name : descriptive-analysis
## Name of programmer(s): Yiwen Zhu
## Dates : created on May 31, 2021
## Purpose of the program: 
##  To perform descriptive analyses: 
## (1) compare analytic sample to ALSPAC entire sample (now Table 1)
## (2) generate prevalence - Figure S1
## (3) calculate correlations across adversity types (Figure S2, Figure S3)
## Input files: 
##  analytic sample: "data/2021-05-31_analytic-sample.Rdata"
##  full ALSPAC phenotype data file: [master data file path]
## Output files: "data/2021-05-31_analytic-sample.Rdata"
##  prevalence plot: "plots/20210531_FigureS1_prevalence_exposure.png"
##  plot for correlations within each adversity:"plots/20210531_FigureS2_correlation_within_adversity.png"
##  plot for correlations across adversity: "plots/20210531_FigureS3_correlation_across_adversity.png"
##  descriptive table: "results/20210531_TableS1-descriptive-comparison.csv"

# set up ------------------------------------------------------------------

## packages 
library(dplyr)
library(tableone)
library(ggplot2)
library(tidyr)
library(ggcorrplot)
library(psych)
library(RColorBrewer)
library(patchwork)

load("data/2021-05-31_analytic-sample.Rdata")

# exclude individuals that did not enter any analysis ---------------------

covars <- c("WHITE", "Female", "mom_birthage", 
            "ppregnum", "birth.weight",
            "sustained.smoke", "ed_momgest")
# prefix for adversity variables
advers <- c("abuse", 
            "Fscore_fixed",  
            "oneadult", 
            "r_faminst", 
            "nbhqual",
            "parc",
            "mompsy")

for (adver in advers){
  adver.names <- grep(paste0("^", adver, "_bin_"), colnames(df), value = T)
  df1 <- df[, c("ID", adver.names, covars, 
                grep("_pc", colnames(df), value = T))]
  df1 <- df1[complete.cases(df1), ]
  df[, paste0(adver, "_complete")] <- ifelse(df$ID %in% df1$ID, 1, 0)
}

df$any_complete <- ifelse(rowSums(df[, grepl("_complete", colnames(df))]) == 0, 0, 1)
table(df$any_complete, useNA = "always") # n=785
df <- df[df$any_complete == 1, ]


# prevalence of adversity -------------------------------------------------

# long format
df <- do.call("rbind", lapply(advers, function(adver){
  adver.names <- grep(paste0("^", adver, "_bin_"), colnames(df_analytic), value = T)
  df_sub <- df_analytic[, adver.names]
  df_sub$adversity <- adver
  colnames(df_sub) <- c("VE", "E", "M", "adversity")
  return(df_sub)
}))

df <- gather(df, key = "time_period", value = "exposure", 1:3)
for_plot <- df %>% 
  group_by(adversity, time_period) %>% 
  summarise(prev = mean(exposure, na.rm = T))
for_plot$adversity <- factor(for_plot$adversity)
levels(for_plot$adversity) <- c("Physical or \nsexual abuse", 
                                "Financial \nhardship",
                                "Maternal \npsychopathology",
                                "Neighborhood \ndisadvantage", 
                                "One adult in the \nhousehold",
                                "Caregiver physical \nor emotional abuse",
                                "Family \ninstability")
for_plot$time_period <- factor(for_plot$time_period, levels = c("VE", "E", "M"))

pd <- position_dodge2(0.5)
(p <- ggplot(aes(x = adversity, y = prev*100, fill = time_period, group = time_period), data = for_plot)+
  geom_bar(stat = "identity", position = pd, color = "black") + 
  scale_fill_brewer(type = "seq", palette = 1, 
                    labels = c("Very early childhood",
                              "Early childhood",
                              "Middle childhood"),
                    name = "Exposure time period") + 
  geom_text(aes(label = paste0(round(prev*100, 1), "%")), position = position_dodge(1), vjust = -0.5) + 
  ylim(0, 30) + 
  ylab("Prevalence of exposure (%)\n") + 
  xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, size = 11),
        legend.position = "bottom"))
## save figure
png("plots/20210531_FigureS1_prevalence_exposure.png", width = 10, height = 5, 
    unit = "in", res = 300)
p
dev.off()
        

# correlation between adversities -----------------------------------------

mat <- polychoric(df_analytic[, exposure_vars], delete = F)$rho

## edit labels
adver_labels <- c("Physical or \nsexual abuse", 
                  "Financial \nhardship",
                  "Maternal \npsychopathology",
                  "Neighborhood \ndisadvantage", 
                  "One adult in the \nhousehold",
                  "Caregiver physical \nor emotional abuse",
                  "Family \ninstability")


rownames(mat) <- sapply(adver_labels, function(adver) paste0(adver, 
                                                             c("\nVEC", 
                                                               "\nEC",
                                                               "\nMC")))
colnames(mat) <- rownames(mat)

## start with one plot per adversity (3x3)
cor_p_list <- lapply(adver_labels, function(adver){
  mat_adver <- mat[grepl(adver, rownames(mat)), grepl(adver, colnames(mat))]
  rownames(mat_adver) <- c("VEC", 
                           "EC",
                           "MC")
  colnames(mat_adver) <- rownames(mat_adver)
  cor_p <- ggcorrplot(mat_adver, type = "lower",
                      outline.col = "white",
                      ggtheme = ggpubr::theme_pubclean,
                      colors = brewer.pal(n = 3, name = "RdYlBu"), 
                      lab = TRUE, 
                      lab_size = 5, 
                      show.diag = TRUE, 
                      show.legend = FALSE,
                      title = adver
                      # title = gsub("\n","", adver)
                      ) + 
    xlab("") + 
    ylab("")
})

png("plots/20210531_FigureS2_correlation_within_adversity.png", 
    width = 12, height = 7, 
    unit = "in", res = 300)
wrap_plots(cor_p_list, ncol = 4)
dev.off()

## then one plot across for each time point

cor_p_list <- lapply(c("VEC", "EC", "MC"), function(tp){
  mat_tp <- mat[grepl(paste0("\n", tp), rownames(mat)), grepl(paste0("\n", tp), colnames(mat))]
  rownames(mat_tp) <- gsub(tp, "", rownames(mat_tp))
  colnames(mat_tp) <- rownames(mat_tp)
  cor_p <- ggcorrplot(mat_tp, type = "lower",
                      outline.col = "white",
                      ggtheme = ggpubr::theme_pubclean,
                      colors = brewer.pal(n = 3, name = "RdYlBu"), 
                      tl.srt = 45, 
                      lab = TRUE, 
                      lab_size = 4, 
                      show.diag = TRUE, 
                      show.legend = FALSE,
                      title = tp
                      # title = gsub("\n","", adver)
  ) + 
    xlab("") + 
    ylab("")
})

png("plots/20210531_FigureS3_correlation_across_adversity.png", 
    width = 12, height = 12, 
    unit = "in", res = 300)
wrap_plots(cor_p_list, ncol = 2)
dev.off()

# Table 1 comparison with ALSPAC full sample ------------------------------
advers <- c("abuse", 
            ## fixed version
            "Fscore_fixed",  
            "oneadult", 
            "r_faminst", 
            "nbhqual",
            "parc",
            "mompsy")

load("[file path for ALSPAC master data file]")
df_analytic <- df
df <- beast
rm(beast)

## remove continuous scores that we do not use in analysis
df <- df %>% dplyr::select(-nbhqual_21m, -nbhqual_33m, -nbhqual_61m, -nbhqual_7y)
## code exposure variables
for (adver in advers){
  print(adver)
  adver.names <- grep(paste0("^", adver), colnames(df), value = T)
  df[, adver.names] <- lapply(df[, adver.names], function(x) as.numeric(as.character(x)))
  
  ## drop = FALSE retains matrix format even if only one column is selected
  ## If ‘TRUE’ (default) the result is coerced to the lowest possible dimension
  any_expos <- ifelse(rowSums(df[, adver.names, drop = FALSE], na.rm = T) > 0, 1, 0)
  all_unexpos <- ifelse(rowSums(df[, adver.names, drop = FALSE], na.rm = F) == 0, 1, 0)
  df[, paste0(adver, "_bin_any")] <- ifelse(any_expos == 1, 1, 
                                            ifelse(all_unexpos == 1, 0, NA))
}
lapply(df[, grep("bin_any", colnames(df), value = T)], table)
exposure_vars <- sapply(advers, function(adver){
  adver.names <- grep(paste0("^", adver, "_bin_any"), colnames(df), value = T)
})

df$birth.weight <- df$birthweight
df <- df %>% dplyr::select(one_of(c("ID", covars, exposure_vars)))
df$analytic <- ifelse(df$ID %in% df_analytic$ID, "Yes", "No")

# generate descriptives in the analytic sample
# as well as p-values for tests comparing distributions in the analytic sample
# vs. in everyone else (ALSPAC participants not included in the analytic sample)
tab1 <- CreateTableOne(vars = c(covars, exposure_vars), 
                       strata = "analytic",
                       data = df, 
                       factorVars = c(covars, exposure_vars)[c(covars, exposure_vars) !="birth.weight"])
tab1_all_alspac <- CreateTableOne(vars = c(covars, exposure_vars), 
                                  # strata = "analytic",
                                  data = df, 
                                  factorVars = c(covars, exposure_vars)[c(covars, exposure_vars) !="birth.weight"])
## add overall descriptives
tab1_alspac <- print(tab1_all_alspac, noSpaces = TRUE, printToggle = FALSE, quote = FALSE, 
                     showAllLevels = TRUE) %>% as.data.frame(.)
tab1_export <- print(tab1, noSpaces = TRUE, printToggle = FALSE, quote = FALSE, 
                     showAllLevels = TRUE) %>% as.data.frame(.)
tab1_export <- cbind(tab1_alspac, tab1_export)

write.csv(tab1_export, "results/20210531_TableS1-descriptive-comparison.csv")
rm(df)
