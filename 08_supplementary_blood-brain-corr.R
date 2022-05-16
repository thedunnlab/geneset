# script info -------------------------------------------------------------
## Program name : supplementary_blood-brain-corr
## Name of programmer(s): Yiwen Zhu
## Dates : created on August 22, 2021
## Purpose of the program: 
##  To create an additional figure to show blood-brain correlations
##  data quried from IMAGE-CpG (Braun et al., 2019)
##  Ref: Braun PR, Han S, Hing B, et al. Genome-wide DNA methylation comparison 
##. between live human brain and peripheral tissues within individuals. 
##  Translational Psychiatry. 2019;9(1):47. doi:10.1038/s41398-019-0376-y
## Input files: 
##  promoter CpG annotations (FEM): "data/annotations_promoter_FEM.Rdata"
##  IMAGE-CpG database: downloaded from https://han-lab.org/methylation/default/imageCpG
## Output files:
##  results for supplementary table: "results/20210601_supp-analysis-combined-adversity-TableS5-LDA.csv"


# set up ------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)


# query using IMAGE-CpG ---------------------------------------------------

load("/data/annotations_promoter_FEM.Rdata")
senspgenes.promoter <- senspgenes.promoter %>%
  select(probeID, group = `Functional group for analysis`, NCBI.Gene.ID)

# Braun 2019, Translational Psychiatry
# load downloaded data from https://han-lab.org/methylation/default/imageCpG
Illumina_450K_covar <- read_table2("Illumina_450K_covar")
braun2019 <- Illumina_450K_covar %>% 
  filter(cgid %in% senspgenes.promoter$probeID) %>% 
  select(probeID = cgid, cor = rho.br.b)
  
braun2019 <- left_join(braun2019, senspgenes.promoter, by = "probeID")
braun2019$group <- factor(braun2019$group, levels = c("opening", "closing", "expression"))

ggplot(aes(x = cor, fill = group), data = braun2019) + 
  geom_histogram() +
  facet_wrap(~group) + 
  ylab("") + 
  xlab("Correlation quried using IMAGE-CpG (Braun et al., 2019)") + 
  ggtitle("") + 
  theme_bw() + 
  # theme(text = element_text(family = "Gill Sans")) + 
  scale_fill_brewer(type = "qual", palette = 3)

ggsave(filename = "plots/20210822_blood-brain-cor_Braun2019.png", 
       height = 3,
       width = 6, 
       units = 'in', 
       dpi = 300)

sum(braun2019$cor > 0)/nrow(braun2019)
sum(braun2019$cor >= 0.6)/nrow(braun2019)

braun2019 %>% 
  group_by(group) %>% 
  summarise(n_strong_r = sum(cor >= 0.6),
            perc_strong_r = sum(cor >= 0.6)/n(),
            n_pos_r = sum(cor > 0),
            perc_pos_r = sum(cor > 0)/n())

strong_genes <- braun2019 %>% 
  filter(cor >= 0.6) %>% 
  .$NCBI.Gene.ID %>% 
  unique(.)
length(strong_genes)/length(unique(braun2019$NCBI.Gene.ID))
