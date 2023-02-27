library(tidyverse)
library(ggplot2)
library(dbplyr)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)

  #COMPARING SIGNIFICANCE FOR SOMATIC LC DRIVER GENES AND DIFFERENT LEVELS (REST-OF-BODY PPIN) 

#All genes from gnomad EXCLUDING somatic genes 
#File name 'LOEUF_gene_score_gnomad_all_genes.xlsx' are LOEUF scores for all genes in gnomad 

library (readxl)
LOEUF_input_all_genes <- read_excel("LOEUF_gene_score_gnomad_all_genes.xlsx") 

#List of all 62 somatic LC driver genes 
all_somatic <- c("NRG1", "CD74", "EGFR", 
                 "ERBB2", "ERBB4", "PTPN13", 
                 "RAD21", "SMARCA4", "TPM3", "AK1", "BAP1", "BRAF", "DDR2", 
                 "DROSHA", "EML4", "EZR", "HIF1A",
                 "HIP1", "KDR", "KEAP1", "KIF5B",
                 "KRAS", "LRIG3", "MAP2K1", "NFE2L2", "NOTCH1", 
                 "PIK3CB","RB1", "RBM10", "RET", "SDC4", 
                 "STK11", "STRN", "TFG", "TP53", "TPR", "CUL3", 
                 "EED", "EPHA3","FAM135B", "LEPROTL1", "N4BP2", 
                 "RFWD3", "SIRPA", "USP44", "ALK", "CCDC6", "FGFR2", "MAP2K2", 
                 "MYCL", "NKX2-1", "PTPRT", "RB1", "ROS1", "SLC34A2", 
                 "SOX2", "CPEB3", "CSMD3", "FAM47C", "GPC5", "LEPROTL1", 
                 "MB21D2", "PTPRD", "RFWD3", "ZNF479") 

#Excluding somatic LC driver genes from all genes 
LOEUF_input_all_genes_exc_somatic <- LOEUF_input_all_genes %>% 
  dplyr::filter(!gene %in% all_somatic)

#Difference in LOEUF scores of 'all genes in gnomAD excluding 
#somatic LC driver genes' and 'somatic LC driver genes'
length(LOEUF_input_all_genes_exc_somatic$oe_lof_upper) #19642
length(LOEUF_input_all_genes$oe_lof_upper) #19704 (DIFFERENCE OF 62 LOEUF scores values, explained by 62 somatic LC driver genes)

#Isolating LOUEF scores of somatic genes 

#Getting LOUEF scores from "LOEUF_gene_score_gnomad_all_genes.xlsx"
library(readxl)
LOEUF_input <- read_excel("LOEUF_gene_score_gnomad_all_genes.xlsx") %>% 
  dplyr::select(gene, oe_lof_upper)

#level 0, 1, and 2 genes 
somatic_LC_driver_genes <- read_excel("somatic_LC_driver_genes.xlsx")

#Merging both datasets (LOEUF scores and somatic LC driver genes) and renaming the file 
all_somatic_LOEUF_scores <- merge(LOEUF_input, somatic_LC_driver_genes, by="gene", how="INNER") 

#Comparison of significance for LOEUF score values for all genes in gnomAD excluding 
#somatic LC driver genes vs somatic LC driver genes
sig_somatic <- wilcox.test(LOEUF_input_all_genes_exc_somatic$oe_lof_upper,  all_somatic_LOEUF_scores$oe_lof_upper)
sig_somatic 

#UNAVAILABLE COMPARISON FOR LEVEL 0 GENES EXCLUDING SOMATIC LC DRIVER GENES 
#AND LEVEL 0 SOMATIC LC DRIVER GENES 
#(Only 1 LOEUF score value for level 0 somatic LC driver gene)  

#Isolating level 1 rest-of-body PPIN  
Level_1_LOEUF_whole_body <- read_excel("Level_1_LOEUF_whole_body.xlsx")
length(Level_1_LOEUF_whole_body$oe_lof_upper) #64 LOEUF SCORES in level 1 rest-of-body PPIN

#Isolating level 1 somatic LC driver genes in level 1 rest-of-body PPIN  
level_1_somatic <- c("CD74", "ERBB2") #2 level 1 rest-of-body PPIN somatic LC driver genes 

#Filtering out level 1 rest-of-body PPIN somatic LC driver genes from level 1 rest-of-body PPIN 
Level_1_LOEUF_whole_body_exc_som <- Level_1_LOEUF_whole_body %>% 
  dplyr::filter(!gene %in% level_1_somatic)

length(Level_1_LOEUF_whole_body_exc_som$oe_lof_upper) #62 LOEUF scores values in level 1 rest-of-body PPIN
#excluding level 1 rest-of-body PPIN somatic LC driver genes 

#LOEUF scores for level 1 rest-of-body PPIN somatic LC driver genes 
level_1_somatic_LOEUF <- Level_1_LOEUF_whole_body %>%
  filter(gene %in% level_1_somatic)  

#Making files numeric 
Level_1_LOEUF_whole_body_exc_som$oe_lof_upper <- as.numeric(Level_1_LOEUF_whole_body_exc_som$oe_lof_upper)
level_1_somatic_LOEUF$oe_lof_upper <- as.numeric(level_1_somatic_LOEUF$oe_lof_upper)

#Comparison of significance for LOEUF scores values for level 1 rest-of-body PPIN excluding 
#level 1 rest-of-body PPIN somatic LC driver genes vs level 1 rest-of-body PPIN somatic LC driver genes
sig_level1 <- wilcox.test(Level_1_LOEUF_whole_body_exc_som$oe_lof_upper, level_1_somatic_LOEUF$oe_lof_upper)
sig_level1 

#Isolating level 2 rest-of-body PPIN 
Level_2_LOEUF_whole_body <- read_excel("Level_2_LOEUF_whole_body.xlsx")
length(Level_2_LOEUF_whole_body$oe_lof_upper) #191 LOEUF scores values in level 2 rest-of-body PPIN 

#Isolating level 2 somatic LC driver genes in level 2 rest-of-body PPIN 
level_2_somatic <- c("TP53", "FAM135B") #2 level 2 rest-of-body somatic LC driver genes  

#Filtering out level 2 rest-of-body PPIN somatic LC driver genes from level 2 rest-of-body PPIN 
Level_2_LOEUF_whole_body_exc_som <- Level_2_LOEUF_whole_body %>% 
  dplyr::filter(!gene %in% level_2_somatic)

length(Level_2_LOEUF_whole_body_exc_som$oe_lof_upper) #189 LOEUF scores values in level 2 rest-of-body PPIN
#excluding level 2 rest-of-body PPIN somatic LC driver genes 

#LOEUF score for level 2 rest-of-body PPIN somatic LC driver genes 
level_2_somatic <- Level_2_LOEUF_whole_body %>%
  filter(gene %in% level_2_somatic)  

#Making files numeric 
Level_2_LOEUF_whole_body_exc_som$oe_lof_upper <- as.numeric(Level_2_LOEUF_whole_body_exc_som$oe_lof_upper)
level_2_somatic$oe_lof_upper <- as.numeric(level_2_somatic$oe_lof_upper)

#Comparison of significance for LOEUF scores values for level 2 rest-of-body PPIN excluding 
#level 2 rest-of-body PPIN somatic LC driver genes vs level 2 rest-of-body PPIN somatic LC driver genes
sig_level2 <- wilcox.test(Level_2_LOEUF_whole_body_exc_som$oe_lof_upper, level_2_somatic$oe_lof_upper)
sig_level2 

#Visualization on violin plot 
my_comparisons <- list(c("0","1"), c("3", "4"), c("5", "6"), c("7", "8"))
library("vioplot")
vioplot(LOEUF_input_all_genes_exc_somatic$oe_lof_upper, all_somatic_LOEUF_scores$oe_lof_upper,
        level0_genes_LOEUF_exc_som$oe_lof_upper, level_0_somatic$oe_lof_upper,
        Level_1_LOEUF_whole_body_exc_som$oe_lof_upper, level_1_somatic$oe_lof_upper,
        Level_2_LOEUF_whole_body_exc_som$oe_lof_upper, level_2_somatic$oe_lof_upper,
        names=c("all genes excluding somatic LC driver genes",
                "somatic LC driver genes", 
                "level 0 genes excluding somatic LC driver genes", "level 0 somatic LC driver genes", 
                "level 1 genes excluding somatic LC driver genes", "level 1 somatic LC driver genes", 
                "level 2 genes excluding somatic LC driver genes", "level 2 somatic LC driver genes"),
        col=c(brewer.pal(3,"Dark2")[1], 
              brewer.pal(3,"Dark2")[1], 
              brewer.pal(3,"Dark2")[2], 
              brewer.pal(3,"Dark2")[2], 
              brewer.pal(3,"Dark2")[3], 
              brewer.pal(3,"Dark2")[3],
              brewer.pal(4,"Dark2")[1],
              brewer.pal(4,"Dark2")[1]),
        ylab = "LOEUF score", xlab = "level")  #Colours were edited in illustrator accordingly 


  #COMPARING SIGNIFICANCE FOR SOMATIC LC DRIVER GENES AND DIFFERENT LEVELS (LUNG PPIN) 
  
  #Isolating level 1 lung PPIN 
  library(readxl)
Level_1_LOEUF_lung <- read_excel("Level_1_LOEUF_lung.xlsx")
length(Level_1_LOEUF_lung$oe_lof_upper) #140 LOEUF score values for level 1 lung PPIN   

#Isolating level 1 somatic LC driver genes in level 1 lung PPIN 
level_1_somatic_lung <- c("EGFR", "ERBB4", "PTPN13", "SMARCA4") #4 level 1 lung PPIN somatic LC driver genes 

#Filtering out level 1 lung PPIN somatic LC driver genes from level 1 lung PPIN 
Level_1_LOEUF_lung_exc_level_1_somatic_lung <- Level_1_LOEUF_lung  %>% 
  dplyr::filter(!gene %in% level_1_somatic_lung)

length(Level_1_LOEUF_lung_exc_level_1_somatic_lung$oe_lof_upper) #136 LOEUF scores values in level 1 lung PPIN 
#excluding level 1 lung PPIN somatic LC driver genes 

#LOEUF score values for level 2 lung PPIN somatic LC driver genes 
level_1_somatic_lung_LOEUF <- Level_1_LOEUF_lung %>%
  filter(gene %in% level_1_somatic_lung)  

#Making files numeric 
Level_1_LOEUF_lung_exc_level_1_somatic_lung$oe_lof_upper <- as.numeric(Level_1_LOEUF_lung_exc_level_1_somatic_lung$oe_lof_upper)
level_1_somatic_lung_LOEUF$oe_lof_upper <- as.numeric(level_1_somatic_lung_LOEUF$oe_lof_upper)

#Comparing of significance for LOEUF scores values for level 1 lung PPIN excluding 
#level 1 lung PPIN somatic LC driver genes vs level 1 lung PPIN somatic LC driver genes
sig_level1_lung <- wilcox.test(Level_1_LOEUF_lung_exc_level_1_somatic_lung$oe_lof_upper, level_1_somatic_lung_LOEUF$oe_lof_upper)
sig_level1_lung

#Isolating level 2 lung PPIN 
Level_2_LOEUF_lung <- read_excel("Level_2_LOEUF_lung.xlsx")
length(Level_2_LOEUF_lung$oe_lof_upper) #347 LOEUF score values in level 2 lung PPIN  

#Isolating level 2 somatic LC driver genes in level 2 lung PPIN 
level_2_somatic_lung <- c("BRAF", "DROSHA", "EZR", "HIF1A", "NOTCH1", "STRN") #6 level 2 lung PPIN somatic LC driver genes 

#Filtering out level 2 lung PPIN somatic LC driver genes in level 2 lung PPIN 
Level_2_LOEUF_lung_exc_level_2_somatic_lung <- Level_2_LOEUF_lung  %>% 
  dplyr::filter(!gene %in% level_2_somatic_lung)

length(Level_2_LOEUF_lung_exc_level_2_somatic_lung$oe_lof_upper) #341 LOEUF scores values in level 2 lung PPIN 
#excluding level 2 lung PPIN somatic LC driver genes 

#LOEUF score values for level 2 lung PPIN somatic LC driver genes 
level_2_somatic_lung_LOEUF <- Level_2_LOEUF_lung %>%
  filter(gene %in% level_2_somatic_lung)  

#Making files numeric 
Level_2_LOEUF_lung_exc_level_2_somatic_lung$oe_lof_upper <- as.numeric( Level_2_LOEUF_lung_exc_level_2_somatic_lung$oe_lof_upper)
level_2_somatic_lung_LOEUF$oe_lof_upper <- as.numeric(level_2_somatic_lung_LOEUF$oe_lof_upper)

#Comparing of significance for LOEUF scores values for level 2 lung PPIN excluding 
#level 2 lung PPIN somatic LC driver genes vs level 2 lung PPIN somatic LC driver genes
sig_level2_lung <- wilcox.test(Level_2_LOEUF_lung_exc_level_2_somatic_lung$oe_lof_upper, level_2_somatic_lung_LOEUF$oe_lof_upper)
sig_level2_lung

#Visualization in violin plot
my_comparisons <- list(c("0","1"), c("3", "4"), c("5", "6"), c("7", "8"))
library("vioplot")
vioplot(LOEUF_input_all_genes_exc_somatic$oe_lof_upper, all_somatic_LOEUF_scores$oe_lof_upper,
        level0_genes_LOEUF_exc_som$oe_lof_upper, level_0_somatic$oe_lof_upper,
        Level_1_LOEUF_lung_exc_level_1_somatic_lung$oe_lof_upper, level_1_somatic_lung_LOEUF$oe_lof_upper,
        Level_2_LOEUF_lung_exc_level_2_somatic_lung$oe_lof_upper, level_2_somatic_lung_LOEUF$oe_lof_upper,
        names=c("all genes excluding somatic LC driver genes",
                "somatic LC driver genes", 
                "level 0 genes excluding somatic LC driver genes", "level 0 somatic LC driver genes", 
                "level 1 genes excluding somatic LC driver genes", "level 1 somatic LC driver genes", 
                "level 2 genes excluding somatic LC driver genes", "level 2 somatic LC driver genes"),
        col=c(brewer.pal(3,"Dark2")[1], 
              brewer.pal(3,"Dark2")[1], 
              brewer.pal(3,"Dark2")[2], 
              brewer.pal(3,"Dark2")[2], 
              brewer.pal(3,"Dark2")[3], 
              brewer.pal(3,"Dark2")[3],
              brewer.pal(4,"Dark2")[1],
              brewer.pal(4,"Dark2")[1]),
        ylab = "LOEUF score", xlab = "level")  #Colours were edited in illutrator accordingly 