library(tidyverse)
library(ggplot2)
library(dbplyr)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)

    #REST-OF-BODY PPIN VIOLIN PLOT   

#Getting the original LOEUF score file (manipulated in excel into LOEUF_score_refined.xlsx)
library(readxl)
LOEUF_input <- read_excel("LOEUF_score_refined.xlsx") %>% 
  dplyr::select(gene, oe_lof_upper)

#Level 0, 1, and 2 STRING rest-of-body PPIN  
#FILE NAME: STRING_PPIN_whole_body
STRING_PPIN_whole_body <- read_excel("STRING_PPIN_whole_body.xlsx")%>% 
  dplyr::select(gene, level)
STRING_PPIN_whole_body$level <- as.factor(STRING_PPIN_whole_body$level)
my_comparisons <- list(c("0","1"), c("1", "2"), c("0", "2")) 

#Renaming the file 
LOEUF_PPIN_whole_body <- merge(LOEUF_input, STRING_PPIN_whole_body, by="gene", how="INNER") 

#Making the violin plot of level 0, 1, and 2 rest-of-body PPIN 
ggplot(data = LOEUF_PPIN_whole_body, aes(x=level, y=oe_lof_upper, fill=level)) + 
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(2,2.1,2.2) ) +
  stat_compare_means(label.y = 2.5) +
  
  geom_boxplot(width=0.1, color="white") + 
  
  #X and Y axis labels 
  labs(x="Level of gene interaction", y="LOUEF score") + 
  
  #White background + no axis lines  
  theme_bw() + theme_classic() + 
  
  #Adding the LOUEF score cut off 
  geom_hline(yintercept = 0.35)  

##STORE 'LOEUF_PPIN_whole_body' AS A FILE 
LOEUF_PPIN_whole_body  <- data.frame(LOEUF_PPIN_whole_body) 

save(LOEUF_PPIN_whole_body , file = "df.Rdata")
load(file = "df.Rdata")
write.csv(LOEUF_PPIN_whole_body , "LOEUF_PPIN_whole_body.csv", row.names = FALSE)

#Manipulated excel file 'LOEUF_score_refined.xlsx'
#(1) add an extra column titled 'Level' labelled 'all_genes'
#(2) Removed column 'transcript'
#(3) Added 'oe_lof_upper' and 'Level' column for level 1 and 2 STRING proteins
#NEW FILE NAME: All_genes_LOEUF_PPIN_WHOLE_BODY 

library(readxl)
All_genes_LOEUF_PPIN_WHOLE_BODY  <- read_excel("All_Genes_LOEUF_PPIN_WHOLE_BODY.xlsx")
my_comparisons <- list(c("0","1"), c("1", "2"), c("0", "2"), c("0", "all_genes"), c("1", "all_genes"),c("2", "all_genes") ) 

#Making the violin plot of all genes in gnomAD, level 0, 1, and 2 rest-of-body PPIN 
ggplot(data = All_genes_LOEUF_PPIN_WHOLE_BODY, aes(x=Level, y=oe_lof_upper, fill=Level)) + 
  geom_violin() +
  
  scale_x_discrete(limits = c("all_genes", "0", "1", "2")) + 
  
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(2, 2.1, 2.2, 2.3, 2.4, 2.1) ) +
  stat_compare_means(label.y = 2.8) +
  
  geom_boxplot(width=0.1, color="white") + 
  
  #X and Y axis labels 
  labs(x="Level of gene interaction", y="LOUEF score") + 
  
  #White background + no axis lines  
  theme_bw() + theme_classic() + 
  
  #Adding the LOUEF score cut off 
  geom_hline(yintercept = 0.35)

#Kruskal Wallis with post hoc Dunn test for comparisons made between all genes in gnomAD, level 0, level 1, and 2 
#of the rest-of-body PPIN (6 comparisons in total)
kruskal.test(All_genes_LOEUF_PPIN_WHOLE_BODY$oe_lof_upper ~ All_genes_LOEUF_PPIN_WHOLE_BODY$Level)
install.packages("dunn.test")
library(dunn.test)
install.packages("FSA")
library(FSA)
dunnTest(All_genes_LOEUF_PPIN_WHOLE_BODY$oe_lof_upper,All_genes_LOEUF_PPIN_WHOLE_BODY$Level, method="bonferroni")
  

  #LUNG PPIN VIOLIN PLOT 

#Getting the original LOEUF score file (manipulated in excel to into LOEUF_score_refined.xlsx )
  library(readxl)
LOEUF_input <- read_excel("LOEUF_score_refined.xlsx") %>% 
  dplyr::select(gene, oe_lof_upper)

#Level 0, 1, and 2 STRING lung (PPIN)
#FILE NAME: STRING_PPIN_lung
STRING_PPIN_lung <- read_excel("STRING_PPIN_lung.xlsx")%>% 
  dplyr::select(gene, level)
STRING_PPIN_lung$level <- as.factor(STRING_PPIN_lung$level)
my_comparisons <- list(c("0","1"), c("1", "2"), c("0", "2")) 

#Renaming the file 
LOEUF_PPIN_lung <- merge(LOEUF_input, STRING_PPIN_lung, by="gene", how="INNER") 

#Making the violin plot of level 0, 1, and 2 lung PPIN 
ggplot(data = LOEUF_PPIN_lung, aes(x=level, y=oe_lof_upper, fill=level)) + 
  geom_violin() +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(2,2.1,2.2) ) +
  stat_compare_means(label.y = 2.5) +
  
  geom_boxplot(width=0.1, color="white") + 
  
  #X and Y axis labels 
  labs(x="Level of gene interaction", y="LOUEF score") + 
  
  #White background + no axis lines  
  theme_bw() + theme_classic() + 
  
  #Adding the LOUEF score cut off 
  geom_hline(yintercept = 0.35)  

##STORE 'LOEUF_PPIN_lung' AS A FILE 
LOEUF_PPIN_lung<- data.frame(LOEUF_PPIN_lung) 

save(LOEUF_PPIN_lung , file = "df.Rdata")
load(file = "df.Rdata")
write.csv(LOEUF_PPIN_lung, "LOEUF_PPIN_lung.csv", row.names = FALSE)

#Manipulated excel file 'LOEUF_score_refined.xlsx'
#(1) add an extra column titled 'Level' labelled 'all_genes'
#(2) Removed column 'transcript'
#(3) Added 'oe_lof_upper' and 'Level' column for level 1 and 2 STRING proteins (lung)
#NEW FILE NAME: All_Genes_LOEUF_PPIN_lung

library(readxl)
All_Genes_LOEUF_PPIN_lung <- read_excel("All_Genes_LOEUF_PPIN_lung.xlsx")
my_comparisons <- list(c("0","1"), c("1", "2"), c("0", "2"), c("0", "all_genes"), c("1", "all_genes"),c("2", "all_genes") ) 

#Making the violin plot of all genes in gnomAD, level 0, 1, and 2 lung PPIN 
ggplot(data = All_Genes_LOEUF_PPIN_lung, aes(x=Level, y=oe_lof_upper, fill=Level)) + 
  geom_violin() +
  
  scale_x_discrete(limits = c("all_genes", "0", "1", "2")) + 
  
  stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(2, 2.1, 2.2, 2.3, 2.4, 2.1) ) +
  stat_compare_means(label.y = 2.8) +
  
  geom_boxplot(width=0.1, color="white") + 
  
  #X and Y axis labels 
  labs(x="Level of gene interaction", y="LOUEF score") + 
  
  #White background + no axis lines  
  theme_bw() + theme_classic() + 
  
  #Adding the LOUEF score cut off 
  geom_hline(yintercept = 0.35)

#Kruskal Wallis with post hoc Dunn test for comparisons made between all genes in gnomAD, level 0, level 1, and 2 
#of the lung PPIN (6 comparisons in total)
kruskal.test(All_Genes_LOEUF_PPIN_lung$oe_lof_upper ~ All_Genes_LOEUF_PPIN_lung$Level)
install.packages("dunn.test")
library(dunn.test)
install.packages("FSA")
library(FSA)
dunnTest(All_Genes_LOEUF_PPIN_lung$oe_lof_upper,All_Genes_LOEUF_PPIN_lung$Level, method="bonferroni")

  #Test for significance between LOEUF scores for level 1 lung PPIN 
  #and level 1 rest-of-body PPIN (Wilcoxon signed-rank test)
wilcox.test(Level_1_LOEUF_lung$oe_lof_upper, Level_1_LOEUF_whole_body$oe_lof_upper)

  #Test for significance between LOEUF scores for level 2 lung PPIN 
  #and level 2 rest-of-body PPIN (Wilcoxon signed-rank test)
wilcox.test(Level_2_LOEUF_lung$oe_lof_upper, Level_2_LOEUF_whole_body$oe_lof_upper)
