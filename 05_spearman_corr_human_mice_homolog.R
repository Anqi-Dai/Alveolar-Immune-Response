# this script will do:
# Spearman correlation of homolog-mapped log2 fold changes for HDE vs EDE and HDE vs NDE, 
# one scatter plot of log2 fold change from each comparison(HDE vs EDE and HDE vs NDE)
library(tidyverse)
library(ggpubr)
library(AnnotationDbi)
library(hgu133plus2.db)
library(biomaRt)
library(VennDiagram)

########################################################################
## loading data
# load the Jay data result that contains HDE and EDE log2FC (copies and organized from the original Gower result)
mice <- read_csv('JayData/Jay_data_core_result.csv')

# load the Bewley data limma result that contains log2FC
human <- read_csv('output/differential_expression_statistics_healthy_samples_Bewley_ALL.csv')

########################################################################
# for HUMAN
# stick to the entrez ID mapping

# get the human ENTREZID for the human data probesetID
human_feature <- AnnotationDbi::select(hgu133plus2.db, human$probesetID, c('ENTREZID')) %>%
  filter(!duplicated(PROBEID)) %>%
  dplyr::select(probesetID =PROBEID,  Human_Entrez_ID = ENTREZID )

human_table <- human %>%
  left_join(human_feature, by = 'probesetID')  %>%
  group_by(Human_Entrez_ID) %>%
  summarise_all(funs(.[which.max(abs(logFC))]))  %>%
  dplyr::select(Human_Entrez_ID, human_FC = logFC, human_FDR = adj.P.Val, human_t = t) %>%
  mutate(Human_Entrez_ID = as.numeric(Human_Entrez_ID))

########################################################################
# for Mouse 

multiple_mouse_id_genes <- mice %>%
  filter(!is.na(Human_Entrez_ID)) %>%
  group_by(Human_Entrez_ID) %>%
  summarise(mouseID_num = n()) %>%
  arrange(desc(mouseID_num)) %>%
  filter(mouseID_num > 1)

# excluding the above
mice_table <-  mice %>%
  filter(!is.na(Human_Entrez_ID)) %>%
  filter(! Human_Entrez_ID %in% multiple_mouse_id_genes$Human_Entrez_ID)

# now the mouse and human entrez id should have a one-to-one relationship

########################################################################
# join the mouse and human data together

final <- mice_table %>%
  left_join(human_table , by ='Human_Entrez_ID')%>%
  filter(!is.na(human_FC)) 


write_csv(final, 'output/fold_change_HDE_NDE_EDE.csv')
