# this script will do:
# Spearman correlation of homolog-mapped log2 fold changes for HDE vs EDE and HDE vs NDE, 
# one scatter plot of log2 fold change from each comparison
library(tidyverse)
library(ggpubr)
library(AnnotationDbi)
library(hgu133plus2.db)
library(biomaRt)

########################################################################
## loading data
# load the Jay data result that contains HDE and EDE log2FC (copies and organized from the original Gower result)
mice <- read_csv('JayData/Jay_data_core_result.csv')

# load the Bewley data limma result that contains log2FC
human <- read_csv('output/differential_expression_statistics_healthy_samples_Bewley_ALL.csv')

########################################################################
# mapping homolog
# finding homolog between the mice and human and only keeping those
# using the ensembl gene ID since it's more unique and won't have same thing for several genes

# get the ensembl_gene_id for the mice entrez id
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
mouse_entrez  <- mice$Mouse_Entrez_ID
mice_genes <- getBM(attributes=c('ensembl_gene_id','entrezgene',
                          'external_gene_name'),
             filters = 'entrezgene',
             values = mouse_entrez,
             mart = ensembl)

Mice <- mice_genes %>%
  dplyr::select(mice_ensembl_id = ensembl_gene_id, Mouse_Entrez_ID = entrezgene)

# get the human ensembl_gene_id for the human data probesetID
feature_human <- AnnotationDbi::select(hgu133plus2.db, human$probesetID, c("ENSEMBL",'ENTREZID')) %>%
  filter(!duplicated(PROBEID)) %>%
  dplyr::select(probesetID =PROBEID, Human_Ensembl_ID = ENSEMBL, Human_Entrez_ID = ENTREZID )

# get the entrezID and probesetID aligned together in a same table
human <- human %>%
  left_join(feature_human, by = 'probesetID') 


# filter to get the genes that have homologs in mice and human
mice_fil <- mice %>%
  filter(!is.na(Mouse_Entrez_ID) & !is.na(Human_Entrez_ID))


# I can't use Adam's entrez mapping information since too many same human entrez ids
mapping <- mice_fil %>%
  dplyr::select(Mouse_Entrez_ID, Human_Entrez_ID, NDE_fold_change, EDE_fold_change, NDE_FDR, EDE_FDR)


Mice <- Mice %>%
  filter(!duplicated(Mouse_Entrez_ID))

mapping <- mapping %>%
  left_join(Mice, by = 'Mouse_Entrez_ID')

# get the human ensembl gene id

mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

mice_ids <- mapping$mice_ensembl_id

# human / mouse
human_mice <- getLDS(attributes=c("ensembl_gene_id"),
       filters="ensembl_gene_id", values=mice_ids, mart=mart2,
       attributesL=c("ensembl_gene_id"), martL=mart1)

human_mice <- human_mice %>%
  dplyr::select(mice_ensembl_id = Gene.stable.ID, Human_Ensembl_ID = Gene.stable.ID.1) %>%
  filter(!duplicated(mice_ensembl_id))

human_FC <- human_mice %>%
  filter(!duplicated(Human_Ensembl_ID)) %>%
  left_join(human %>%
              dplyr::select(Human_Ensembl_ID, logFC, HDE_FDR = adj.P.Val) %>%
              filter(!duplicated(Human_Ensembl_ID)), by = 'Human_Ensembl_ID')

final <- mapping %>%
  inner_join(human_FC, by = 'mice_ensembl_id') %>%
  dplyr::select(mice_ensembl_id, Human_Ensembl_ID, HDE_fold_change = logFC,HDE_FDR, NDE_fold_change, NDE_FDR, EDE_fold_change, EDE_FDR) %>%
  filter(!is.na(HDE_fold_change))

write_csv(final, 'output/fold_change_HDE_NDE_EDE.csv')

########################################################################
# Spearman correlation of homolog-mapped log2 fold changes for HDE vs EDE and HDE vs NDE

cor.test(final$HDE_fold_change, final$EDE_fold_change, method = "spearman", exact=FALSE)
  
cor.test(final$HDE_fold_change, final$NDE_fold_change, method = "spearman", exact=FALSE)

########################################################################
# scatter plot of log2 fold change from each comparison

final %>%
  mutate(log10P_HDE = -log10(HDE_FDR) ) %>%
  ggscatter( x  ='HDE_fold_change', y = 'log10P_HDE',alpha = 0.2,
             title = 'HDE volcano plot', 
             xlab = 'log2 Fold Change',
             ylab = '-log10(Padj)')


final %>%
  mutate(log10P_EDE = -log10(EDE_FDR) ) %>%
  ggscatter( x  ='EDE_fold_change', y = 'log10P_EDE',alpha = 0.2,
             title = 'EDE volcano plot', 
             xlab = 'log2 Fold Change',
             ylab = '-log10(Padj)')

final %>%
  mutate(log10P_NDE = -log10(NDE_FDR) ) %>%
  ggscatter( x  ='NDE_fold_change', y = 'log10P_NDE',alpha = 0.2,
             title = 'NDE volcano plot', 
             xlab = 'log2 Fold Change',
             ylab = '-log10(Padj)')
