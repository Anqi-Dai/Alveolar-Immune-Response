# this script will work on QC of the Bewley *raw* data downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6491/, the AM expression from healthy human subjects --- 3 exposed vs 3 control according to the original manuscript's methods. 

date()

# Installing and loading required packages ========================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("affyPLM", "affy", "AnnotationDbi","hgu133plus2.db"), version = "3.8")

library(affyPLM)
library(affy)
library(tidyverse)
library(limma)
library(annotate)

# Make an expression set object with the healthy human subjects

## load the pheno data file
pheno <- read_tsv('BewleyData/E-MTAB-6491.sdrf.txt') %>%
  rename_all(
    funs(
        stringr::str_replace_all(., ' ', '_') %>%
        stringr::str_replace_all(., 'Factor_Value', '') %>%
        stringr::str_replace_all(., '\\[|\\]', '') 
      
  )) %>%
  dplyr::select(Source_Name, infect, disease) %>%
  rename(Source = Source_Name, Infect = infect, Disease = disease) %>%
  filter(Disease == 'normal') %>%
  mutate(Infect =  stringr::str_replace_all(Infect, 'Streptococcus pneumoniae ','')) %>%
  mutate(Infect = factor(Infect, levels = c('CGSP14','none')))


## read those cel files in the Source column
data_path <- 'BewleyData/E-MTAB-6491.raw.1'

dat <- pheno %>%
  dplyr::select(Source) %>%
  mutate(contents = map(Source, ~ ReadAffy(filenames = file.path(data_path, .)), compress = TRUE)) %>% 
  unnest()

dat[1,2]

test <- ReadAffy(filenames = file.path(data_path, '2MI.CEL'), compress = TRUE)
head(test@assayData)
(pData(test))
