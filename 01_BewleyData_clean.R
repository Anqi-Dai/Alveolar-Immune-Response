# this script will work on QC of the Bewley *raw* data downloaded from https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6491/, the AM expression from healthy human subjects --- 3 exposed vs 3 control according to the original manuscript's methods. 

# Installing and loading required packages ========================================================
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("affyPLM", "affy", "AnnotationDbi","hgu133plus2.db"), version = "3.8")

library(affyPLM)
library(affy)
library(tidyverse)
library(limma)
library(annotate)
library(Biobase)
library(AnnotationDbi)
library(hgu133plus2.db)


# Make an expression set object with the healthy human subjects

## read those cel files in the Source column (AffyBatch extends eset class)

data_path <- 'BewleyData/E-MTAB-6491.raw.1'
fns <- list.files(data_path, pattern = '^[1-6]M|^[1-6]D', full.names = T)

# use readAffy function to read in the data
# and the resulting object is just like an expressionset with 6 columns and XXX rows with rows being the probesets.
dat <- ReadAffy(filenames = fns, verbose=T)

# re-assign the feature info

feature_info <- AnnotationDbi::select(hgu133plus2.db, featureNames(dat), c("SYMBOL","ENTREZID", "GENENAME")) %>%
  filter(!duplicated(PROBEID))%>%
  column_to_rownames(var = 'PROBEID')  

featureData <- new('AnnotatedDataFrame', data = feature_info)

featureData(dat) <- featureData

# re-assign the pheno info

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
  mutate(Infect = factor(Infect, levels = c('CGSP14','none'))) %>%
  dplyr::select(-Disease)

pheno <- pheno %>%
  mutate(rowname = stringr::str_replace_all(Source, '.CEL$', '') ) %>%
  mutate(rowname = paste0('S',rowname) ) %>%
  arrange(Source) %>%
  column_to_rownames(var = 'rowname')  


phenoData <- new("AnnotatedDataFrame", data = pheno)

# change the sample names in the dat
sampleNames(dat) <- rownames(pheno)

# re-assign
phenoData(dat) <- phenoData

# double check the current dat (affyBatch obj)
head(pData(dat))
head(fData(dat))


## save the eset
write_rds(dat, 'processedData/bewley_data_healthy_subjects.RDS')


