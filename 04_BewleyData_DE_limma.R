# this script will do a differential expression using limma
library(tidyverse)
library(limma)
library(Biobase)

##########################################################################
# load the data (filtered)
fdat <- read_rds('processedData/bewley_data_healthy_subjects_rma_FILTERED.RDS')
pData(fdat)$Infect <- as.character(pData(fdat)$Infect)

##########################################################################
# run Limma
infect_status <- pData(fdat)$Infect
design <- model.matrix(~ 0 + infect_status)
colnames(design) <- c('CGSP14','none')
fit <- lmFit(fdat, design)

command_str <- paste("makeContrasts(","(", 'CGSP14' , "-", 'none', ")", ",levels = design)", sep = "")
contrast.matrix <- eval(parse(text =command_str)) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

limma1 <- topTable(fit = fit2, coef = 1,adjust="BH", number = nrow(fit2)) %>%
  rownames_to_column 

##########################################################################
# output the limma result  
colnames(limma1)[1] <- 'probesetID'

limma1 %>%
  write_csv('output/differential_expression_statistics_healthy_samples_Bewley_ALL.csv')



limma1  %>%
  filter(adj.P.Val < 0.05) %>%
  write_csv('output/differential_expression_statistics_healthy_samples_Bewley_SIG.csv')
