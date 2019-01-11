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
  rownames_to_column %>%
  #rename(probesetID = rowname) %>%
  filter(adj.P.Val < 0.05)

# using the intercept and no contrast
infect_status <- pData(fdat)$Infect
infect_status <- factor(infect_status, levels= c('none','CGSP14'))
design <- model.matrix(~infect_status)
colnames(design) <- c("Intercept","CGSP14-none")
fit <- lmFit(fdat, design)

fit2 <- eBayes(fit)
limma2 <- topTable(fit = fit2, coef = 1,adjust="BH", number = nrow(fit2)) %>%
  filter(adj.P.Val < 0.05)

##########################################################################
# output the limma result  
colnames(limma1)[1] <- 'probesetID'

limma1 %>%
  write_csv('output/differential_expression_statistics_healthy_samples_Bewley.csv')
