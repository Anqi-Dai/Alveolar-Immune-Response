# this script will do a differential expression using limma

library(tidyverse)
library(limma)
library(Biobase)

# load the data (filtered)
fdat <- read_rds('processedData/bewley_data_healthy_subjects_rma_FILTERED.RDS')
pData(fdat)$Infect <- as.character(pData(fdat)$Infect)


# Limma
infect_status <- pData(fdat)$Infect
design <- model.matrix(~ 0 + infect_status)
colnames(design) <- c('CGSP14','none')
fit <- lmFit(fdat, design)

command_str <- paste("makeContrasts(","(", 'CGSP14' , "-", 'none', ")", ",levels = design)", sep = "")
contrast.matrix <- eval(parse(text =command_str)) 

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# organize the limma result
limma.res <- data.frame(
  t_stats = fit2$t %>% as.data.frame,
  pval = fit2$p.value %>% as.data.frame
) %>%
  rename(t_stats = names(.)[1],
         pval= names(.)[2]) %>%
  rownames_to_column %>%
  rename(probesetID = rowname) %>%
  mutate(padj = p.adjust(pval, method = 'BH')) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

