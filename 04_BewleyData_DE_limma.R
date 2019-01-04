# this script will do a differential expression using limma

library(tidyverse)
library(limma)

# load the data (filtered)
fdat <- read_rds('processedData/bewley_data_healthy_subjects_rma_FILTERED.RDS')

## wrapper for limma, assumes data is already log2 normalized
run_limma <- function(eset, class_id, control, treatment)
{
  control_inds <- which(pData(eset)[, class_id] == control)
  treatment_inds <- which(pData(eset)[, class_id] == treatment)
  
  eset.control <- eset[, control_inds]
  eset.treatment <- eset[, treatment_inds]
  eset.compare <- eset[, c(control_inds, treatment_inds)]
  condition <- as.character(pData(eset.compare)[, class_id])
  colData <- data.frame(condition=as.character(pData(eset.compare)[, class_id]))
  
  design <- model.matrix(~ 0 + factor(condition))
  colnames(design) <- levels( factor(condition))
  fit <- lmFit(eset.compare, design)
  command_str <- paste("makeContrasts(",
                       "(", treatment , "-", control, ")", 
                       ",levels = design)", sep = "")
  contrast.matrix <- eval(parse(text =command_str)) 
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, coef=1, adjust="BH", number =length(fit2) ,
                  sort.by = "none")
  return(res)
}