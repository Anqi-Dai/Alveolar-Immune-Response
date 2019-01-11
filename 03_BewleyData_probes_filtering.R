# this script will do the probes filtering according to the methods in the original paper :
# Probes whose expression was in the 20th lowest percentile were removed. 

library(tidyverse)

# load the normalized data
mdat <- readRDS('processedData/bewley_data_healthy_subjects_gcrma.RDS')

# aggregation: get the 20% quantile of each sample and select the max of that to be the threshold 
quant20 <- apply(exprs(mdat), 2, function(Sample){
  quantile(Sample,  probs = 0.2)
})

thre <- max(quant20)

# removing those probesets that have the max expression across the samples less than or equal to the thre

max_each_row <- apply(exprs(mdat), 1, max)

to_be_removed <- names(max_each_row[max_each_row <= thre])

# filter out those to_be_removed
mdat.f <- mdat[!featureNames(mdat) %in% to_be_removed,]

# write out the data
write_rds(mdat.f, 'processedData/bewley_data_healthy_subjects_rma_FILTERED.RDS')
