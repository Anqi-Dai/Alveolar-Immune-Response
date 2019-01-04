# this script will do the probes filtering according to the methods in the original paper :
# Probes whose expression was in the 20th lowest percentile were removed. 

library(tidyverse)

# load the normalized data
mdat <- readRDS('processedData/bewley_data_healthy_subjects_rma.RDS')

# aggregation: get the median value of one probe expression across the samples
median.exprs <- apply(exprs(mdat), 1, median)

# the threshold of being at the 20% quantile of the array of median values
thre <- quantile(median.exprs, probs = 0.2)

# the names of the probes that pass the filter
pass <- names(median.exprs[median.exprs > thre])

# filter out the probes that the median value is at the lowest 20%
mdat.f <- mdat[featureNames(mdat) %in% pass,]

# write out the data
write_rds(mdat.f, 'processedData/bewley_data_healthy_subjects_rma_FILTERED.RDS')