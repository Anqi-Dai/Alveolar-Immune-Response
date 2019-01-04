# this script will do the normalization and QC of the Bewley data

library(tidyverse)
library(affyPLM)
library(ggpubr)

# loading the data (the AffyBatch object)
bdat <- read_rds('processedData/bewley_data_healthy_subjects.RDS')

####### background correction and normalization. ###########
bdat.d <- rma(bdat)

saveRDS(bdat.d, 'processedData/bewley_data_healthy_subjects_rma.RDS')

###### QC using the raw data################################

plm.fit <- fitPLM(bdat)

# compute the relative log expression (RLE) of each gene and plot the median values per sample as a barplot
median.rle <- RLE(plm.fit,main="RLE",type="stats")[1,]


median.rle.df <- data.frame(
  Sample = names(median.rle),
  med_RLE = median.rle
) 

pdf('figs/median_rle.pdf')
median.rle.df %>%
  ggbarplot(x = 'Sample', y = 'med_RLE')
dev.off()
 
# compute the normalized unscaled standard error (NUSE) of each gene and plot the median values per sample as a histogram
median.nuse <- NUSE(plm.fit,main="NUSE",type="stats")[1,]

median.nuse.df <- data.frame(
  Sample = names(median.nuse),
  med_NUSE = median.nuse
) 

pdf('figs/median_NUSE.pdf')
median.nuse.df %>%
  ggbarplot(x = 'Sample', y = 'med_NUSE')
dev.off()


# the PCA to see how the data separates (using the normalized data.)

# scaling and centering the expression matrix
m <- t(scale(t(exprs(bdat.d))))

# performing PCA
data.pca <- prcomp(m, center=FALSE)

pca.df <- data.pca$rotation

# draw the pca plot
pdf('figs/PCA_pc1VSpc2.pdf')
pca.df %>%
  as.data.frame %>%
  rownames_to_column %>%
  left_join(pData(bdat.d) %>%
              rownames_to_column, by = 'rowname') %>%
  ggscatter(x = 'PC1', y = 'PC2', color = 'Infect',
            title = 'PC2 VS PC1 of the normalized data')
dev.off()

# It appears that the raw data is of high quality.



