# this script will do the normalization and QC of the Bewley data
library(tidyverse)
library(affyPLM)
library(ggpubr)
library(gcrma)
library(Biobase) 

####### background correction and normalization. ###########
data_path <- 'BewleyData/E-MTAB-6491.raw.1'
fns <- list.files(data_path, pattern = '^[1-6]M|^[1-6]D', full.names = T)
bdat <- ReadAffy(filenames = fns, verbose=T)
bdat.d <- gcrma(bdat)

# add the pheno info
pheno <- data.frame(
  row.names = sampleNames(bdat.d),
  Infect = factor(c(rep('none',3), rep('CGSP14',3)), levels = c('none','CGSP14')),
  stringsAsFactors = F
)

pData(bdat.d) <- pheno

saveRDS(bdat.d, 'processedData/bewley_data_healthy_subjects_gcrma.RDS')

###### QC using the raw data################################

plm.fit <- fitPLM(bdat)

# compute the relative log expression (RLE) of each gene and plot the median values per sample as a barplot
median.rle <- RLE(plm.fit,main="RLE",type="stats")[1,]


median.rle.df <- data.frame(
  Sample = names(median.rle),
  med_RLE = median.rle
) 


median.rle.df %>%
  ggbarplot(x = 'Sample', y = 'med_RLE', title = 'Median RLE per sample',
            fill = '#00468B', color = '#00468B', ylab = 'Median RLE ') +
  ggsave('figs/median_rle.jpg', width = 6, height = 6, dpi = 300)


# compute the normalized unscaled standard error (NUSE) of each gene and plot the median values per sample as a histogram
median.nuse <- NUSE(plm.fit,main="NUSE",type="stats")[1,]

median.nuse.df <- data.frame(
  Sample = names(median.nuse),
  med_NUSE = median.nuse
) 


median.nuse.df %>%
  ggbarplot(x = 'Sample', y = 'med_NUSE', title = 'Median NUSE per sample',fill = '#00468B', color = '#00468B' ) +
  ggsave('figs/median_NUSE.jpg', width = 6, height = 6, dpi = 300)
  

# It appears that the every sample in the data is of high quality.

