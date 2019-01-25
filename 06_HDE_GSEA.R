# GSEA statistics for HDE using the MsigDB gene set database in tabular form
# use fgsea for all of the DEs.

library(tidyverse)
library(ggpubr)
library(kableExtra)
library(fgsea)

# using MsigDB 6.2
####################################################################################

# load the results that have human data for statistics ranking

final <- read_csv('output/fold_change_HDE_NDE_EDE.csv') %>%
  arrange(desc(human_t))

# prepare the ranks statistics (using the t and add names that are the entrez ID)

t_ranks <- final  %>%
  pull(human_t)

names(t_ranks) <- final$Human_Entrez_ID

####################################################################################

# prepare the gmt gsets
gmtPath <- '../../../MsigDB_geneset/entrez_6.2'
  
gsets <- list(
  h.all = gmtPathways(file.path(gmtPath, 'h.all.v6.2.entrez.gmt')),
  c1.all = gmtPathways(file.path(gmtPath, 'c1.all.v6.2.entrez.gmt')),
  c2.all = gmtPathways(file.path(gmtPath, 'c2.all.v6.2.entrez.gmt')),
  c2.cp.biocarta = gmtPathways(file.path(gmtPath, 'c2.cp.biocarta.v6.2.entrez.gmt')),
  c2.cp.kegg = gmtPathways(file.path(gmtPath, 'c2.cp.kegg.v6.2.entrez.gmt')),
  c2.cp.reactome = gmtPathways(file.path(gmtPath, 'c2.cp.reactome.v6.2.entrez.gmt')),
  c3.all = gmtPathways(file.path(gmtPath, 'c3.all.v6.2.entrez.gmt')),
  c4.all = gmtPathways(file.path(gmtPath, 'c4.all.v6.2.entrez.gmt')),
  c5.all = gmtPathways(file.path(gmtPath, 'c5.all.v6.2.entrez.gmt')),
  c6.all = gmtPathways(file.path(gmtPath, 'c6.all.v6.2.entrez.gmt')),
  c7.all = gmtPathways(file.path(gmtPath, 'c7.all.v6.2.entrez.gmt'))
)

####################################################################################
# running fgsea

res <- lapply(gsets, function(gs){
  fgsea(pathways = gs, 
        stats = t_ranks,
        minSize=15,
        maxSize=500,
        nperm=10000)
}) %>%
  bind_rows %>%
  arrange(desc(NES)) %>%
  mutate(leadingEdge = as.character(leadingEdge))

# output the result table
write_csv(res, 'output/HDE_GSEA_results.csv')
  
  