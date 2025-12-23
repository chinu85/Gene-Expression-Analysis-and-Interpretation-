# 03_pathway_enrichment.R

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

# load the sig genes
sig_genes <- readRDS("results/DEG_significant.rds")

# we need entrez ids for clusterprofiler
genes <- sig_genes$gene

ids <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
valid_ids <- na.omit(ids)

# run GO enrichment
ego <- enrichGO(gene = valid_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "BH",
                readable = TRUE)

if (!is.null(ego)) {
  write.csv(as.data.frame(ego), "results/GO_enrichment_DEGs.csv", row.names = FALSE)
}

# run KEGG
kk <- enrichKEGG(gene = valid_ids, organism = 'hsa')

if (!is.null(kk)) {
  write.csv(as.data.frame(kk), "results/KEGG_enrichment_DEGs.csv", row.names = FALSE)
}
