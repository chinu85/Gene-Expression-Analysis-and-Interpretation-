# 04_vst_PCA_heatmap.R

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggplot2)

dds <- readRDS("results/dds.rds")
de_res <- readRDS("results/DEG_significant.rds")

# vst transform for plots
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
saveRDS(vst_mat, "results/vst_mat.rds")

# 1. PCA Plot
# use top 500 variable genes
rv <- apply(vst_mat, 1, var)
select <- order(rv, decreasing = TRUE)[1:500]

pca <- prcomp(t(vst_mat[select, ]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)

# prep data for ggplot
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  HER2 = colData(dds)$HER2_status
)

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = HER2)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "%")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "%")) +
  coord_fixed() +
  theme_bw()

pdf("results/PCA_plot.pdf", width = 7, height = 6)
print(p)
dev.off()

# 2. Heatmap
# taking the top 50 significant genes
top_genes <- de_res %>%
  dplyr::arrange(padj) %>%
  head(50) %>%
  pull(gene)

# filtering the vst matrix
heat_mat <- vst_mat[intersect(top_genes, rownames(vst_mat)), ]

# annotations for the heatmap
anno <- data.frame(HER2 = colData(dds)$HER2_status)
rownames(anno) <- colnames(vst_mat)

pdf("results/Heatmap_top50_DEGs.pdf", width = 8, height = 10)
pheatmap(heat_mat,
         annotation_col = anno,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE)
dev.off()
