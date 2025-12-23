library(tidyverse)
library(DESeq2)
count_dat <- readRDS("results/expr_mat_round.rds")
meta <- readRDS("results/meta2.rds")
meta$HER2_status <- factor(meta$HER2_status)
meta$HER2_status <- relevel(meta$HER2_status, ref = "NotAmplified")
dds <- DESeqDataSetFromMatrix(countData = count_dat, colData = meta, design = ~HER2_status)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep, ]
dds <- DESeq(dds)
saveRDS(dds, "results/dds.rds")
res <- results(dds, contrast = c("HER2_status", "Amplified", "NotAmplified"))
res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$gene <- rownames(res_df)
res_df <- res_df %>% dplyr::relocate(gene)
write.csv(res_df, "results/DE_results_full.csv", row.names = FALSE)
sig_genes <- res_df %>% dplyr::filter(padj < 0.05)
saveRDS(sig_genes, "results/DEG_significant.rds")
top10 <- sig_genes %>%
  dplyr::mutate(abs_fc = abs(log2FoldChange)) %>%
  dplyr::arrange(desc(abs_fc)) %>%
  dplyr::slice(1:10) %>%
  dplyr::select(-abs_fc)
write.csv(top10, "results/DE_top10_FC.csv", row.names = FALSE)
