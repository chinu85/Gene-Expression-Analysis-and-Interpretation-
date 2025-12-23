library(tidyverse)
library(data.table)
data_dir <- "brca_tcga_pan_can_atlas_2018"
dir.create("results", showWarnings = FALSE)
rna <- fread(file.path(data_dir, "data_mrna_seq_v2_rsem.txt"), data.table = FALSE)
rna <- rna[!is.na(rna$Hugo_Symbol) & rna$Hugo_Symbol != "", ]
if (any(duplicated(rna$Hugo_Symbol))) {
  rna <- rna[!duplicated(rna$Hugo_Symbol), ]
}
rownames(rna) <- rna$Hugo_Symbol
expr_mat <- as.matrix(rna[, -c(1, 2)])
clinical <- fread(file.path(data_dir, "data_clinical_patient.txt"), data.table = FALSE, skip = 4)
cna <- fread(file.path(data_dir, "data_cna.txt"), data.table = FALSE)
cna <- cna[!is.na(cna$Hugo_Symbol) & cna$Hugo_Symbol != "", ]
if (any(duplicated(cna$Hugo_Symbol))) {
  cna <- cna[!duplicated(cna$Hugo_Symbol), ]
}
rownames(cna) <- cna$Hugo_Symbol
cna_mat <- as.matrix(cna[, -c(1, 2)])
erbb2_vals <- cna_mat["ERBB2", ]
her2_status <- ifelse(erbb2_vals > 0, "Amplified", "NotAmplified")
her2_df <- data.frame(patient_id = substr(names(erbb2_vals), 1, 12), HER2_status = her2_status, stringsAsFactors = FALSE)
samples <- colnames(expr_mat)
keep <- substr(samples, 14, 15) == "01"
expr_mat <- expr_mat[, keep]
patients <- substr(colnames(expr_mat), 1, 12)
meta <- data.frame(sample_id = colnames(expr_mat), patient_id = patients, stringsAsFactors = FALSE)
meta2 <- meta %>%
  left_join(her2_df, by = "patient_id") %>%
  left_join(clinical, by = c("patient_id" = "PATIENT_ID"))
meta2 <- meta2 %>% filter(!is.na(HER2_status))
common <- intersect(colnames(expr_mat), meta2$sample_id)
expr_mat <- expr_mat[, common]
meta2 <- meta2[match(common, meta2$sample_id), ]
count_mat <- round(expr_mat)
saveRDS(count_mat, "results/expr_mat_round.rds")
saveRDS(meta2, "results/meta2.rds")
