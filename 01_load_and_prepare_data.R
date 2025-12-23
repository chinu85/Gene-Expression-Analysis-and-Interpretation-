# 01_load_and_prepare_data.R

library(tidyverse)
library(data.table)

# setting up paths
# assuming the data is in this subfolder
data_dir <- "brca_tcga_pan_can_atlas_2018"

# create the output folder if it's not there
dir.create("results", showWarnings = FALSE)

# 1. load the rna seq data
# this file is pretty big so fread is faster
rna <- fread(file.path(data_dir, "data_mrna_seq_v2_rsem.txt"), data.table = FALSE)

# cleanup the duplicate genes if possible
# usually the first column is gene symbol, check for dupes
# remove rows with empty or NA gene symbol
rna <- rna[!is.na(rna$Hugo_Symbol) & rna$Hugo_Symbol != "", ]

if (any(duplicated(rna$Hugo_Symbol))) {
  rna <- rna[!duplicated(rna$Hugo_Symbol), ]
}

# set row names to genes and drop the first two text columns
rownames(rna) <- rna$Hugo_Symbol
expr_mat <- as.matrix(rna[, -c(1, 2)])

# 2. clinical data
# skip the first 4 lines cause they are header descriptions
clinical <- fread(file.path(data_dir, "data_clinical_patient.txt"), data.table = FALSE, skip = 4)

# 3. cna data for her2 status
cna <- fread(file.path(data_dir, "data_cna.txt"), data.table = FALSE)

# check for and remove duplicates/NAs in CNA data
cna <- cna[!is.na(cna$Hugo_Symbol) & cna$Hugo_Symbol != "", ]
if (any(duplicated(cna$Hugo_Symbol))) {
  cna <- cna[!duplicated(cna$Hugo_Symbol), ]
}

rownames(cna) <- cna$Hugo_Symbol
cna_mat <- as.matrix(cna[, -c(1, 2)])

# we need ERBB2 to determine status
erbb2_vals <- cna_mat["ERBB2", ]

# define status based on amplification positive > 0
her2_status <- ifelse(erbb2_vals > 0, "Amplified", "NotAmplified")

# make a dataframe we can merge later
her2_df <- data.frame(
  patient_id = substr(names(erbb2_vals), 1, 12),
  HER2_status = her2_status,
  stringsAsFactors = FALSE
)

# 4. merging it all
# only keep the primary tumor samples (01) from the expression data
samples <- colnames(expr_mat)
keep <- substr(samples, 14, 15) == "01"
expr_mat <- expr_mat[, keep]

# match patient IDs
patients <- substr(colnames(expr_mat), 1, 12)

meta <- data.frame(
  sample_id = colnames(expr_mat),
  patient_id = patients,
  stringsAsFactors = FALSE
)

# merge her2 and clinical info
meta2 <- meta %>%
  left_join(her2_df, by = "patient_id") %>%
  left_join(clinical, by = c("patient_id" = "PATIENT_ID"))

# dropping samples that don't have her2 status
meta2 <- meta2 %>% filter(!is.na(HER2_status))

# make sure expression matrix matches meta
common <- intersect(colnames(expr_mat), meta2$sample_id)
expr_mat <- expr_mat[, common]
meta2 <- meta2[match(common, meta2$sample_id), ]

# round to integers for deseq2 later
count_mat <- round(expr_mat)

saveRDS(count_mat, "results/expr_mat_round.rds")
saveRDS(meta2, "results/meta2.rds")
