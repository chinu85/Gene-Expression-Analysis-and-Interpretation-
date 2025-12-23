# 05_lasso_cox_glmnet.R

library(tidyverse)
library(survival)
library(glmnet)

# loading inputs
vst_mat <- readRDS("results/vst_mat.rds")
dds <- readRDS("results/dds.rds")
deg <- readRDS("results/DEG_significant.rds")

# prep metadata from the deseq object
meta <- as.data.frame(colData(dds))

# making sure survival columns are numeric
meta$OS_MONTHS <- as.numeric(meta$OS_MONTHS)

# handling the event column
# sometimes it's text like "Deceased", so converting to 1/0
meta$OS_EVENT <- ifelse(grepl("DECEASED|Dead|1", meta$OS_STATUS, ignore.case = TRUE), 1, 0)

# remove samples with missing survival info
valid <- !is.na(meta$OS_MONTHS) & meta$OS_MONTHS > 0 & !is.na(meta$OS_EVENT)
meta_surv <- meta[valid, ]

# match up the vst matrix
common <- intersect(colnames(vst_mat), rownames(meta_surv))
vst_surv <- vst_mat[, common]
meta_surv <- meta_surv[common, ]

# feature selection
# just taking the top 50 significant genes to use as predictors
features <- deg %>%
  dplyr::arrange(padj) %>%
  head(100) %>%
  pull(gene)

# intersect with available genes in vst
features <- intersect(features, rownames(vst_surv))

# set up the matrix for glmnet
X <- t(vst_surv[features, ])
y <- Surv(meta_surv$OS_MONTHS, meta_surv$OS_EVENT)

# run lasso cox
set.seed(1998)
cv_fit <- cv.glmnet(X, y, family = "cox", alpha = 1)

# plotting the cv curve locally
pdf("results/Lasso_CV_curve.pdf")
plot(cv_fit)
dev.off()

# get the coefficients for the best lambda
coefs <- coef(cv_fit, s = "lambda.min")
active_idx <- which(coefs != 0)
active_genes <- rownames(coefs)[active_idx]
active_coefs <- coefs[active_idx]

results_df <- data.frame(gene = active_genes, coef = active_coefs)

write.csv(results_df, "results/glmnet_selected_genes.csv", row.names = FALSE)
