library(variancePartition)
library(edgeR)
library(tidyverse)

# 1. Load Data
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]
colnames(dge) <- gsub("\\.", "_", colnames(dge))
rownames(dge$samples) <- gsub("\\.", "_", rownames(dge$samples))

# Use the ORIGINAL summed groups file
cellfreqs <- readRDS("../../data/analysis_out/fardeep/fardeep_lm22_batch_corrected_cpm_nolog_add_cbc_mat_summed_groups.rds")
rownames(cellfreqs) <- gsub("\\.", "_", rownames(cellfreqs))
cell_vars <- colnames(cellfreqs)

overlapping_subj <- intersect(colnames(dge), rownames(cellfreqs))
dge_subset <- dge[, overlapping_subj]
meta <- dge_subset$samples
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

cf_df <- as.data.frame(cellfreqs[overlapping_subj, ])
meta <- cbind(meta, cf_df)

cont_vars <- c("Age.months", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", cell_vars)
for (v in cont_vars) meta[[v]] <- scale(as.numeric(meta[[v]]))[,1]

meta$cellfreqs <- as.matrix(meta[, cell_vars])

# 2. Compare HVG vs Random
cpm_full <- cpm(dge_subset, log=TRUE)
rv <- matrixStats::rowVars(cpm_full)
hvg <- names(sort(rv, decreasing=TRUE)[1:100])
set.seed(42)
rand_genes <- sample(rownames(dge_subset), 100)

form_with_subj <- ~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + cellfreqs

run_vp_test <- function(genes, label) {
  v <- voom(dge_subset[genes,], model.matrix(~Age.months, data=meta))
  res <- fitExtractVarPartModel(v, form_with_subj, meta)
  
  cell_cols <- grep("^cellfreqs", colnames(res))
  cat("\n---", label, "---\n")
  cat("Mean Age Variance: ", mean(res$Age.months), "\n")
  cat("Mean Subject Variance: ", mean(res$Subject.ID), "\n")
  cat("Mean Summed CellFreq Variance: ", mean(rowSums(res[, cell_cols])), "\n")
  cat("Mean Residual Variance: ", mean(res$Residuals), "\n")
}

run_vp_test(hvg, "Top 100 HVGs (With Subject)")
run_vp_test(rand_genes, "Random 100 Genes (With Subject)")

# 3. Effect of removing Subject.ID on Random Genes
form_no_subj <- ~ sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + cellfreqs
v_rand <- voom(dge_subset[rand_genes,], model.matrix(~Age.months, data=meta))
res_no_subj <- fitExtractVarPartModel(v_rand, form_no_subj, meta)
cell_cols <- grep("^cellfreqs", colnames(res_no_subj))
cat("\n--- Random 100 Genes (WITHOUT Subject) ---\n")
cat("Mean Age Variance: ", mean(res_no_subj$Age.months), "\n")
cat("Mean Summed CellFreq Variance: ", mean(rowSums(res_no_subj[, cell_cols])), "\n")
cat("Mean Residual Variance: ", mean(res_no_subj$Residuals), "\n")
