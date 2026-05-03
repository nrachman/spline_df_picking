library(variancePartition)
library(edgeR)
library(tidyverse)
library(lme4)
library(matrixStats)

# 1. Load and Sync Data ----------------------------------------------------
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]
colnames(dge) <- gsub("\\.", "_", colnames(dge))
rownames(dge$samples) <- gsub("\\.", "_", rownames(dge$samples))

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
meta$cellfreqs_mat <- as.matrix(meta[, cell_vars])

# 2. Gene Selection --------------------------------------------------------
cpm_full <- cpm(dge_subset, log=TRUE)
rv <- rowVars(cpm_full)
hvg <- names(sort(rv, decreasing=TRUE)[1:250])
set.seed(42)
rand_genes <- sample(setdiff(rownames(dge_subset), hvg), 250)

# 3. Model Definitions -----------------------------------------------------
form_full <- ~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + cellfreqs_mat

# Helper to run VP and extract summed cellfreq
run_vp_batch <- function(genes, label) {
  v <- voom(dge_subset[genes,], model.matrix(~Age.months, data=meta))
  res <- fitExtractVarPartModel(v, form_full, meta)
  
  cell_cols <- grep("^cellfreqs", colnames(res))
  df <- as.data.frame(res) %>%
    mutate(gene = rownames(res),
           GeneSet = label,
           Summed_CellFreq = rowSums(res[, cell_cols])) %>%
    select(gene, GeneSet, Subject.ID, Age.months, Summed_CellFreq, Residuals)
  return(df)
}

cat("Running VP on HVGs...\n")
hvg_res <- run_vp_batch(hvg, "Top 250 HVGs")

cat("Running VP on Random Genes...\n")
rand_res <- run_vp_batch(rand_genes, "Random 250 Genes")

final_res <- bind_rows(hvg_res, rand_res)

dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)
saveRDS(final_res, "output/rds/vp_discrepancy_data.rds")
cat("Data saved to output/rds/vp_discrepancy_data.rds\n")
