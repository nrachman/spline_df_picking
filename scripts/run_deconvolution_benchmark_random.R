library(variancePartition)
library(edgeR)
library(tidyverse)
library(lme4)
library(matrixStats)
library(BiocParallel)

# 1. Load and Sync Data ----------------------------------------------------
options(warn = -1)
cat("Loading data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

# Load CBC data for Mid.pcnt
cat("Loading CBC data...\n")
cbc_data <- read.table("../../data/inputs/cbc/cbc_cleaned.tsv", header=TRUE, sep="\t")

# Load deconvolution data
cat("Loading Deconvolution data...\n")
cellfreqs <- readRDS("../../data/analysis_out/fardeep/fardeep_lm22_batch_corrected_cpm_nolog_add_cbc_mat_summed_groups.rds")

# Strip _S suffixes for matching
strip_s <- function(x) gsub("_S[0-9]+$", "", x)
colnames(dge) <- strip_s(colnames(dge))
rownames(cellfreqs) <- strip_s(rownames(cellfreqs))
# CBC IDs use dots (6895.10PX) while others use underscores (6895_10PX)
cbc_data$Sample_ID <- gsub("\\.", "_", strip_s(cbc_data$Sample_ID))

# Sync all three
overlapping_samples <- intersect(intersect(colnames(dge), rownames(cellfreqs)), cbc_data$Sample_ID)
cat("Overlapping samples found:", length(overlapping_samples), "\n")
if(length(overlapping_samples) == 0) stop("Zero overlapping samples! Check ID formats.")
dge_subset <- dge[, overlapping_samples]
cellfreqs <- cellfreqs[overlapping_samples, ]
cbc_subset <- cbc_data[match(overlapping_samples, cbc_data$Sample_ID), ]

meta <- dge_subset$samples
meta$Subject.ID <- factor(meta$Subject.ID)
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

# Add cell frequencies and CBC variables to metadata
cf_df <- as.data.frame(cellfreqs)
cbc_vars_to_add <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
meta <- cbind(meta, cf_df)
for(v in cbc_vars_to_add) {
  meta[[v]] <- cbc_subset[[v]]
}

# Define variable sets
cbc_vars <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt") # Use the core 3 or 5? Original Analysis 3 used Lymph, Mid, Gran
deconv_vars <- setdiff(colnames(cellfreqs), c("Lymph.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt", "RBC", "PLT"))

# Scale all continuous covariates
cont_vars <- c("Age.months", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
               cbc_vars_to_add, colnames(cf_df))
cont_vars <- unique(intersect(cont_vars, colnames(meta)))

for (v in cont_vars) {
  meta[[v]] <- scale(as.numeric(meta[[v]]))[,1]
}

# 2. Gene Selection (1000 Random Genes) ------------------------------------
cat("Selecting 1000 Random Genes...\n")
cpm_full <- cpm(dge_subset, log=TRUE)
rv <- rowVars(cpm_full)
hvg_excluded <- names(sort(rv, decreasing=TRUE)[1:2000])
set.seed(42)
rand_genes <- sample(setdiff(rownames(dge_subset), hvg_excluded), 1000)
v_rand <- voom(dge_subset[rand_genes,], model.matrix(~Age.months, data=meta))

# 3. Benchmark Formulas ----------------------------------------------------
base_form <- "(1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent"

# Use backticks for cell variables
cbc_mat_str <- paste0("`", cbc_vars, "`", collapse = " + ")
deconv_mat_str <- paste0("`", deconv_vars, "`", collapse = " + ")

form_cbc <- as.formula(paste("~", base_form, "+", cbc_mat_str))
form_deconv <- as.formula(paste("~", base_form, "+", deconv_mat_str))

ctrl <- lmerControl(optimizer = "bobyqa", 
                    optCtrl = list(maxfun = 100000),
                    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))

# 4. Run VP Benchmark ------------------------------------------------------
cat("Running CBC Model Variance Partitioning...\n")
res_cbc <- suppressMessages(suppressWarnings(
  fitExtractVarPartModel(v_rand, form_cbc, meta, BPPARAM=SerialParam(), control=ctrl)
))

cat("Running Deconvolution Model Variance Partitioning...\n")
res_deconv <- suppressMessages(suppressWarnings(
  fitExtractVarPartModel(v_rand, form_deconv, meta, BPPARAM=SerialParam(), control=ctrl)
))

# 5. Aggregate and Save ----------------------------------------------------
summarize_vp <- function(res, label, cell_vars) {
  df <- as.data.frame(res)
  cell_cols <- intersect(colnames(df), cell_vars)
  df$Summed_CellFreq <- rowSums(df[, cell_cols, drop=FALSE])
  df$Model <- label
  df$gene <- rownames(df)
  return(df %>% select(gene, Model, Subject.ID, Age.months, Summed_CellFreq, Residuals))
}

results_cbc <- summarize_vp(res_cbc, "CBC (3-type)", cbc_vars)
results_deconv <- summarize_vp(res_deconv, "Deconvolution (11-type)", deconv_vars)

final_results <- bind_rows(results_cbc, results_deconv)
saveRDS(final_results, "output/rds/cell_freq_comparison_results_random.rds")

cat("Benchmark complete. Results saved to output/rds/cell_freq_comparison_results_random.rds\n")
