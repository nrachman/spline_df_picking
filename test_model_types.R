library(variancePartition)
library(edgeR)
library(dplyr)
library(lme4)

cat("Loading data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]
colnames(dge) <- gsub("\\.", "_", colnames(dge))
rownames(dge$samples) <- gsub("\\.", "_", rownames(dge$samples))

cellfreqs <- readRDS("../../data/analysis_out/fardeep/fardeep_lm22_batch_corrected_cpm_nolog_add_cbc_mat_summed_groups.rds")
rownames(cellfreqs) <- gsub("\\.", "_", rownames(cellfreqs))
overlapping_subj <- intersect(colnames(dge), rownames(cellfreqs))

dge <- dge[, overlapping_subj]
cellfreqs <- cellfreqs[overlapping_subj, ]

meta <- dge$samples
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)
meta <- cbind(meta, as.data.frame(cellfreqs))

cont_vars <- c("Age.months", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", colnames(cellfreqs))
for(v in cont_vars) meta[[v]] <- scale(as.numeric(meta[[v]]))[,1]

set.seed(42)
rand_genes <- sample(rownames(dge), 250)
v <- voom(dge[rand_genes,], model.matrix(~Age.months, data=meta))

base_form <- "(1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent"

# 1. Paper Model (18 cell types: CBC + Deconv)
all_18_vars <- paste0("`", colnames(cellfreqs), "`", collapse=" + ")
form_18 <- as.formula(paste("~", base_form, "+", all_18_vars))

# 2. Deconv Only Model (11 cell types)
deconv_only <- setdiff(colnames(cellfreqs), c("Lymph.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt", "RBC", "PLT"))
deconv_11_vars <- paste0("`", deconv_only, "`", collapse=" + ")
form_11 <- as.formula(paste("~", base_form, "+", deconv_11_vars))

# 3. CBC Only Model (5 cell types)
cbc_only <- c("Lymph.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt")
cbc_5_vars <- paste0("`", cbc_only, "`", collapse=" + ")
form_5 <- as.formula(paste("~", base_form, "+", cbc_5_vars))

run_vp <- function(form, cell_vars) {
  res <- fitExtractVarPartModel(v, form, meta, control=lmerControl(check.conv.singular="ignore"))
  cell_cols <- intersect(colnames(res), cell_vars)
  med <- median(rowSums(res[, cell_cols, drop=FALSE]))
  return(med)
}

cat("Median Cell Freq Variance Explained:\n")
cat("18-Type (Original Paper Model):", run_vp(form_18, colnames(cellfreqs)), "\n")
cat("11-Type (Deconvolution Benchmark):", run_vp(form_11, deconv_only), "\n")
cat(" 5-Type (CBC Benchmark):", run_vp(form_5, cbc_only), "\n")
