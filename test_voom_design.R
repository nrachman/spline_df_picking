library(variancePartition)
library(edgeR)
library(tidyverse)

# 1. Load Data
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

# 2. Genes
set.seed(42)
rand_genes <- sample(rownames(dge_subset), 250)

# 3. Model
form_full <- ~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + cellfreqs_mat

cat("Testing voom with ~Age.months\n")
v1 <- voom(dge_subset[rand_genes,], model.matrix(~Age.months, data=meta))
res1 <- fitExtractVarPartModel(v1, form_full, meta)
cell_cols1 <- grep("^cellfreqs", colnames(res1))
med1 <- median(rowSums(res1[, cell_cols1]))

cat("Testing voom with ~Lib_prep_batches + Age.months\n")
v2 <- voom(dge_subset[rand_genes,], model.matrix(~Lib_prep_batches + Age.months, data=meta))
res2 <- fitExtractVarPartModel(v2, form_full, meta)
cell_cols2 <- grep("^cellfreqs", colnames(res2))
med2 <- median(rowSums(res2[, cell_cols2]))

cat("Median CellFreq Variance:\n")
cat("~Age.months: ", med1, "\n")
cat("~Lib_prep_batches + Age.months: ", med2, "\n")
