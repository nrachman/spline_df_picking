library(variancePartition)
library(edgeR)

# Load data
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

set.seed(42)
rand_genes <- sample(rownames(dge), 50)
meta_unscaled <- cbind(meta, cellfreqs)
names(meta_unscaled)[(ncol(meta)+1):ncol(meta_unscaled)] <- paste0("cellfreqs", names(cellfreqs))
cell_vars_unscaled <- paste0("cellfreqs", names(cellfreqs))

form_normal <- paste("~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent +", paste(cell_vars_unscaled, collapse=" + "))

form_dup <- paste("~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + mk_dup.PERCENT_DUPLICATION +", paste(cell_vars_unscaled, collapse=" + "))

v <- voom(dge[rand_genes, ], model.matrix(~Lib_prep_batches + Age.months, data = meta_unscaled))
res_normal <- as.data.frame(fitExtractVarPartModel(v, as.formula(form_normal), meta_unscaled))
med_normal <- median(rowSums(res_normal[, grep("cellfreq", colnames(res_normal)), drop=FALSE]))

res_dup <- as.data.frame(fitExtractVarPartModel(v, as.formula(form_dup), meta_unscaled))
med_dup <- median(rowSums(res_dup[, grep("cellfreq", colnames(res_dup)), drop=FALSE]))

cat("Normal Median: ", med_normal, "\n")
cat("Dup Median:    ", med_dup, "\n")
