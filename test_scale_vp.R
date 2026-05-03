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

# Random Genes
set.seed(42)
rand_genes <- sample(rownames(dge), 50)

# Unscaled Meta
meta_unscaled <- cbind(meta, cellfreqs)
# Need to name them carefully
names(meta_unscaled)[(ncol(meta)+1):ncol(meta_unscaled)] <- paste0("cellfreqs", names(cellfreqs))
cell_vars_unscaled <- paste0("cellfreqs", names(cellfreqs))
form_str_unscaled <- paste("~ (1|Subject.ID) + sex.numeric + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + Age.months + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent +", paste(cell_vars_unscaled, collapse=" + "))

v_unscaled <- voom(dge[rand_genes, ], model.matrix(~Lib_prep_batches + Age.months, data = meta_unscaled))
res_unscaled <- as.data.frame(fitExtractVarPartModel(v_unscaled, as.formula(form_str_unscaled), meta_unscaled))
med_unscaled <- median(rowSums(res_unscaled[, grep("cellfreq", colnames(res_unscaled)), drop=FALSE]))


# Scaled Meta
meta_scaled <- meta_unscaled
cont_vars <- c("Age.months", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", cell_vars_unscaled)
for (v in cont_vars) {
  meta_scaled[[v]] <- scale(as.numeric(meta_scaled[[v]]))[,1]
}

v_scaled <- voom(dge[rand_genes, ], model.matrix(~Lib_prep_batches + Age.months, data = meta_scaled))
res_scaled <- as.data.frame(fitExtractVarPartModel(v_scaled, as.formula(form_str_unscaled), meta_scaled))
med_scaled <- median(rowSums(res_scaled[, grep("cellfreq", colnames(res_scaled)), drop=FALSE]))

cat("Unscaled Median: ", med_unscaled, "\n")
cat("Scaled Median:   ", med_scaled, "\n")
