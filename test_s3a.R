vp <- readRDS("../../data/analysis_out/variancePartition/varPart_include_cbc_fardeep_lm22_summed.rds")
vp <- data.frame(vp)
cell_cols <- startsWith(colnames(vp), "cellfreq")
vp_sub <- vp[, !cell_cols]
vp_sub$CellFreqs <- rowSums(vp[, cell_cols])
cat("Medians in paper:\n")
print(apply(vp_sub, 2, median))
