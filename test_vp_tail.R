library(dplyr)

cat("1. Original Paper (16,599 genes):\n")
vp <- readRDS("../../data/analysis_out/variancePartition/varPart_include_cbc_fardeep_lm22_summed.rds")
vp <- data.frame(vp)
cell_cols <- startsWith(colnames(vp), "cellfreq")
cell_vp <- rowSums(vp[, cell_cols])
cat(">50% variance: ", mean(cell_vp > 0.5) * 100, "%\n")
cat(">25% variance: ", mean(cell_vp > 0.25) * 100, "%\n")

cat("\n2. Random Genes Benchmark (1000 genes):\n")
df_bench <- readRDS("output/rds/cell_freq_comparison_results_random.rds") %>%
  filter(Model == "Deconvolution (11-type)")
cat(">50% variance: ", mean(df_bench$Summed_CellFreq > 0.5) * 100, "%\n")
cat(">25% variance: ", mean(df_bench$Summed_CellFreq > 0.25) * 100, "%\n")

cat("\n3. Discrepancy Random Genes (250 genes from voom):\n")
df_discrep <- readRDS("output/rds/vp_discrepancy_data.rds") %>%
  filter(GeneSet == "Random 250 Genes")
cat(">50% variance: ", mean(df_discrep$Summed_CellFreq > 0.5) * 100, "%\n")
cat(">25% variance: ", mean(df_discrep$Summed_CellFreq > 0.25) * 100, "%\n")
