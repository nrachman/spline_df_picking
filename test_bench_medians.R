library(dplyr)
df <- readRDS("output/rds/cell_freq_comparison_results_random.rds")
df %>% group_by(Model) %>% summarise(
  median_CellFreq = median(Summed_CellFreq),
  median_Subject = median(Subject.ID),
  median_Residuals = median(Residuals)
)
