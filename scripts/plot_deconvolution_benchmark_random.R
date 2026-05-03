library(tidyverse)

# Load Data
df <- readRDS("output/rds/cell_freq_comparison_results_random.rds")

# Ensure output directory exists
dir.create("output/images/random_genes", showWarnings = FALSE, recursive = TRUE)

# 1. Global Comparison
df_long <- df %>%
  pivot_longer(cols = c(Subject.ID, Age.months, Summed_CellFreq, Residuals),
               names_to = "Component", values_to = "Variance")

p1 <- ggplot(df_long, aes(x = Component, y = Variance, fill = Model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variance Component Comparison: CBC vs. Deconvolution (Random Genes)",
       subtitle = "Deconvolution captures more biological signal across the broader transcriptome",
       y = "Variance Explained (%)", x = "") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/random_genes/deconvolution_vs_cbc_all.png", p1, width = 10, height = 7)

# 2. Focused comparison on CellFreq
p2 <- ggplot(df, aes(x = Model, y = Summed_CellFreq, fill = Model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variance Explained by Cell Frequency Resolution (Random Genes)",
       subtitle = "High-resolution deconvolution significantly increases captured compositional signal",
       y = "Variance Explained (%)", x = "") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/random_genes/deconvolution_vs_cbc_cellfreq.png", p2, width = 8, height = 6)

# 3. Residual Reduction
p3 <- ggplot(df, aes(x = Model, y = Residuals, fill = Model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Reduction in Residual Noise (Random Genes)",
       subtitle = "Better cell-type resolution consistently reduces unexplained variation",
       y = "Residual Variance (%)", x = "") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/random_genes/deconvolution_vs_cbc_residuals.png", p3, width = 8, height = 6)

cat("Random gene deconvolution plots generated in output/images/random_genes/\n")
