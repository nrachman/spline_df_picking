library(tidyverse)

# Load Data
df <- readRDS("output/rds/vp_discrepancy_data.rds")

# Prepare for plotting
df_long <- df %>%
  pivot_longer(cols = c(Subject.ID, Age.months, Summed_CellFreq, Residuals),
               names_to = "Component", values_to = "Variance") %>%
  mutate(Component = factor(Component, levels = c("Subject.ID", "Summed_CellFreq", "Age.months", "Residuals")))

# 1. Comparison by Gene Set
p1 <- ggplot(df_long, aes(x = Component, y = Variance, fill = GeneSet)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variance Partitioning Discrepancy: HVGs vs. Random Genes",
       subtitle = "Highly Variable Genes are dominated by Subject Identity (individuality)",
       y = "Variance Explained (%)",
       x = "") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/discrepancy_hvg_vs_rand.png", p1, width = 10, height = 7)

# 2. Focused comparison on CellFreq
p2 <- ggplot(df, aes(x = GeneSet, y = Summed_CellFreq, fill = GeneSet)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Cell Frequency Signal across Gene Sets",
       subtitle = "The compositional signal is much stronger in the broader transcriptome",
       y = "Variance Explained by Cell Frequencies (%)",
       x = "") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/discrepancy_cellfreq_focus.png", p2, width = 8, height = 6)

cat("Discrepancy plots generated in output/images/\n")
