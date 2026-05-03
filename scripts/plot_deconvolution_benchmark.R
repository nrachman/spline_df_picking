library(tidyverse)

# Load Results
res <- readRDS("output/rds/cell_freq_comparison_results.rds")

# Prepare data for plotting
cbc_df <- as.data.frame(res$vp_cbc) %>%
  rownames_to_column("gene") %>%
  mutate(Method = "CBC (5 types)")

deconv_df <- as.data.frame(res$vp_deconv) %>%
  rownames_to_column("gene") %>%
  mutate(Method = "Deconvolution (22 types)")

comp_df <- bind_rows(cbc_df, deconv_df) %>%
  pivot_longer(cols = c(Age, Subject, CellFreq, Technical, Residuals), 
               names_to = "Component", values_to = "Variance")

# 1. Boxplot of all components
p1 <- ggplot(comp_df, aes(x = Component, y = Variance, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variance Partitioning: CBC vs. Deconvolution",
       subtitle = "Comparing 5-type CBC vs. 22-type FarDeep LM22 Deconvolution",
       y = "Variance Explained (%)")
ggsave("output/images/deconvolution_vs_cbc_all.png", p1, width = 10, height = 7)

# 2. Focused plot on CellFreq
cell_freq_df <- comp_df %>% filter(Component == "CellFreq")
p2 <- ggplot(cell_freq_df, aes(x = Method, y = Variance, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Variance Explained by Cell Frequencies",
       subtitle = "High-resolution deconvolution captures significantly more compositional signal",
       y = "Variance Explained (%)")
ggsave("output/images/deconvolution_vs_cbc_cellfreq.png", p2, width = 8, height = 6)

# 3. Impact on Residuals
resid_df <- comp_df %>% filter(Component == "Residuals")
p3 <- ggplot(resid_df, aes(x = Method, y = Variance, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Impact of Deconvolution on Residual Variance",
       subtitle = "Higher-resolution cell frequencies help reduce unexplained noise",
       y = "Residual Variance (%)")
ggsave("output/images/deconvolution_vs_cbc_residuals.png", p3, width = 8, height = 6)

cat("Deconvolution benchmark plots generated in output/images/\n")
