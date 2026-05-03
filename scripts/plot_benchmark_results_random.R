library(tidyverse)

# Load data from the correct path
results <- readRDS("output/rds/spline_benchmark_results_random.rds")

# Convert results list to tidy data frames
all_metrics <- list()
all_vp <- list()
all_pvals <- list()

for (name in names(results)) {
  res <- results[[name]]
  fam <- res$family
  comp <- res$complexity
  
  # Metrics (AIC/BIC)
  metrics <- as.data.frame(res$metrics) %>%
    mutate(gene = rownames(res$metrics), family = fam, complexity = comp)
  all_metrics[[name]] <- metrics
  
  # Variance Partition (Subject.ID)
  vp <- as.data.frame(res$vp) %>%
    mutate(gene = rownames(res$vp), family = fam, complexity = comp) %>%
    select(gene, family, complexity, Subject.ID)
  all_vp[[name]] <- vp
  
  # P-values (Age)
  pvals <- data.frame(gene = names(res$age_pvals), pval = res$age_pvals) %>%
    mutate(family = fam, complexity = comp)
  all_pvals[[name]] <- pvals
}

df_metrics <- bind_rows(all_metrics)
df_vp <- bind_rows(all_vp)
df_pvals <- bind_rows(all_pvals)

# Ensure output directory exists
dir.create("output/images/random_genes", showWarnings = FALSE, recursive = TRUE)

# 1. AIC Tournament Bar Chart ---------------------------------------------
df_min_aic <- df_metrics %>%
  group_by(gene, family) %>%
  filter(AIC == min(AIC, na.rm = TRUE)) %>%
  ungroup() %>%
  count(family, complexity) %>%
  group_by(family) %>%
  mutate(percent = n / sum(n) * 100)

p1 <- ggplot(df_min_aic, aes(x = factor(complexity), y = percent, fill = family)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "AIC Tournament: Winning Complexity (Random Genes)",
       x = "Complexity (df/n)", y = "Percent of Genes (%)",
       fill = "Spline Family")
ggsave("output/images/random_genes/aic_tournament.png", p1, width = 8, height = 5)

# 2. Subject Variance Stability -------------------------------------------
df_vp_summary <- df_vp %>%
  group_by(family, complexity) %>%
  summarise(median_subject = median(Subject.ID, na.rm = TRUE))

p2 <- ggplot(df_vp_summary, aes(x = complexity, y = median_subject, color = family)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Subject Variance Stability (Random Genes)",
       x = "Complexity (df/n)", y = "Median Subject Variance (%)") +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent)
ggsave("output/images/random_genes/subject_variance_stability.png", p2, width = 8, height = 5)

# 3. Age Inference Power --------------------------------------------------
df_sig <- df_pvals %>%
  group_by(family, complexity) %>%
  summarise(n_sig = sum(pval < 0.05, na.rm = TRUE))

p3 <- ggplot(df_sig, aes(x = complexity, y = n_sig, color = family)) +
  geom_line() + geom_point() +
  theme_minimal() +
  labs(title = "Age Inference Power (Random Genes)",
       x = "Complexity (df/n)", y = "Number of Significant Genes (p < 0.05)")
ggsave("output/images/random_genes/age_inference_power.png", p3, width = 8, height = 5)

# 4. Subject Variance Distribution ----------------------------------------
p4 <- ggplot(df_vp, aes(x = factor(complexity), y = Subject.ID, fill = family)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "Subject Variance Distribution (Random Genes)",
       x = "Complexity (df/n)", y = "Subject Variance (%)") +
  scale_y_continuous(labels = scales::percent)
ggsave("output/images/random_genes/subject_variance_distribution.png", p4, width = 10, height = 6)

# 5. Subject Variance Correlation (between complexities) ------------------
# Focus on ns for simplicity
df_vp_wide <- df_vp %>%
  filter(family == "ns") %>%
  select(gene, complexity, Subject.ID) %>%
  pivot_wider(names_from = complexity, values_from = Subject.ID) %>%
  select(-gene)

cor_mat <- cor(df_vp_wide, use = "complete.obs")
p5 <- reshape2::melt(cor_mat) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.5) +
  theme_minimal() +
  labs(title = "Subject Variance Correlation Heatmap (Random Genes)",
       x = "Complexity", y = "Complexity", fill = "Pearson R")
ggsave("output/images/random_genes/subject_variance_correlation.png", p5, width = 10, height = 5)

# 6. AIC Comparison (ns vs lspline) ---------------------------------------
df_aic_diff <- df_metrics %>%
  select(gene, family, complexity, AIC) %>%
  pivot_wider(names_from = family, values_from = AIC) %>%
  mutate(AIC_Diff = ns - lspline) # Negative means ns is better

p6 <- ggplot(df_aic_diff, aes(x = factor(complexity), y = AIC_Diff)) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "AIC Comparison: ns vs lspline (Random Genes)",
       subtitle = "Negative values indicate Natural Splines (ns) provide better fit",
       x = "Complexity (df/n)", y = "AIC(ns) - AIC(lspline)")
ggsave("output/images/random_genes/aic_comparison_ns_vs_lspline.png", p6, width = 8, height = 5)

# 7. Global Variance Partitioning (at df=4) -------------------------------
# This requires running a full VP for one configuration, we can approximate 
# from the results if we had all components. For now skip or use dummy if not available.
# The benchmark results only saved Subject.ID and Age p-vals.
# We will use the discrepancy audit data for the global view if needed.

cat("Random gene benchmark plots generated in output/images/random_genes/\n")
