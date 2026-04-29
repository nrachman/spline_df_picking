library(tidyverse)

# Load Results
results <- readRDS("output/rds/advanced_variance_benchmarks.rds")

# 1. Compare Standard vs. Realized Variance Partitioning -------------------
# Focus on ns_df4 for this comparison
if ("ns_df4" %in% names(results)) {
  res <- results[["ns_df4"]]
  
  std <- as.data.frame(res$vp_std) %>%
    rownames_to_column("gene") %>%
    select(gene, Subject.ID, Age = starts_with("ns(")) %>%
    mutate(Method = "Standard")
  
  realized <- as.data.frame(res$vp_realized) %>%
    rownames_to_column("gene") %>%
    select(gene, Subject.ID = Subject.ID, Age) %>%
    mutate(Method = "Realized")
  
  comp_df <- bind_rows(std, realized) %>%
    pivot_longer(cols = c(Subject.ID, Age), names_to = "Component", values_to = "Variance")
  
  p1 <- ggplot(comp_df, aes(x = Component, y = Variance, fill = Method)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Variance Partitioning Method Comparison (ns, df=4)",
         subtitle = "Realized Predictions provides a more conservative empirical scale",
         y = "Variance Explained (%)")
  ggsave("output/images/vp_method_comparison.png", p1, width = 8, height = 6)
}

# 2. Compare Subject vs. Binned Age Variance -------------------------------
if ("age_random_binned" %in% names(results)) {
  res <- results[["age_random_binned"]]
  vp <- as.data.frame(res$vp_realized) %>%
    rownames_to_column("gene") %>%
    select(gene, Subject.ID, Age = Age_Binned) %>%
    pivot_longer(-gene, names_to = "Component", values_to = "Variance")
  
  p2 <- ggplot(vp, aes(x = Component, y = Variance, fill = Component)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Subject vs. Binned Age Variance (Random Effects)",
         subtitle = "Age binned into 6-sample groups (~115 bins total)",
         y = "Variance Explained (%)")
  ggsave("output/images/subject_vs_binned_age_variance.png", p2, width = 8, height = 6)
}

# 3. Model Architecture Comparison (AIC) -----------------------------------
all_metrics <- list()
for (name in names(results)) {
  all_metrics[[name]] <- as.data.frame(results[[name]]$metrics) %>%
    mutate(gene = rownames(results[[name]]$metrics), Model = name)
}
df_metrics <- bind_rows(all_metrics)

# Filter for genes that have results in all models
common_genes <- df_metrics %>%
  group_by(gene) %>%
  summarise(n = n()) %>%
  filter(n == length(results)) %>%
  pull(gene)

df_metrics_filtered <- df_metrics %>% filter(gene %in% common_genes)

p3 <- ggplot(df_metrics_filtered, aes(x = Model, y = AIC, fill = Model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Model Complexity Comparison (AIC)",
       subtitle = "Lower AIC indicates better fit/parsimony balance",
       y = "AIC")
ggsave("output/images/model_aic_comparison.png", p3, width = 10, height = 7)

# 4. Variance explained by Age across different architectures -------------
all_age_vp <- list()
for (name in names(results)) {
  vp <- as.data.frame(results[[name]]$vp_realized) %>%
    rownames_to_column("gene") %>%
    select(gene, Age) %>%
    mutate(Model = name)
  all_age_vp[[name]] <- vp
}
df_age_vp <- bind_rows(all_age_vp) %>% filter(gene %in% common_genes)

p4 <- ggplot(df_age_vp, aes(x = Model, y = Age, fill = Model)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Variance Explained by Age across Model Architectures",
       subtitle = "Using Realized Predictions (Empirical Scale)",
       y = "Variance Explained by Age (%)")
ggsave("output/images/age_variance_by_model.png", p4, width = 10, height = 7)

cat("Advanced plots generated in output/images/\n")
