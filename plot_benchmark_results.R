library(tidyverse)

results <- readRDS("spline_benchmark_results.rds")

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
  labs(title = "AIC Tournament: Winning Complexity",
       x = "Complexity (df/n)", y = "Percent of Genes (%)",
       fill = "Spline Family")
ggsave("aic_tournament.png", p1, width = 8, height = 5)

# 2. Subject Variance Robustness -----------------------------------------
df_vp_median <- df_vp %>%
  group_by(family, complexity) %>%
  summarise(median_subject_vp = median(Subject.ID, na.rm = TRUE), .groups = "drop")

p2 <- ggplot(df_vp_median, aes(x = complexity, y = median_subject_vp, color = family)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  labs(title = "Subject Variance Stability",
       x = "Complexity (df/n)", y = "Median Subject Variance (%)",
       color = "Spline Family")
ggsave("subject_variance_stability.png", p2, width = 8, height = 5)

# 3. Age Inference Robustness --------------------------------------------
df_sig_age <- df_pvals %>%
  group_by(family, complexity) %>%
  summarise(sig_count = sum(p.adjust(pval, method = "fdr") < 0.05, na.rm = TRUE), .groups = "drop")

p3 <- ggplot(df_sig_age, aes(x = complexity, y = sig_count, color = family)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:10) +
  theme_minimal() +
  labs(title = "Age Inference Power",
       x = "Complexity (df/n)", y = "Significant Age Genes (FDR < 0.05)",
       color = "Spline Family")
ggsave("age_inference_power.png", p3, width = 8, height = 5)

# 4. Distribution of Subject Variance Explained -----------------------------
p4 <- ggplot(df_vp, aes(x = factor(complexity), y = Subject.ID, fill = family)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Subject Variance Explained",
       x = "Complexity (df/n)", y = "Variance Explained by Subject.ID (%)",
       fill = "Spline Family")
ggsave("subject_variance_distribution.png", p4, width = 10, height = 6)

# 5. Correlation of Subject Variance Across Degrees of Freedom --------------
df_vp_wide_ns <- df_vp %>% filter(family == "ns") %>% select(gene, complexity, Subject.ID) %>%
  pivot_wider(names_from = complexity, values_from = Subject.ID, names_prefix = "df_")
cor_mat_ns <- cor(df_vp_wide_ns %>% select(-gene), use = "pairwise.complete.obs")
cor_df_ns <- as.data.frame(as.table(cor_mat_ns)) %>% mutate(family = "ns")

df_vp_wide_lspline <- df_vp %>% filter(family == "lspline") %>% select(gene, complexity, Subject.ID) %>%
  pivot_wider(names_from = complexity, values_from = Subject.ID, names_prefix = "n_")
cor_mat_lspline <- cor(df_vp_wide_lspline %>% select(-gene), use = "pairwise.complete.obs")
cor_df_lspline <- as.data.frame(as.table(cor_mat_lspline)) %>% mutate(family = "lspline")

cor_df <- bind_rows(cor_df_ns, cor_df_lspline)
cor_df$Var1 <- factor(cor_df$Var1, levels = unique(c(paste0("df_", 1:10), paste0("n_", 1:10))))
cor_df$Var2 <- factor(cor_df$Var2, levels = unique(c(paste0("df_", 1:10), paste0("n_", 1:10))))

p5 <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  facet_wrap(~ family, scales = "free") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0.8, limit = c(0, 1), name = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation of Subject Variance Across Complexities",
       x = "Complexity 1", y = "Complexity 2")
ggsave("subject_variance_correlation.png", p5, width = 10, height = 5)

# 6. AIC Comparison between ns and lspline --------------------------------
df_aic_wide <- df_metrics %>%
  select(gene, family, complexity, AIC) %>%
  pivot_wider(names_from = family, values_from = AIC) %>%
  mutate(delta_AIC = ns - lspline)

p6 <- ggplot(df_aic_wide, aes(x = factor(complexity), y = delta_AIC)) +
  geom_boxplot(outlier.size = 0.5, fill = "lightblue", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "AIC Comparison (ns vs lspline)",
       subtitle = "Negative Delta AIC indicates Natural Spline (ns) fits better",
       x = "Complexity (df/n)", y = "Delta AIC (ns - lspline)")
ggsave("aic_comparison_ns_vs_lspline.png", p6, width = 8, height = 5)

# 7. Global Variance Partitioning (Representative Model: ns_3) ----------------
if ("ns_3" %in% names(results)) {
  vp_ns3 <- as.data.frame(results[["ns_3"]]$vp) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "covariate", values_to = "variance")
  
  # Clean up covariate names for plotting
  vp_ns3$covariate <- gsub("ns\\(Age.months, df = 3\\)", "Age (Spline)", vp_ns3$covariate)
  
  p7 <- ggplot(vp_ns3, aes(x = reorder(covariate, variance, FUN = median), y = variance)) +
    geom_boxplot(fill = "steelblue", outlier.size = 0.5, alpha = 0.7) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Global Variance Partitioning (ns, df=3)",
         subtitle = "Helps identify redundant covariates (those with near-zero variance)",
         x = "Covariate", y = "Variance Explained (%)")
  ggsave("global_variance_partitioning.png", p7, width = 10, height = 7)
}

cat("Plots generated: aic_tournament.png, subject_variance_stability.png, age_inference_power.png, subject_variance_distribution.png, subject_variance_correlation.png, aic_comparison_ns_vs_lspline.png, global_variance_partitioning.png\n")
