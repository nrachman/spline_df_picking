library(variancePartition)
library(lme4)
library(splines)
library(lspline)
library(tidyverse)
library(edgeR)

# 1. Load Data (Synced) ----------------------------------------------------
cat("Loading and syncing data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

cbc <- read_tsv("../../data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
cbc <- cbc %>% mutate(Sample_ID = gsub("\\.", "_", Sample_ID))

keep_samples <- intersect(dge$samples$Sample.name, cbc$Sample_ID)
dge <- dge[, dge$samples$Sample.name %in% keep_samples]

keep_cells <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt")
meta <- dge$samples %>%
  left_join(cbc %>% select(Sample_ID, all_of(keep_cells)), by = c("Sample.name" = "Sample_ID"))

rownames(meta) <- colnames(dge)

# Factors and numeric conversion
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)
meta$Age_Years_Cat <- factor(round(meta$Age.months / 12, 0))

# Sync order
meta <- meta %>% arrange(Age.months)
dge <- dge[, rownames(meta)]

# Scale continuous covariates (EXCEPT Age.months for plotting raw scale)
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

# Normalized expression (log2-CPM)
expr_mat <- cpm(dge, log = TRUE)

# 2. Select Genes ---------------------------------------------------------
# Based on common age-related genes or previous benchmark results
# We'll pick some that are likely to show diverse patterns
results <- readRDS("output/rds/spline_benchmark_results.rds")
res_ns4 <- results[["ns_4"]]
pvals <- p.adjust(res_ns4$age_pvals, method = "fdr")
candidate_genes <- names(sort(pvals))[1:50]

# Pick 5 interesting ones (some high signal, some maybe less so but non-linear)
target_genes <- c(candidate_genes[1], candidate_genes[5], candidate_genes[10], candidate_genes[15], candidate_genes[20])
# Ensure they are in expr_mat
target_genes <- target_genes[target_genes %in% rownames(expr_mat)]

cat("Plotting genes:", paste(target_genes, collapse=", "), "\n")

# 3. Fit Models and Plot --------------------------------------------------
dir.create("output/images/gene_examples", showWarnings = FALSE, recursive = TRUE)

# Create prediction data
age_range <- seq(min(meta$Age.months), max(meta$Age.months), length.out = 100)
pred_df_base <- data.frame(
  Age.months = age_range,
  sex.numeric = 0, RIN = 0, mk_dup.PERCENT_DUPLICATION = 0, 
  star.uniquely_mapped_percent = 0, Lymph.pcnt = 0, Mid.pcnt = 0, 
  Gran.pcnt = 0, MON.pcnt = 0, EOS.pcnt = 0,
  RNA.isolation.Batch = levels(meta$RNA.isolation.Batch)[1],
  Lib_prep_batches = levels(meta$Lib_prep_batches)[1],
  Year.Drawn = levels(meta$Year.Drawn)[1],
  Subject.ID = levels(meta$Subject.ID)[1] # Note: we will use population-level prediction
)
pred_df_base$Age_Years_Cat <- factor(round(pred_df_base$Age.months / 12, 0), levels = levels(meta$Age_Years_Cat))

ctrl <- lmerControl(optimizer = "bobyqa")

for (gene in target_genes) {
  cat("  Processing", gene, "...\n")
  y <- expr_mat[gene, ]
  
  # ns(df=4)
  m_ns <- lmer(y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + ns(Age.months, df=4) + 
                sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
                Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt, 
              data = meta, control = ctrl)
  
  # lspline(n=4)
  m_ls <- lmer(y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + elspline(Age.months, n=4) + 
                sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
                Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt, 
              data = meta, control = ctrl)
  
  # age_fixed_cat
  m_cat <- lmer(y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + Age_Years_Cat + 
                sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
                Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt, 
              data = meta, control = ctrl)
  
  # Predictions (Fixing all other covariates to mean=0)
  p_ns <- predict(m_ns, newdata = pred_df_base, re.form = NA)
  p_ls <- predict(m_ls, newdata = pred_df_base, re.form = NA)
  p_cat <- predict(m_cat, newdata = pred_df_base, re.form = NA)
  
  preds <- data.frame(
    Age.months = age_range,
    ns_df4 = p_ns,
    lspline_n4 = p_ls,
    Categorical = p_cat
  ) %>% pivot_longer(-Age.months, names_to = "Model", values_to = "Expression")
  
  # Adjust raw data for visualization (subtracting other fixed effects except intercept)
  # This helps see the "age-only" trend in the raw points
  X <- getME(m_ns, "X")
  beta <- fixef(m_ns)
  age_terms <- grep("ns\\(Age.months", colnames(X), value = TRUE)
  non_age_fixed <- setdiff(colnames(X), c("(Intercept)", age_terms))
  
  y_adj <- y - as.matrix(meta[, c("sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")]) %*% beta[c("sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")]
  
  plot_data <- data.frame(Age.months = meta$Age.months, Expression = as.numeric(y_adj))
  
  p <- ggplot() +
    geom_point(data = plot_data, aes(x = Age.months, y = Expression), alpha = 0.3, size = 1) +
    geom_line(data = preds, aes(x = Age.months, y = Expression, color = Model, linetype = Model), size = 1) +
    theme_minimal() +
    labs(title = paste("Age Trajectory Comparison:", gene),
         subtitle = "Expression adjusted for sex, RIN, technical factors, and cell counts",
         x = "Age (months)", y = "Adjusted Expression (log2-CPM)") +
    scale_color_manual(values = c("ns_df4" = "#E41A1C", "lspline_n4" = "#377EB8", "Categorical" = "#4DAF4A"))
  
  ggsave(sprintf("output/images/gene_examples/age_trend_%s.png", gene), p, width = 8, height = 5)
}
