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
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)
meta$Age_Years_Cat <- factor(round(meta$Age.months / 12, 0))

meta <- meta %>% arrange(Age.months)
dge <- dge[, rownames(meta)]

continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

expr_mat <- cpm(dge, log = TRUE)

# 2. Target Genes (Identified in previous step) ----------------------------
target_genes <- c("HBG2", "IGF2BP3", "MID2", "CFH", "IGFBP3", "CD177", "NAP1L2", "TMTC1", "MPP2", "NEFL")

# 3. Fit Models and Plot --------------------------------------------------
dir.create("output/images/model_differences", showWarnings = FALSE, recursive = TRUE)

age_range <- seq(min(meta$Age.months), max(meta$Age.months), length.out = 100)
pred_df_base <- data.frame(
  Age.months = age_range,
  sex.numeric = 0, RIN = 0, mk_dup.PERCENT_DUPLICATION = 0, 
  star.uniquely_mapped_percent = 0, Lymph.pcnt = 0, Mid.pcnt = 0, 
  Gran.pcnt = 0, MON.pcnt = 0, EOS.pcnt = 0,
  RNA.isolation.Batch = levels(meta$RNA.isolation.Batch)[1],
  Lib_prep_batches = levels(meta$Lib_prep_batches)[1],
  Year.Drawn = levels(meta$Year.Drawn)[1],
  Subject.ID = levels(meta$Subject.ID)[1]
)
pred_df_base$Age_Years_Cat <- factor(round(pred_df_base$Age.months / 12, 0), levels = levels(meta$Age_Years_Cat))

ctrl <- lmerControl(optimizer = "bobyqa")

# Get p-values for labeling
results <- readRDS("output/rds/spline_benchmark_results.rds")
ns4_p <- results[["ns_4"]]$age_pvals
ls4_p <- results[["lspline_4"]]$age_pvals

for (gene in target_genes) {
  if (!(gene %in% rownames(expr_mat))) next
  cat("  Processing", gene, "...\n")
  y <- expr_mat[gene, ]
  
  m_ns <- lmer(y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + ns(Age.months, df=4) + 
                sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
                Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt, 
              data = meta, control = ctrl)
  
  m_ls <- lmer(y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + elspline(Age.months, n=4) + 
                sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
                Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt, 
              data = meta, control = ctrl)
  
  p_ns_pred <- predict(m_ns, newdata = pred_df_base, re.form = NA)
  p_ls_pred <- predict(m_ls, newdata = pred_df_base, re.form = NA)
  
  preds <- data.frame(
    Age.months = age_range,
    ns_df4 = p_ns_pred,
    lspline_n4 = p_ls_pred
  ) %>% pivot_longer(-Age.months, names_to = "Model", values_to = "Expression")
  
  # Adjusted raw data
  beta <- fixef(m_ns)
  y_adj <- y - as.matrix(meta[, c("sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")]) %*% beta[c("sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")]
  plot_data <- data.frame(Age.months = meta$Age.months, Expression = as.numeric(y_adj))
  
  p <- ggplot() +
    geom_point(data = plot_data, aes(x = Age.months, y = Expression), alpha = 0.3, size = 1) +
    geom_line(data = preds, aes(x = Age.months, y = Expression, color = Model, linetype = Model), size = 1) +
    theme_minimal() +
    labs(title = paste("Model Divergence:", gene),
         subtitle = sprintf("ns p=%.2e | lspline p=%.2e\nAdjusted for technical factors and cell counts", ns4_p[gene], ls4_p[gene]),
         x = "Age (months)", y = "Adjusted Expression (log2-CPM)") +
    scale_color_manual(values = c("ns_df4" = "#E41A1C", "lspline_n4" = "#377EB8"))
  
  ggsave(sprintf("output/images/model_differences/diff_%s.png", gene), p, width = 8, height = 5)
}
