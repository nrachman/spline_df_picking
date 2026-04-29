library(variancePartition)
library(lme4)
library(splines)
library(lspline)
library(tidyverse)
library(edgeR)
library(BiocParallel)
library(matrixStats)

# Set working directory to project root relative to this script
setwd("../../")

# Parallel processing
register(SerialParam())

# 1. Load Data ------------------------------------------------------------
cat("Loading data...\n")
dge <- readRDS("data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

cat("Loading CBC data...\n")
cbc <- read_tsv("data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
cbc <- cbc %>% mutate(Sample_ID = gsub("\\.", "_", Sample_ID))

keep_samples <- intersect(dge$samples$Sample.name, cbc$Sample_ID)
dge <- dge[, dge$samples$Sample.name %in% keep_samples]

keep_cells <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt")
meta <- dge$samples %>%
  left_join(cbc %>% select(Sample_ID, all_of(keep_cells)), by = c("Sample.name" = "Sample_ID"))

# Set rownames to match colnames(dge)
rownames(meta) <- colnames(dge)

# Factors and numeric conversion
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

cols_to_check <- c("Subject.ID", "RNA.isolation.Batch", "Lib_prep_batches", "Year.Drawn", 
                  "sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent",
                  "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")

# Drop samples with any NAs
keep_final <- complete.cases(meta[, cols_to_check])
meta <- meta[keep_final, ]
dge <- dge[, keep_final]

# Scale continuous covariates (except Age.months which we might bin or use in splines)
cat("Scaling continuous covariates...\n")
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

# 2. Advanced Age Transformations ----------------------------------------

# A. Age Binned (6 samples per bin)
cat("Creating Age Bins (6 samples per bin)...\n")
n_bins <- floor(nrow(meta) / 6)
meta <- meta %>%
  arrange(Age.months) %>%
  mutate(Age_Binned = factor(rep(1:n_bins, each = 6, length.out = n())))

# B. Age Years Categorical
cat("Creating Age Years Categorical...\n")
meta$Age_Years_Cat <- factor(round(meta$Age.months / 12, 0))

# 3. Gene Selection (1000 HVGs) -------------------------------------------
cat("Selecting 1000 HVGs...\n")
cpm_data <- cpm(dge, log = TRUE)
rv <- rowVars(cpm_data)
hvg_genes <- names(sort(rv, decreasing = TRUE)[1:1000])
dge_hvg <- dge[hvg_genes, ]

# 4. Helper Function: Realized Predictions Variance Partitioning ---------
calc_realized_vp <- function(model, age_term_names) {
  # Get design matrices
  X <- getME(model, "X")
  beta <- fixef(model)
  Z <- getME(model, "Z")
  u <- getME(model, "b")
  
  # 1. Total Age Effect (Combined spline columns)
  age_cols <- which(colnames(X) %in% age_term_names)
  if (length(age_cols) == 0) {
     # Check random effects if Age was modeled as random
     age_cols <- character(0)
  }
  
  # We'll do this more generally: 
  # Variance of each component's realized prediction
  
  # Fixed Effects Components
  # Group fixed effects if they belong to the same concept (like Spline)
  # For this benchmark, we'll focus on Age, Subject, and "Other Fixed"
  
  # Age Prediction
  if (length(age_cols) > 0) {
    pred_age <- as.numeric(X[, age_cols, drop=FALSE] %*% beta[age_cols])
  } else {
    pred_age <- rep(0, nrow(X))
  }
  
  # Random Effects Components
  gp <- getME(model, "Gp")
  re_names <- names(ranef(model))
  
  re_vars <- list()
  for (i in seq_along(re_names)) {
    start <- gp[i] + 1
    end <- gp[i+1]
    Z_sub <- Z[, start:end, drop=FALSE]
    u_sub <- u[start:end]
    pred_re <- as.numeric(Z_sub %*% u_sub)
    
    # If Age is a random effect (Age_Binned or Age_Years_Cat)
    if (re_names[i] %in% c("Age_Binned", "Age_Years_Cat")) {
      pred_age <- pred_age + pred_re
    } else {
      re_vars[[re_names[i]]] <- var(pred_re)
    }
  }
  
  var_age <- var(pred_age)
  var_resid <- sigma(model)^2
  
  vars <- c(Age = var_age, unlist(re_vars), Residuals = var_resid)
  total <- sum(vars)
  return(vars / total)
}

# 5. Benchmarking ---------------------------------------------------------
cat("Starting Advanced Benchmarking...\n")

# Base Formula RHS
base_covariates <- "sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn)"

# Models to Test
models_to_test <- list(
  "ns_df4" = list(
    formula = paste("~", base_covariates, "+ (1|Subject.ID) + ns(Age.months, df=4)"),
    age_terms = c("ns(Age.months, df = 4)1", "ns(Age.months, df = 4)2", "ns(Age.months, df = 4)3", "ns(Age.months, df = 4)4")
  ),
  "lspline_n4" = list(
    formula = paste("~", base_covariates, "+ (1|Subject.ID) + elspline(Age.months, n=4)"),
    age_terms = c("elspline(Age.months, n = 4)1", "elspline(Age.months, n = 4)2", "elspline(Age.months, n = 4)3", "elspline(Age.months, n = 4)4")
  ),
  "age_random_binned" = list(
    formula = paste("~", base_covariates, "+ (1|Subject.ID) + (1|Age_Binned)"),
    age_terms = character(0)
  ),
  "age_fixed_cat" = list(
    formula = paste("~", base_covariates, "+ (1|Subject.ID) + Age_Years_Cat"),
    age_terms = paste0("Age_Years_Cat", levels(meta$Age_Years_Cat)[-1])
  ),
  "age_random_cat" = list(
    formula = paste("~", base_covariates, "+ (1|Subject.ID) + (1|Age_Years_Cat)"),
    age_terms = character(0)
  )
)

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

advanced_results <- list()

for (mod_name in names(models_to_test)) {
  cat(sprintf("\n--- Processing Model: %s ---\n", mod_name))
  spec <- models_to_test[[mod_name]]
  form <- as.formula(spec$formula)
  
  # voomWithDreamWeights
  vobj <- voomWithDreamWeights(dge_hvg$counts, form, meta, BPPARAM=SerialParam(), control=ctrl)
  
  # A. Standard Variance Partitioning
  cat("  Calculating Standard Variance Partitioning...\n")
  vp_std <- fitExtractVarPartModel(vobj, form, meta, BPPARAM=SerialParam(), control=ctrl)
  
  # B. Realized Predictions Variance Partitioning
  cat("  Calculating Realized Predictions Variance Partitioning...\n")
  expr_mat <- vobj$E
  realized_vp_list <- lapply(1:nrow(expr_mat), function(i) {
    if (i %% 100 == 0) cat(sprintf("    Gene %d/1000...\n", i))
    tryCatch({
      m <- lmer(as.formula(paste("expr_mat[i, ]", spec$formula)), data = meta, REML = FALSE, control = ctrl)
      calc_realized_vp(m, spec$age_terms)
    }, error = function(e) return(NULL))
  })
  
  # Filter out nulls and combine
  valid_indices <- which(!sapply(realized_vp_list, is.null))
  realized_vp_df <- do.call(rbind, realized_vp_list[valid_indices])
  rownames(realized_vp_df) <- rownames(expr_mat)[valid_indices]
  
  # C. AIC/BIC
  cat("  Calculating AIC/BIC...\n")
  metrics_list <- lapply(1:nrow(expr_mat), function(i) {
    tryCatch({
      m <- lmer(as.formula(paste("expr_mat[i, ]", spec$formula)), data = meta, REML = FALSE, control = ctrl)
      c(AIC = AIC(m), BIC = BIC(m))
    }, error = function(e) return(c(AIC = NA, BIC = NA)))
  })
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- rownames(expr_mat)
  
  advanced_results[[mod_name]] <- list(
    vp_std = vp_std,
    vp_realized = realized_vp_df,
    metrics = metrics_df
  )
  
  # Save intermediate
  saveRDS(advanced_results, "output/rds/advanced_variance_benchmarks.rds")
}

cat("\nAdvanced Benchmarking Complete. Results saved to output/rds/advanced_variance_benchmarks.rds\n")
