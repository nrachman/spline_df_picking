library(variancePartition)
library(lme4)
library(lspline)
library(tidyverse)
library(edgeR)
library(matrixStats)
library(BiocParallel)

register(SerialParam())

# 1. Load and Prep Data ----------------------------------------------------
cat("Loading DGE data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
cat("DGE data loaded.\n")
# Use consistent ID format (underscore)
colnames(dge) <- gsub("\\.", "_", colnames(dge))
rownames(dge$samples) <- gsub("\\.", "_", rownames(dge$samples))

dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

# Create Join_Key for DGE (Prefix before _S)
dge_meta <- dge$samples %>% 
  rownames_to_column("DGE_Full_ID") %>%
  mutate(Join_Key = gsub("_S[0-9]+$", "", DGE_Full_ID))

# Load CBC
cat("Loading CBC data...\n")
cbc <- read_tsv("../../data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
cat("CBC data loaded.\n")
cbc <- cbc %>% 
  mutate(Sample_ID = gsub("\\.", "_", Sample_ID)) %>%
  mutate(Join_Key = gsub("_S[0-9]+$", "", Sample_ID))

# Load Deconvolution (LM22)
cat("Loading Deconvolution data...\n")
fardeep <- readRDS("../../data/analysis_out/fardeep/fardeep_lm22_batch_corrected_cpm_nolog_add_cbc_mat.rds")
cat("Deconvolution data loaded.\n")
lm22_cols <- colnames(fardeep)[8:29]
fardeep_df <- as.data.frame(fardeep) %>% 
  rownames_to_column("FD_Full_ID") %>%
  mutate(Join_Key = gsub("_S[0-9]+$", "", FD_Full_ID)) %>%
  select(Join_Key, all_of(lm22_cols))

# Merge
meta <- dge_meta %>%
  inner_join(cbc %>% select(Join_Key, Lymph.pcnt, Mid.pcnt, Gran.pcnt, MON.pcnt, EOS.pcnt), by = "Join_Key") %>%
  inner_join(fardeep_df, by = "Join_Key")

# Ensure unique and consistent indexing with the DGE object
rownames(meta) <- meta$DGE_Full_ID
dge_subset <- dge[, meta$DGE_Full_ID]

cat(sprintf("Using %d samples common to DGE, CBC, and Deconvolution.\n", nrow(meta)))

# Factors
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

# Scale all continuous covariates
cbc_vars <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Age.months",
                     cbc_vars, lm22_cols)

for (v in continuous_vars) {
  if (v %in% colnames(meta)) {
    meta[[v]] <- scale(as.numeric(meta[[v]]))[,1]
  }
}

# 2. Gene Selection (1000 HVGs) -------------------------------------------
cat("Selecting 1000 HVGs...\n")
cpm_data <- cpm(dge_subset, log = TRUE)
rv <- rowVars(cpm_data)
names(rv) <- rownames(cpm_data)
hvg_genes <- names(sort(rv, decreasing = TRUE)[1:1000])
expr_mat <- cpm_data[hvg_genes, ]

cat(sprintf("Expression matrix dimensions: %d x %d\n", nrow(expr_mat), ncol(expr_mat)))

# 3. Helper for Variance Partitioning ------------------------------------
calc_vp_aggregated <- function(model, age_terms, cell_terms) {
  re_vars_df <- as.data.frame(VarCorr(model))
  re_vars <- setNames(re_vars_df$vcov, re_vars_df$grp)
  var_resid <- sigma(model)^2
  
  X <- getME(model, "X")
  beta <- fixef(model)
  
  # Total fixed variance calculation
  X_no_int <- X[, colnames(X) != "(Intercept)", drop=FALSE]
  
  # Age Bucket
  age_cols <- which(colnames(X) %in% age_terms)
  var_age <- if(length(age_cols) > 0) var(as.numeric(X[, age_cols, drop=FALSE] %*% beta[age_cols])) else 0
  
  # Cell Bucket (Handling potential backticks in model colnames)
  cell_terms_clean <- paste0("`", cell_terms, "`")
  cell_cols <- which(colnames(X) %in% cell_terms | colnames(X) %in% cell_terms_clean)
  var_cells <- if(length(cell_cols) > 0) var(as.numeric(X[, cell_cols, drop=FALSE] %*% beta[cell_cols])) else 0
  
  # Subject Bucket
  var_subject <- if("Subject.ID" %in% names(re_vars)) re_vars[["Subject.ID"]] else 0
  
  # Other components (Batch, etc.)
  re_others <- re_vars[!(names(re_vars) %in% c("Subject.ID", "Residual"))]
  var_re_others <- sum(unlist(re_others))
  
  fixed_others_names <- setdiff(colnames(X_no_int), c(age_terms, cell_terms, cell_terms_clean))
  fixed_others_cols <- which(colnames(X) %in% fixed_others_names)
  var_fixed_others <- if(length(fixed_others_cols) > 0) var(as.numeric(X[, fixed_others_cols, drop=FALSE] %*% beta[fixed_others_cols])) else 0
  
  vars <- c(Age = var_age, Subject = var_subject, CellFreq = var_cells, 
            Technical = var_re_others + var_fixed_others, Residuals = var_resid)
  total <- sum(vars)
  return(vars / total)
}

# 4. Benchmarking ---------------------------------------------------------
cat("Starting Benchmark...\n")
ctrl <- lmerControl(optimizer = "bobyqa")
age_terms_ls4 <- c("elspline(Age.months, n = 4)1", "elspline(Age.months, n = 4)2", "elspline(Age.months, n = 4)3", "elspline(Age.months, n = 4)4")

# A. CBC-based Model
cat("Running CBC-based model...\n")
form_cbc_str <- paste("y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + elspline(Age.months, n=4) +", paste(cbc_vars, collapse=" + "))

vp_cbc_list <- list()
for(i in 1:nrow(expr_mat)) {
  if(i %% 200 == 0) cat(sprintf("  Gene %d/1000...\n", i))
  y <- expr_mat[i, ]
  tryCatch({
    m <- suppressMessages(suppressWarnings(lmer(as.formula(form_cbc_str), data = meta, REML = FALSE, control = ctrl)))
    vp_cbc_list[[i]] <- calc_vp_aggregated(m, age_terms_ls4, cbc_vars)
  }, error = function(e) { vp_cbc_list[[i]] <<- rep(NA, 5) })
}
vp_cbc_df <- do.call(rbind, vp_cbc_list)
rownames(vp_cbc_df) <- rownames(expr_mat)

# B. Deconvolution-based Model
cat("Running Deconvolution-based model...\n")
# Backtick deconv variables for formula
lm22_vars_bt <- paste0("`", lm22_cols, "`")
form_deconv_str <- paste("y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + elspline(Age.months, n=4) +", paste(lm22_vars_bt, collapse=" + "))

vp_deconv_list <- list()
for(i in 1:nrow(expr_mat)) {
  if(i %% 200 == 0) cat(sprintf("  Gene %d/1000...\n", i))
  y <- expr_mat[i, ]
  tryCatch({
    m <- suppressMessages(suppressWarnings(lmer(as.formula(form_deconv_str), data = meta, REML = FALSE, control = ctrl)))
    vp_deconv_list[[i]] <- calc_vp_aggregated(m, age_terms_ls4, lm22_cols)
  }, error = function(e) { vp_deconv_list[[i]] <<- rep(NA, 5) })
}
vp_deconv_df <- do.call(rbind, vp_deconv_list)
rownames(vp_deconv_df) <- rownames(expr_mat)

# Save
results <- list(vp_cbc = vp_cbc_df, vp_deconv = vp_deconv_df)
saveRDS(results, "output/rds/cell_freq_comparison_results.rds")
cat("Done. Results saved to output/rds/cell_freq_comparison_results.rds\n")
