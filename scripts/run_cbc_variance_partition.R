library(variancePartition)
library(lme4)
library(lspline)
library(tidyverse)
library(edgeR)
library(matrixStats)
library(BiocParallel)

register(SerialParam())

# 1. Load and Prep Data ----------------------------------------------------
cat("Loading and syncing data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
colnames(dge) <- gsub("\\.", "_", colnames(dge))
rownames(dge$samples) <- gsub("\\.", "_", rownames(dge$samples))

dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

# Create Join_Key for DGE (Prefix before _S)
dge_meta <- dge$samples %>% 
  rownames_to_column("DGE_Full_ID") %>%
  mutate(Join_Key = gsub("_S[0-9]+$", "", DGE_Full_ID))

# Load CBC
cbc <- read_tsv("../../data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
cbc <- cbc %>% 
  mutate(Sample_ID = gsub("\\.", "_", Sample_ID)) %>%
  mutate(Join_Key = gsub("_S[0-9]+$", "", Sample_ID))
cbc_vars <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")

# Merge
meta <- dge_meta %>%
  inner_join(cbc %>% select(Join_Key, all_of(cbc_vars)), by = "Join_Key")

# Ensure unique and consistent indexing
rownames(meta) <- meta$DGE_Full_ID
dge_subset <- dge[, meta$DGE_Full_ID]

cat(sprintf("Using %d samples common to DGE and CBC.\n", nrow(meta)))

# Factors
meta$Subject.ID <- factor(meta$Subject.ID)
meta$sex.numeric <- as.numeric(factor(meta$sex))
meta$Year.Drawn <- factor(meta$Year.Drawn)
meta$RNA.isolation.Batch <- factor(meta$RNA.isolation.Batch)
meta$Lib_prep_batches <- factor(meta$Lib_prep_batches)

# Scale all continuous covariates
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Age.months", cbc_vars)
for (v in continuous_vars) meta[[v]] <- scale(as.numeric(meta[[v]]))[,1]

# 2. Gene Selection (1000 HVGs) -------------------------------------------
cat("Selecting 1000 HVGs...\n")
cpm_data <- cpm(dge_subset, log = TRUE)
rv <- rowVars(cpm_data)
names(rv) <- rownames(cpm_data)
hvg_genes <- names(sort(rv, decreasing = TRUE)[1:1000])
expr_mat <- cpm_data[hvg_genes, ]

# 3. Helper for Aggregated Variance Partitioning -------------------------
calc_vp_aggregated <- function(model, age_terms, cell_terms) {
  re_vars_df <- as.data.frame(VarCorr(model))
  re_vars <- setNames(re_vars_df$vcov, re_vars_df$grp)
  var_resid <- sigma(model)^2
  
  X <- getME(model, "X")
  beta <- fixef(model)
  X_no_int <- X[, colnames(X) != "(Intercept)", drop=FALSE]
  
  # Age Bucket
  age_cols <- which(colnames(X) %in% age_terms)
  var_age <- if(length(age_cols) > 0) var(as.numeric(X[, age_cols, drop=FALSE] %*% beta[age_cols])) else 0
  
  # Cell Bucket
  cell_cols <- which(colnames(X) %in% cell_terms)
  var_cells <- if(length(cell_cols) > 0) var(as.numeric(X[, cell_cols, drop=FALSE] %*% beta[cell_cols])) else 0
  
  # Subject Bucket
  var_subject <- if("Subject.ID" %in% names(re_vars)) re_vars[["Subject.ID"]] else 0
  
  # Technical
  re_others <- sum(unlist(re_vars[!(names(re_vars) %in% c("Subject.ID", "Residual"))]))
  fixed_others_names <- setdiff(colnames(X_no_int), c(age_terms, cell_terms))
  fixed_others_cols <- which(colnames(X) %in% fixed_others_names)
  var_fixed_others <- if(length(fixed_others_cols) > 0) var(as.numeric(X[, fixed_others_cols, drop=FALSE] %*% beta[fixed_others_cols])) else 0
  
  vars <- c(Age = var_age, Subject = var_subject, CellFreq = var_cells, 
            Technical = re_others + var_fixed_others, Residuals = var_resid)
  return(vars / sum(vars))
}

# 4. Modeling ------------------------------------------------------------
cat("Starting Variance Partitioning (lspline n=4)...\n")
ctrl <- lmerControl(optimizer = "bobyqa")
age_terms <- c("elspline(Age.months, n = 4)1", "elspline(Age.months, n = 4)2", "elspline(Age.months, n = 4)3", "elspline(Age.months, n = 4)4")

form_str <- paste("y ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + elspline(Age.months, n=4) +", 
                  paste(cbc_vars, collapse=" + "))

vp_list <- list()
for(i in 1:nrow(expr_mat)) {
  if(i %% 200 == 0) cat(sprintf("  Gene %d/1000...\n", i))
  y <- expr_mat[i, ]
  tryCatch({
    m <- suppressMessages(suppressWarnings(lmer(as.formula(form_str), data = meta, REML = FALSE, control = ctrl)))
    vp_list[[i]] <- calc_vp_aggregated(m, age_terms, cbc_vars)
  }, error = function(e) { vp_list[[i]] <<- rep(NA, 5) })
}

vp_df <- do.call(rbind, vp_list)
rownames(vp_df) <- rownames(expr_mat)

dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)
saveRDS(vp_df, "output/rds/cbc_full_vp_results.rds")
cat("Results saved to output/rds/cbc_full_vp_results.rds\n")
