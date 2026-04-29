library(variancePartition)
library(splines)
library(lspline)
library(tidyverse)
library(edgeR)
library(BiocParallel)
library(lme4)
library(matrixStats)

# Set working directory to project root relative to this script
setwd("../../")

# Parallel processing - switching to SerialParam for stability
register(SerialParam())

# 1. Load Data ------------------------------------------------------------
cat("Loading data...\n")
dge <- readRDS("data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

# 2. CBC Data Integration -------------------------------------------------
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

# Check for NAs in metadata for formula components
cols_to_check <- c("Subject.ID", "RNA.isolation.Batch", "Lib_prep_batches", "Year.Drawn", 
                  "sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent",
                  "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")

cat("Checking for NAs in metadata columns:\n")
na_counts <- colSums(is.na(meta[, cols_to_check]))
print(na_counts)

# Drop samples with any NAs in these columns
keep_final <- complete.cases(meta[, cols_to_check])
cat(sprintf("Dropping %d samples due to NAs\n", sum(!keep_final)))
meta <- meta[keep_final, ]
dge <- dge[, keep_final]

# Scale continuous covariates to prevent numeric instability and dropped variance components
cat("Scaling continuous covariates...\n")
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

# 3. Gene Selection (1000 HVGs for benchmark) --------------------------------------------
cat("Selecting 1000 HVGs...\n")
cpm_data <- cpm(dge, log = TRUE)
rv <- rowVars(cpm_data)
hvg_genes <- names(sort(rv, decreasing = TRUE)[1:1000])
dge_hvg <- dge[hvg_genes, ]
v_hvg <- voom(dge_hvg) # Basic voom for log-counts

# 4. Helper Functions -----------------------------------------------------

extract_spline_pvalues <- function(fit, term_pattern) {
  all_coefs <- colnames(coef(fit))
  spline_coefs <- all_coefs[grepl(term_pattern, all_coefs)]
  if (length(spline_coefs) == 0) {
    # Try exact match for linear case
    if ("Age.months" %in% all_coefs) spline_coefs <- "Age.months"
    else return(rep(NA, nrow(fit)))
  }
  res <- topTable(fit, coef = spline_coefs, number = Inf, sort.by = "none")
  pvals <- res[rownames(fit), "P.Value"]
  names(pvals) <- rownames(fit)
  return(pvals)
}

# Function to fit lmer and extract AIC/BIC for a single gene
fit_gene_metrics <- function(gene_expr, form, metadata, ctrl) {
  tryCatch({
    m <- lmer(as.formula(paste("gene_expr", form)), data = metadata, REML = FALSE, control = ctrl)
    return(c(AIC = AIC(m), BIC = BIC(m)))
  }, error = function(e) return(c(AIC = NA, BIC = NA)))
}

# 5. Modeling Tournament --------------------------------------------------
cat("Starting Modeling Tournament...\n")

# Use bobyqa optimizer to prevent variance collapse with spline matrices
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))

partial_file <- "review_comments_address/spline_degree_freedom/spline_benchmark_results_partial.rds"
if (file.exists(partial_file)) {
  cat("Loading partial results from", partial_file, "\n")
  results <- readRDS(partial_file)
} else {
  results <- list()
}
base_formula_rhs <- "~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt"

families <- list(
  ns = function(df) paste0("ns(Age.months, df=", df, ")"),
  lspline = function(n) if(n == 1) "Age.months" else paste0("elspline(Age.months, n=", n, ")")
)

for (fam_name in names(families)) {
  for (complexity in 1:10) {
    cat(sprintf("\n--- %s complexity %d ---\n", fam_name, complexity))
    
    res_name <- paste0(fam_name, "_", complexity)
    if (res_name %in% names(results)) {
      cat("Skipping", res_name, "as it already exists in partial results\n")
      next
    }
    
    spline_term <- families[[fam_name]](complexity)
    form_str <- paste(base_formula_rhs, "+", spline_term)
    form <- as.formula(form_str)
    
    # A. Variance Partition (via fitExtractVarPartModel)
    cat("Running fitExtractVarPartModel...\n")
    # Wrap in try to avoid halting
    vp <- tryCatch({
      vobj <- voomWithDreamWeights(dge_hvg$counts, form, meta, BPPARAM=SerialParam(), control=ctrl)
      fitExtractVarPartModel(vobj, form, meta, BPPARAM=SerialParam(), control=ctrl)
    }, error = function(e) { cat("Error in VarPart: ", e$message, "\n"); return(NULL) })
    
    if (is.null(vp)) next
    
    # B. P-values (via dream)
    cat("Running dream for P-values...\n")
    fit <- dream(vobj, form, meta, BPPARAM=SerialParam(), control=ctrl)
    
    term_pattern <- if(fam_name == "ns") "ns\\(Age.months" else if(complexity == 1) "Age.months" else "elspline\\(Age.months"
    age_pvals <- extract_spline_pvalues(fit, term_pattern)
    
    # C. AIC/BIC (via lmer)
    cat(sprintf("Calculating AIC/BIC for %d genes...\n", nrow(dge_hvg)))
    expr_mat <- v_hvg$E # log2-CPM from voom
    metrics_list <- lapply(1:nrow(expr_mat), function(i) {
      gene_name <- rownames(expr_mat)[i]
      if (i %% 10 == 0) cat(sprintf("  Gene %d/%d (%s)...\n", i, nrow(dge_hvg), gene_name))
      fit_gene_metrics(expr_mat[i, ], paste(form_str), meta, ctrl)
    })
    metrics_df <- do.call(rbind, metrics_list)
    rownames(metrics_df) <- rownames(expr_mat)
    
    # Store results
    results[[paste0(fam_name, "_", complexity)]] <- list(
      vp = vp,
      age_pvals = age_pvals,
      metrics = metrics_df,
      complexity = complexity,
      family = fam_name
    )
    
    # Save intermediate results
    saveRDS(results, "review_comments_address/spline_degree_freedom/spline_benchmark_results_partial.rds")
  }
}

# Save final results
saveRDS(results, "review_comments_address/spline_degree_freedom/spline_benchmark_results.rds")
cat("\nBenchmark complete. Results saved to review_comments_address/spline_degree_freedom/spline_benchmark_results.rds\n")
