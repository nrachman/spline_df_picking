library(variancePartition)
library(splines)
library(lspline)
library(tidyverse)
library(edgeR)
library(BiocParallel)
library(lme4)
library(matrixStats)

# 1. Load Data ------------------------------------------------------------
options(warn = -1)
cat("Loading data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

# 2. CBC Data Integration -------------------------------------------------
cat("Loading CBC data...\n")
cbc <- read_tsv("../../data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
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

# Scale continuous covariates
cat("Scaling continuous covariates...\n")
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

# 3. Gene Selection (1000 Random Genes for benchmark) --------------------------------------------
cat("Selecting 1000 Random Genes...\n")
cpm_data <- cpm(dge, log = TRUE)
rv <- rowVars(cpm_data)
# Exclude the top 2000 HVGs to ensure these are truly "random/average" genes
hvg_excluded <- names(sort(rv, decreasing = TRUE)[1:2000])
set.seed(123)
rand_genes <- sample(setdiff(rownames(dge), hvg_excluded), 500)
dge_subset <- dge[rand_genes, ]

# 4. Helper Functions -----------------------------------------------------

extract_spline_pvalues <- function(fit, term_pattern) {
  all_coefs <- colnames(coef(fit))
  spline_coefs <- all_coefs[grepl(term_pattern, all_coefs)]
  if (length(spline_coefs) == 0) {
    if ("Age.months" %in% all_coefs) spline_coefs <- "Age.months"
    else return(rep(NA, nrow(fit)))
  }
  res <- topTable(fit, coef = spline_coefs, number = Inf, sort.by = "none")
  pvals <- res[rownames(fit), "P.Value"]
  names(pvals) <- rownames(fit)
  return(pvals)
}

fit_gene_metrics <- function(gene_expr, form, metadata, ctrl) {
  tryCatch({
    m <- lmer(as.formula(paste("gene_expr", form)), data = metadata, REML = FALSE, control = ctrl)
    return(c(AIC = AIC(m), BIC = BIC(m)))
  }, error = function(e) return(c(AIC = NA, BIC = NA)))
}

# 5. Modeling Tournament --------------------------------------------------
cat("Starting Modeling Tournament...\n")

# Use bobyqa optimizer and ignore singular fits to suppress chatter
ctrl <- lmerControl(optimizer = "bobyqa", 
                    optCtrl = list(maxfun = 100000),
                    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))

dir.create("output/rds", showWarnings = FALSE, recursive = TRUE)
partial_file <- "output/rds/spline_benchmark_results_random_partial.rds"
if (file.exists(partial_file)) {
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
    
    # A. Variance Partition
    cat("Running fitExtractVarPartModel...\n")
    vobj <- tryCatch({
      suppressMessages(suppressWarnings(
        voomWithDreamWeights(dge_subset$counts, form, meta, BPPARAM=SerialParam(), control=ctrl)
      ))
    }, error = function(e) { cat("Error in voom: ", e$message, "\n"); return(NULL) })
    
    if (is.null(vobj)) next
    
    vp <- tryCatch({
      suppressMessages(suppressWarnings(
        fitExtractVarPartModel(vobj, form, meta, BPPARAM=SerialParam(), control=ctrl)
      ))
    }, error = function(e) { cat("Error in VarPart: ", e$message, "\n"); return(NULL) })
    
    if (is.null(vp)) next
    
    # B. P-values
    cat("Running dream for P-values...\n")
    fit <- suppressMessages(suppressWarnings(
      dream(vobj, form, meta, BPPARAM=SerialParam(), control=ctrl)
    ))
    term_pattern <- if(fam_name == "ns") "ns\\(Age.months" else if(complexity == 1) "Age.months" else "elspline\\(Age.months"
    age_pvals <- extract_spline_pvalues(fit, term_pattern)
    
    # C. AIC/BIC
    cat("Calculating AIC/BIC...\n")
    metrics_list <- bplapply(1:nrow(vobj), function(i) {
      fit_gene_metrics(vobj$E[i,], form_str, meta, ctrl)
    }, BPPARAM = SerialParam())
    metrics_mat <- do.call(rbind, metrics_list)
    rownames(metrics_mat) <- rownames(vobj)
    
    results[[res_name]] <- list(
      metrics = metrics_mat,
      vp = vp,
      age_pvals = age_pvals,
      complexity = complexity,
      family = fam_name
    )
    
    saveRDS(results, partial_file)
  }
}

saveRDS(results, "output/rds/spline_benchmark_results_random.rds")
cat("\nBenchmark complete. Results saved to output/rds/spline_benchmark_results_random.rds\n")
