library(variancePartition)
library(lme4)
library(splines)
library(lspline)
library(tidyverse)
library(edgeR)
library(matrixStats)

# 1. Load Data ------------------------------------------------------------
cat("Loading data...\n")
dge <- readRDS("../../data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

cat("Loading CBC data...\n")
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

cols_to_check <- c("Subject.ID", "RNA.isolation.Batch", "Lib_prep_batches", "Year.Drawn", 
                  "sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent",
                  "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")

keep_final <- complete.cases(meta[, cols_to_check])
meta <- meta[keep_final, ]
dge <- dge[, keep_final]

continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

# Reordering meta by Age.months (SCRAMBLING)
cat("Reordering meta by Age.months (scrambling relationship with dge columns)...\n")
meta <- meta %>% arrange(Age.months)
# NOTE: dge is NOT reordered here, which matches the bug in run_advanced_variance_benchmarks.R

# Select top HVG
cpm_data <- cpm(dge, log = TRUE)
rv <- rowVars(cpm_data)
hvg_gene <- names(sort(rv, decreasing = TRUE)[1])
expr <- cpm_data[hvg_gene, ]

# Fit model with scrambled data
cat("Fitting model with scrambled data for gene:", hvg_gene, "\n")
base_covariates <- "sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn)"
form_str <- paste("expr ~", base_covariates, "+ (1|Subject.ID) + ns(Age.months, df=4)")
m_scrambled <- lmer(as.formula(form_str), data = meta, REML = FALSE)
cat("Scrambled Subject Variance:\n")
print(VarCorr(m_scrambled))

# Now try CORRECT order
cat("\nRe-syncing dge columns to meta order...\n")
expr_ordered <- expr[rownames(meta)]
m_ordered <- lmer(as.formula(paste("expr_ordered ~", base_covariates, "+ (1|Subject.ID) + ns(Age.months, df=4)")), data = meta, REML = FALSE)
cat("Ordered Subject Variance:\n")
print(VarCorr(m_ordered))
