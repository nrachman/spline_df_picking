library(variancePartition)
library(lme4)
library(splines)
library(tidyverse)

setwd("../../")
dge <- readRDS("data/processed/limma_dgelist/dgelist.rds")
dge <- dge[, dge$samples$Sample.type == "case" & !is.na(dge$samples$RIN)]

cbc <- read_tsv("data/inputs/cbc/cbc_cleaned.tsv", show_col_types = FALSE)
cbc <- cbc %>% mutate(Sample_ID = gsub("\\.", "_", Sample_ID))
keep_samples <- intersect(dge$samples$Sample.name, cbc$Sample_ID)
dge <- dge[, dge$samples$Sample.name %in% keep_samples]
keep_cells <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt")
meta <- dge$samples %>% left_join(cbc %>% select(Sample_ID, all_of(keep_cells)), by = c("Sample.name" = "Sample_ID"))

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

cat("Scaling continuous covariates...\n")
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

expr <- as.numeric(edgeR::cpm(dge, log=TRUE)[1, keep_final])
form_base <- expr ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt

m_ns1 <- lmer(update(form_base, . ~ . + ns(Age.months, df=1)), data=meta)
m_ls1 <- lmer(update(form_base, . ~ . + Age.months), data=meta)

print("ns_1 VarCorr:")
print(VarCorr(m_ns1))
print("lspline_1 VarCorr:")
print(VarCorr(m_ls1))

expr_mat <- matrix(expr, nrow=1)
rownames(expr_mat) <- "gene1"

vp_ns <- fitExtractVarPartModel(expr_mat, update(form_base, . ~ . + ns(Age.months, df=1)), meta, showWarnings=FALSE)
vp_ls <- fitExtractVarPartModel(expr_mat, update(form_base, . ~ . + Age.months), meta, showWarnings=FALSE)

print("variancePartition ns_1:")
print(vp_ns)
print("variancePartition lspline_1:")
print(vp_ls)
