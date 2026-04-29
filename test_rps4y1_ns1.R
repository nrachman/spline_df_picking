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
keep_final <- complete.cases(meta[, c("Subject.ID", "RNA.isolation.Batch", "Lib_prep_batches", "Year.Drawn", "sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")])
meta <- meta[keep_final, ]

cpm_data <- edgeR::cpm(dge, log=TRUE)
expr <- cpm_data["RPS4Y1", keep_final]
meta$expr <- as.numeric(expr)

form_ns <- expr ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + ns(Age.months, df=1)

cat("\nUNSCALED Age.months:\n")
m_unscaled <- lmer(form_ns, data=meta)
print(VarCorr(m_unscaled))

cat("\nSCALED Age.months:\n")
meta$Age.months <- scale(meta$Age.months)[,1]
m_scaled <- lmer(form_ns, data=meta)
print(VarCorr(m_scaled))
