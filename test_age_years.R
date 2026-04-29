library(variancePartition)
library(lme4)
library(splines)
library(lspline)
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
dge <- dge[, keep_final]
rownames(meta) <- colnames(dge)

# Convert Age to Years
meta$Age.years <- meta$Age.months / 12.0

# ONLY scale the other continuous covariates
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]

counts <- dge$counts[c("RPS4Y1", "KDM5D", "XIST"), ]

form_ns <- ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + ns(Age.years, df=1)

form_ls <- ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + Age.years

vobj_ns <- voomWithDreamWeights(counts, form_ns, meta, BPPARAM=SerialParam())
vp_ns <- fitExtractVarPartModel(vobj_ns, form_ns, meta, BPPARAM=SerialParam())

vobj_ls <- voomWithDreamWeights(counts, form_ls, meta, BPPARAM=SerialParam())
vp_ls <- fitExtractVarPartModel(vobj_ls, form_ls, meta, BPPARAM=SerialParam())

print("ns_1 (Age.years) Subject.ID Variance:")
print(vp_ns$Subject.ID)
print("lspline_1 (Age.years) Subject.ID Variance:")
print(vp_ls$Subject.ID)
