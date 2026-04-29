library(variancePartition)
library(lme4)
library(edgeR)

setwd("../../")
dge <- readRDS("data/processed/limma_dgelist/dgelist.rds")
meta <- dge$samples
keep <- complete.cases(meta[,c("Subject.ID", "Age.months", "sex.numeric", "RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "RNA.isolation.Batch", "Lib_prep_batches", "Year.Drawn")])
meta <- meta[keep,]
meta$Subject.ID <- factor(meta$Subject.ID)

expr <- cpm(dge, log=TRUE)[1:50, rownames(meta)]

form1 <- ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + Age.months

meta$Age.months <- scale(meta$Age.months)[,1]
form2 <- ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + Age.months

vp1 <- fitExtractVarPartModel(expr, form1, meta, showWarnings=F)
vp2 <- fitExtractVarPartModel(expr, form2, meta, showWarnings=F)

print("UNSCALED MEDIAN SUBJECT VARIANCE:")
print(median(vp1$Subject.ID, na.rm=T))
print("SCALED MEDIAN SUBJECT VARIANCE:")
print(median(vp2$Subject.ID, na.rm=T))
