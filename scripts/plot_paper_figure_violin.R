library(tidyverse)
library(variancePartition)

# 1. Load Original Paper Full Transcriptome Data
VP_NICA_IN_PATH <- "../../data/analysis_out/variancePartition/varPart_include_cbc_fardeep_lm22_summed.rds"
vp <- readRDS(VP_NICA_IN_PATH)
vp <- data.frame(vp)
cell_cols <- startsWith(colnames(vp), "cellfreq")
cell_vp <- rowSums(vp[, cell_cols])
vp_sub <- vp[, !cell_cols]
vp_sub <- vp_sub %>%
        mutate(CellFreqs = cell_vp)

# Ensure order of columns for plotting (often variancePartition sorts by median by default)
# plotVarPart expects a data.frame where columns are the variance components
p_full <- plotVarPart(vp_sub) + 
        theme_classic() +
        theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in")) +
        labs(title = "Full Transcriptome (Paper Figure)",
             subtitle = "Cell frequencies explain a huge proportion of variance across 16,599 genes")
             
ggsave(plot = p_full, filename = "output/images/paper_violin_full_transcriptome.png", height = 6, width = 7)

# 2. Let's also make plotVarPart violins for our Discrepancy Data for comparison
df_discrep <- readRDS("output/rds/vp_discrepancy_data.rds")
# Split by GeneSet
df_hvg <- df_discrep %>% filter(GeneSet == "Top 250 HVGs") %>%
  select(Subject.ID, Age.months, Summed_CellFreq, Residuals) %>%
  rename(Subject = Subject.ID, Age = Age.months, CellFreqs = Summed_CellFreq) %>%
  as.data.frame()
  
p_hvg <- plotVarPart(df_hvg) +
        theme_classic() +
        theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in")) +
        labs(title = "Top 250 HVGs",
             subtitle = "Subject identity dominates; cell frequencies are minimized")
ggsave(plot = p_hvg, filename = "output/images/paper_violin_hvg.png", height = 6, width = 7)

df_rand <- df_discrep %>% filter(GeneSet == "Random 250 Genes") %>%
  select(Subject.ID, Age.months, Summed_CellFreq, Residuals) %>%
  rename(Subject = Subject.ID, Age = Age.months, CellFreqs = Summed_CellFreq) %>%
  as.data.frame()

p_rand <- plotVarPart(df_rand) +
        theme_classic() +
        theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in")) +
        labs(title = "Random 250 Genes",
             subtitle = "Mirrors the full transcriptome: High cell frequency variance")
ggsave(plot = p_rand, filename = "output/images/paper_violin_rand.png", height = 6, width = 7)

cat("Violin plots generated!\n")
