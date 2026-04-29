# Plan: Rigorous Spline Selection and Benchmarking

## 1. Objective
The goal of this plan is to establish a rigorous, empirical framework for modeling non-linear age trajectories (ages 0–14) in longitudinal RNA-seq data. 

Specifically, we aim to:
1.  **Benchmark Spline Types:** Compare Natural Cubic Splines (`ns`) against Piecewise Linear Splines (`lspline`) using the full model architecture defined in the established `run_varpart.R` pipeline.
2.  **Determine Optimal Complexity:** Rigorously select the optimal degrees of freedom (df) or number of knots ($n \in [1, 10]$) using the **Global Empirical Prior** approach.
3.  **Assess Robustness:** Evaluate how spline complexity impacts the stability of Subject variance (individuality) and the effect sizes of other critical covariates, including technical factors and cell frequencies.
4.  **Ensure Compatibility:** Guarantee that the winning approach integrates seamlessly with the `variancePartition` and `dream` ecosystem for both variance partitioning and hypothesis testing.

---

## 2. Methodology: The Global Empirical Prior Benchmark

Evaluating thousands of models across 15,000 genes is computationally wasteful and statistically noisy. Instead, we will establish a Global Empirical Prior using a representative subset, executed from the project root: `.../Nicaragua_immune_intrinsicness_bulk/`.

### Step 1: Data Integration & Cleaning
We will load and merge the expression metadata with Complete Blood Count (CBC) data to control for cell type proportions.
```R
# Load CBC data
cbc <- read_tsv("data/inputs/cbc/cbc_cleaned.tsv")
cbc <- cbc %>% mutate(Sample_ID = gsub("\\.", "_", Sample_ID))

# Filter for case samples and valid RIN
keep_samples <- meta$Sample.name[meta$Sample.type == "case" & !is.na(meta$RIN)]
keep_samples <- intersect(keep_samples, cbc$Sample_ID)

# Prepare cell frequency matrix
keep_cells <- c("Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt", "BAS.pcnt")
metadata <- meta %>% 
  filter(Sample.name %in% keep_samples) %>%
  left_join(cbc %>% select(Sample_ID, all_of(keep_cells)), by = c("Sample.name" = "Sample_ID"))
```

### Step 2: The Modeling Tournament
For each of 500 highly variable genes (HVGs), we will fit models using the `variancePartition::dream` framework. The baseline formula follows `run_varpart.R`, augmented with cell frequencies:

`~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + [Age Spline] + sex.numeric + RIN + mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + BAS.pcnt`

We will test two families of splines across 10 levels of complexity:
*   **Family A: Natural Cubic Splines (`splines::ns`)**
    *   Models: `ns(Age.months, df=1)` through `ns(Age.months, df=10)`.
*   **Family B: Piecewise Linear Splines (`lspline::elspline`)**
    *   Models: `elspline(Age.months, n=1)` through `elspline(Age.months, n=10)`.

---

## 3. Evaluation Metrics

### Metric A: Statistical Fit and Parsimony (AIC)
AIC remains our primary metric to ensure we capture developmental dynamics without underfitting, addressing reviewer skepticism regarding the dominance of Subject identity.

### Metric B: Robustness of Subject Identity (Variance and Effect Sizes)
We will monitor "variance stealing" by tracking the Percentage of Variance Explained by `Subject.ID` and the stability of Subject-level BLUPs as spline complexity increases.

### Metric C: Robustness of Inference (Age P-values and Covariate Stability)
We will perform joint F-tests (ANOVA) on the spline terms to extract p-values for the Age effect. Crucially, we will verify that increasing spline complexity does not destabilize the coefficients or p-values for stable covariates (e.g., `sex.numeric` or cell proportions).

---

## 4. Alternative Methods Considered

`ns` remains the strongest candidate because it offers local, smooth biological realism and computational efficiency compatible with the `variancePartition` ecosystem, whereas alternatives like Penalized Splines (`gamm4`) are computationally prohibitive for this scale.

---

## 5. Implementation & Deployment

1.  **Run Benchmark:** Execute the tournament script using the comprehensive model formula.
2.  **Declare Winner:** Identify the spline type and `df` that minimizes AIC while preserving Subject variance stability.
3.  **Final Fit:**
    ```R
    # Example: Natural splines with df=3
    form <- ~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn) + 
              ns(Age.months, df=3) + sex.numeric + RIN + 
              mk_dup.PERCENT_DUPLICATION + star.uniquely_mapped_percent + 
              Lymph.pcnt + Mid.pcnt + Gran.pcnt + MON.pcnt + EOS.pcnt + BAS.pcnt
    
    vobj <- voomWithDreamWeights(counts, form, metadata)
    fit  <- dream(vobj, form, metadata)
    ```
4.  **Extract Results:** Use `extractVarPart(fit)` for variance decomposition and a custom joint F-test wrapper for high-resolution hypothesis testing of the Age effect.
