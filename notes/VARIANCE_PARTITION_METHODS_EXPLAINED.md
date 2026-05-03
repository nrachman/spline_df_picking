# Variance Partitioning Methods: Standard Parameters vs. Realized Predictions

When fitting mixed-effects models to RNA-seq data, the goal of variance partitioning is to determine what percentage of the total variation in gene expression is driven by each covariate (e.g., Age, Subject, Batch). 

There are two primary mathematical approaches to calculate these percentages. While they often yield similar results in simple models, understanding the conceptual difference between them is critical—especially when modeling complex, non-linear trajectories like splines.

---

## 1. The Standard Approach: Variance Parameters ($\sigma^2$)
*(This is the default method used internally by `variancePartition::calcVarPart` for random effects)*

In a linear mixed model (LMM), random effects are assumed to be drawn from a normal distribution. For example, the effect of `Subject.ID` is modeled as:
$u_{subject} \sim N(0, \sigma^2_{subject})$

### How it works:
1.  **Random Effects:** The algorithm directly extracts the estimated variance parameter ($\sigma^2_{subject}$, $\sigma^2_{batch}$) from the fitted model. 
2.  **Residuals:** It extracts the residual variance ($\sigma^2_{residual}$).
3.  **Fixed Effects:** Because fixed effects (like `Sex` or linear `Age`) do not have a "variance parameter" in the same way, `variancePartition` calculates the variance of the linear predictor for the specific sample: $Var(X \hat{\beta})$.
4.  **Total Variance:** It sums these components: $V_{total} = \sigma^2_{subject} + \sigma^2_{batch} + Var(X \hat{\beta}) + \sigma^2_{residual}$
5.  **Fractions:** It divides each component by $V_{total}$.

### Pros & Cons:
*   **Pro:** Extremely fast computationally because it mostly just reads pre-calculated parameters from the model summary.
*   **Con (The Mismatch):** It mixes two different philosophical scales. It compares the *theoretical population variance* of the random effects ($\sigma^2$) against the *empirical sample variance* of the fixed effects ($Var(X \hat{\beta})$). 
*   **Con (The Overestimation due to Shrinkage):** Mixed models calculate individual subject baselines (Best Linear Unbiased Predictors, or BLUPs) by "shrinking" extreme values toward the global mean, especially for noisy or sparse data. Because the individual estimates are mathematically shrunk toward the center, the actual sample variance of your subjects' BLUPs is mathematically forced to be *smaller* than the estimated population variance ($\sigma^2_{subject}$). Therefore, the standard method ($\sigma^2$) will systematically allocate a **larger** percentage of variance to `Subject.ID` than what is actually realized in the physical data points of your experiment.

---

## 2. The Empirical Approach: Realized Predictions ($Var(\hat{y})$)
*(This is the method we must use to properly evaluate complex splines and GAMMs)*

Instead of looking at the abstract population parameters, this approach looks at the actual predictions the model makes for the specific data points in your experiment.

### How it works:
It calculates the specific fitted value ($\hat{y}$) contributed by each component for every single observation, and then takes the sample variance of those predictions.
1.  **Fixed Effects (e.g., Age):** Generate the prediction based only on Age: $\hat{y}_{age} = X_{age} \hat{\beta}_{age}$. Then calculate the sample variance: $Var(\hat{y}_{age})$.
2.  **Random Effects (e.g., Subject):** Extract the actual Best Linear Unbiased Predictors (BLUPs) for each subject in your dataset, and multiply them by the design matrix: $\hat{y}_{subj} = Z_{subj} \hat{u}_{subj}$. Then calculate the sample variance: $Var(\hat{y}_{subj})$.
3.  **Residuals:** Calculate the variance of the actual model residuals: $Var(\epsilon)$.
4.  **Fractions:** Sum these sample variances and calculate the percentages.

### Pros & Cons:
*   **Pro (Scale Consistency):** It puts fixed effects, random effects, and residuals on the exact same scale: the empirical variance of the actual sample data. 
*   **Pro (Spline Compatibility):** This is the **only** mathematically valid way to partition variance for Generalized Additive Mixed Models (GAMMs). In `gamm4`, the $\sigma^2$ parameter for a spline represents a *smoothing penalty*, not biological variance. You *must* use the realized predictions to see how much variance the curve actually explains.
*   **Con:** Computationally heavier because you have to perform matrix multiplications for every observation across every model component.

---

## 3. Why This Matters for Your Study

In standard models without splines, `variancePartition`'s hybrid approach works wonderfully. 

However, you are dealing with a complex overlapping cohort design (ages 0-14) and considering multiple spline bases (`ns` vs `lspline`) to model non-linear development. 

1. **Spline Aggregation:** When you use `ns(Age, df=3)`, it creates three separate columns in the fixed effects matrix. To find the "Total Variance Explained by Age", you cannot just add their individual variances together (because the spline components covary). You must use the Realized Predictions approach: calculate the combined prediction for the whole spline ($\hat{y}_{age} = X_{ns1}\beta_1 + X_{ns2}\beta_2 + X_{ns3}\beta_3$) and take the variance of that total prediction. *(Note: modern versions of `variancePartition::calcVarPart` actually do this aggregation for fixed-effect splines correctly under the hood!)*
2. **The Subject vs. Age Debate:** Your reviewers are skeptical about the high variance attributed to `Subject.ID`. They may (rightfully) suspect that the standard $\sigma^2$ approach is mathematically inflating the importance of Subject compared to the fixed Age effect due to the difference between population parameter estimates and the actual "shrunken" BLUPs in the sample.
3. **The Solution:** By evaluating the **Realized Predictions**, you calculate the variance of the *shrunken, conservative* BLUPs that `Subject.ID` actually contributed to your specific 690 samples. If the Realized Prediction variance for Subject remains massively higher than the Realized Prediction variance for Age (even when using AIC to give Age the most flexible `df` possible), you have an irrefutable, empirical mathematical argument for the reviewers. You are proving the Subject effect is massive even when using the most conservative metric available.

# Variance Partition Discrepancy Analysis

## Overview
Recent benchmarking of cell-frequency variance (CBC vs. Deconvolution) showed lower variance explained by cell types compared to previous results. Investigation revealed that this is primarily due to **gene selection bias** and the inclusion of **Subject ID** as a random effect.

## Key Findings

### 1. HVGs vs. Random Genes
When analyzing the top 100 Highly Variable Genes (HVGs), the variance is dominated by Subject Identity (individuality). When looking at 100 Random Genes (representing the broader transcriptome), the variance explained by Cell Frequencies quadruples.

| Gene Set | Subject Var (Mean) | CellFreq Var (Mean) | Residual Var (Mean) |
| :--- | :--- | :--- | :--- |
| **Top 100 HVGs** | 54.0% | 4.9% | 21.6% |
| **Random 100 Genes** | 10.3% | 20.7% | 60.9% |

**Conclusion:** HVGs are genes that are "highly variable" specifically because they vary between individuals. This biological signal "crowds out" the compositional signal from cell frequencies. Random genes have less individuality, making the cell composition signal relatively larger.

### 2. The Role of Subject ID
Including `Subject.ID` as a random effect captures intrinsic individual differences. If `Subject.ID` is removed from the model for random genes, the CellFreq variance remains stable (~21%), but the Residual variance increases to compensate for the lost Subject signal.

### 3. Realized vs. Computed Variance Explained
Standard variance partitioning relies on theoretical population parameters ($\sigma^2$), which can be sensitive to singular fits or inflated by outliers in small cohorts.
- **Computed (Standard):** Uses the estimated variance components directly from the model object (`VarCorr`).
- **Realized:** Calculates the empirical variance of the Best Linear Unbiased Predictors (BLUPs) for these specific subjects.
- **Result:** The "Realized" method is more conservative and provides a more grounded, empirical scale for comparing the magnitude of Subject vs. Age effects, preventing the overestimation that can occur with theoretical population averages.

### 4. Spline Complexity
Using `lspline(n=4)` for Age partitions developmental variation more granularly than a linear Age model. However, this has a smaller impact on CellFreq variance than the choice of gene set.

## Recommendation for Reporting
- Acknowledge that variance components are relative to the gene set analyzed.
- For global transcriptome summaries, random gene subsets or all-gene averages are more appropriate for showing cell-frequency impacts.
- For high-resolution trajectory modeling (the focus of this benchmark), focusing on HVGs is correct as it captures the most robust biological signals (Subject and Age).
