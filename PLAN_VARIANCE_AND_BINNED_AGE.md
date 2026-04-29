# Plan: Advanced Variance Partitioning and Binned Age Modeling

This plan outlines the next steps for refining the variance partitioning analysis and exploring alternative ways to model Age as a random effect.

## 1. Comparison of Variance Partitioning Methods
The goal is to determine if the "Realized Predictions" approach provides a more conservative and robust estimate of Subject variance compared to the standard population-parameter approach.

- **Task:** Implement the "Realized Predictions" variance partitioning method as described in `notes/VARIANCE_PARTITION_METHODS_EXPLAINED.md`.
- **Scope:** 1,000 highly variable genes (HVGs).
- **Model Configuration:** Use 4 degrees of freedom ($df=4$) for both Natural Splines (`ns`) and Linear Splines (`lspline`).
- **Optimization:** Use `bobyqa` for all fits.
- **Output:** A comparison table and plot showing the % variance explained by Subject, Age, and Residuals for both the "Standard" and "Realized" methods.

## 2. Modeling Age as a Random Effect (Binned Approach)
To address reviewer concerns regarding the scale of Subject vs. Age effects, we will explore treating Age as a categorical random effect.

- **Binning Strategy:**
    - Create a new variable `Age_Binned` by dividing the samples into 6 equally sized bins based on their `Age.months` rank (quantiles).
- **Model Formula:**
    - `~ (1|Subject.ID) + (1|Age_Binned) + (1|RNA.isolation.Batch) + ... [rest of covariates]`
- **Comparison:** Compare the variance partition of `(1|Age_Binned)` vs. `(1|Subject.ID)`.

## 3. Age in Years as Categorical
Compare the spline approach to a simpler categorical approach.

- **Transformation:** Create `Age_Years_Cat = factor(floor(Age.months / 12))`.
- **Scenarios:**
    1. **Fixed Effect:** Model Age as a categorical fixed effect.
    2. **Random Effect:** Model Age as `(1|Age_Years_Cat)`.
- **Benchmarking:** Compare AIC and Variance Fractions between these categorical models and the $df=4$ spline models.

## 4. Implementation Details
- **Script:** Create `scripts/run_advanced_variance_benchmarks.R`.
- **Optimizer:** Enforce `lmerControl(optimizer = "bobyqa")` across all models.
- **Standardization:** Maintain the scaling of all continuous covariates (RIN, Cell Percentages, etc.).

---

## Clarifying Questions for the User
1. **Binning Logic:** For the "equally sized bins of 6 based on their rankings," did you mean 6 bins total (e.g., ~115 samples per bin) or bins that contain 6 samples each (e.g., ~115 bins total)?
Response: 115 bins total, with each bin containing 6 samples (approximately because some samples might be dropped due to incomplete data)

2. **Age Categorical Floor:** For Age in Years as categorical, is `floor(Age.months / 12)` the correct grouping, or should I use a different interval?
Response: can you use `round(Age.months / 12)` instead but use the correct number of digits for the rounding

3. **Spline DF for Comparison:** For the categorical comparison, is $df=4$ the correct spline baseline, or should we use the "winner" from the AIC tournament?

Response: I like the AIC tournament but I think it suggest that we should just use 1 df for age, but there are some non-linear possibilities that I would like to capture so here I am erring on the side of slightly overfitting the age trend for a few genes to capture the correct one for others

4. **GitHub Remote:** Please provide the URL for the GitHub repository so I can push the organized codebase.

Response: https://github.com/nrachman/spline_df_picking