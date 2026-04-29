# Note on Variable Scaling in Mixed-Effects Models

During the evaluation of the spline models in `run_spline_benchmark.R`, an important issue with numerical stability was identified when comparing Natural Cubic Splines (`ns`) with Linear Splines (`lspline`).

## The Problem: Subject Variance Collapse
Initially, the benchmark showed that the variance explained by `Subject.ID` was highly unstable and often collapsed to near zero (e.g., $10^{-12}$) when using `lspline` at lower degrees of freedom, whereas the `ns` models maintained a robust and expected distribution of Subject variance.

## The Cause: Numerical Instability
The underlying cause is a known issue in Mixed-Effects Modeling, specifically within the `lme4` optimizer used by `variancePartition`:

1. **Differing Scales:** The continuous covariates in the model (e.g., `Age.months` ranging from ~8 to 176, or cell percentages like `Lymph.pcnt` ranging from 0 to 100) exist on vastly different scales compared to the small, subtle technical batch effects.
2. **Natural Splines (`ns`):** The `splines::ns()` function automatically standardizes the input variable and projects it onto an orthogonal basis matrix (roughly bounded between -1 and 1) before the model is fit. This inadvertently protects the optimizer from numerical scaling issues.
3. **Linear Splines (`lspline`):** In contrast, `lspline::elspline()` passes the raw, unscaled `Age.months` directly into the formula as a piecewise linear term. 
4. **Optimizer Failure:** When faced with these large, unscaled variances, the penalized solver in `lme4` struggles to converge on the variance components. To mathematically "solve" the equations, it often hits a boundary constraint, effectively dropping the `Subject.ID` variance parameter down to 0.

## The Solution: Standardizing Continuous Covariates
To ensure a fair and robust comparison between the spline families, it is necessary to **standardize all continuous covariates** (mean = 0, standard deviation = 1) before they enter the model.

In `run_spline_benchmark.R`, this is achieved by applying the `scale()` function to `Age.months`, `RIN`, and all cell percentage fractions. 

```R
continuous_vars <- c("RIN", "mk_dup.PERCENT_DUPLICATION", "star.uniquely_mapped_percent", 
                     "Age.months", "Lymph.pcnt", "Mid.pcnt", "Gran.pcnt", "MON.pcnt", "EOS.pcnt")
for (v in continuous_vars) meta[[v]] <- scale(meta[[v]])[,1]
```

### Interpretation
Standardizing predictors does **not** change the biological interpretation or the underlying shape of the models:
- The spline knots still fall at the exact same relative percentiles of the data.
- The variance partitioning fractions are mathematically equivalent to what they would be with perfect, infinite-precision raw data.
- It simply stabilizes the underlying math, preventing boundary collapse and ensuring accurate variance estimates for `Subject.ID`.