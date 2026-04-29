# Note on Boundary (Singular) Fits in the Spline Benchmark

During the execution of the spline benchmark (`run_spline_benchmark.R`), you will likely see a large number of warnings output to the console or log file:

`boundary (singular) fit: see help('isSingular')`

## What is a Singular Fit?
A "singular fit" occurs in Mixed-Effects Models (like those used in `lmer`, `variancePartition`, and `dream`) when the model estimates the variance of one (or more) of your random effects to be exactly **zero**. 

In your comprehensive benchmark formula, you have four random effects:
`~ (1|Subject.ID) + (1|RNA.isolation.Batch) + (1|Lib_prep_batches) + (1|Year.Drawn)`

For many individual genes, one or more of these technical batches (e.g., `Year.Drawn` or `RNA.isolation.Batch`) might not actually contribute any measurable variance above the general background noise. When the model tries to estimate that tiny variance, it hits the lower mathematical boundary (0) and flags a "singular fit".

## Is this an issue for the benchmark?

**No.** It is a common occurrence in high-throughput transcriptomics and is generally not a problem for the goals of this specific benchmark:

1. **Variance Partitioning:** `variancePartition` is designed to handle this gracefully. When a random effect has zero variance, it simply assigns `0%` variance explained to that covariate in the final output. The key metric—`Subject.ID` variance—will still be accurately estimated.
2. **P-values (Age Inference):** The `dream` framework is highly robust. Even if a nuisance random effect (like a sequencing batch) has zero variance, the fixed effects (like your Age splines) and their standard errors are still estimated correctly, yielding valid p-values.
3. **AIC/BIC Comparisons:** A singular fit means the model is technically overparameterized for that *specific gene*. However, the AIC is still mathematically valid and usable for comparing spline complexities (since the random effect structure is held constant across your comparisons).

## Should I try to fix it?

**No action is strictly required.** Because you are running a "Global Empirical Prior" across 1,000 highly variable genes, you *want* to keep the model formula constant so the AIC and variance results are perfectly comparable across all genes and all spline complexities. Trying to dynamically drop random effects gene-by-gene would make the benchmark much harder to interpret and execute at scale.

If you see these warnings, it just means: *"For this specific gene, one of the batch effects didn't contribute to the variance."*