# Bootstrap Standard Errors for Covariate Model Parameters

Estimates standard errors and p-values for the covariate structural
model parameters using nonparametric bootstrap resampling. In each
replicate, a new model is fitted on a bootstrap sample and class labels
are aligned to those of the original model via globally optimal
assignment (solving the linear sum assignment problem by enumeration for
K \<= 8), correcting for label switching. Results can be passed directly
to
[`confint.mixture_model`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md)
and
[`wald_omnibus_test`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md).

## Usage

``` r
bootstrap_covariates(
  model_state,
  X,
  Y,
  n_reps = 100,
  random_state = 123,
  ref_class = 1
)
```

## Arguments

- model_state:

  A fitted `mixture_model` object with a covariate structural model
  (fitted with `structural = "covariate"`).

- X:

  Numeric matrix. The original measurement data used to fit
  `model_state`.

- Y:

  Numeric matrix. The original covariate data used to fit `model_state`.

- n_reps:

  Positive integer. Number of bootstrap replications. Default is `100`.

- random_state:

  Integer seed for reproducibility. Default is `123`.

- ref_class:

  Integer. Reference class for centering bootstrap betas. Should match
  the `ref_class` used in subsequent calls to
  [`confint.mixture_model`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md)
  and
  [`wald_omnibus_test`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md).
  Default is `1`.

## Value

A named list with four elements:

- `standard_errors` K x D numeric matrix of bootstrap standard errors
  (one row per class, one column per predictor).

- `p_values` K x D numeric matrix of two-sided p-values.

- `boot_betas` 3-D array of shape `(n_reps, K, D)` containing aligned,
  reference-centered bootstrap beta estimates. Required by
  [`wald_omnibus_test()`](https://pdvalencia.github.io/mixtureEM/reference/wald_omnibus_test.md).

- `orig_betas` K x D matrix of the original model's beta coefficients,
  centered on `ref_class`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
Z <- matrix(rnorm(100), nrow = 100)
colnames(Z) <- "age"
fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
                   structural = "covariate",
                   n_steps = 3, correction = "ML", n_init = 5)
boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
boot$p_values
confint(fit, boot_results = boot)
} # }
```
