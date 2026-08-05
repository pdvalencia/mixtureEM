# Summarise a Fitted Mixture Model

Prints a detailed summary of the structural model parameters. Depending
on which structural model was fitted, this includes covariate regression
coefficients (as odds ratios with 95\\ distal outcome means, or
class-specific regression effects. If no structural model is present, a
notice is printed directing the user to
[`measurement_summary`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md).

## Usage

``` r
# S3 method for class 'mixture_model'
summary(object, ref_class = NULL, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ref_class:

  Integer. The reference latent class for pairwise contrasts. Defaults
  to the first class (`1`).

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

Invisibly, a list holding the printed numbers in tidy form, ready for
further use: `$coefficients` (one row per class contrast and covariate:
estimate, SE, p, odds ratio and its confidence limits), `$omnibus` (the
per-covariate omnibus Wald tests), and, when a distal outcome is
present, `$outcome` (predicted probabilities or class means/estimates
with their tests). Returns `NULL` when the model has no structural part.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
summary(fit)
#> Notice: No structural model found. Use measurement_summary() for item parameters.
```
