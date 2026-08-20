# Print a Brief Summary of a Fitted Mixture Model

Prints a compact overview of the fitted model including: number of
classes, estimation method, convergence status, the fit indices, and
estimated class proportions. The indices shown are the same ones
[`compare_mixtures`](https://pdvalencia.github.io/mixtureEM/reference/compare_mixtures.md)
tabulates, so one printed model and a table over a range of K can be
read together; lower AIC, BIC and SABIC are better. For full parameter
tables, use
[`summary.mixture_model`](https://pdvalencia.github.io/mixtureEM/reference/summary.mixture_model.md)
(structural parameters) or
[`measurement_summary`](https://pdvalencia.github.io/mixtureEM/reference/measurement_summary.md)
(item parameters).

## Usage

``` r
# S3 method for class 'mixture_model'
print(x, ...)
```

## Arguments

- x:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

Invisibly returns `x`. Called for its printed side-effect.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(300, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
print(fit)
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 136 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -203.34
#>   Parameters     : 7
#>   AIC            : 420.69
#>   BIC            : 438.92
#>   SABIC          : 416.82
#>   Rel. Entropy   : 0.2215
#>   Best solution  : found by 20 of 20 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 57.78%
#>   Class 2: 42.22%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
# or equivalently:
fit
#> =========================================================
#>                   LATENT MIXTURE MODEL                   
#> =========================================================
#> Classes Estimated  : 2
#> Estimation Method  : 1-step
#> Converged          : TRUE (in 136 iterations)
#> ---------------------------------------------------------
#>   Log-Likelihood : -203.34
#>   Parameters     : 7
#>   AIC            : 420.69
#>   BIC            : 438.92
#>   SABIC          : 416.82
#>   Rel. Entropy   : 0.2215
#>   Best solution  : found by 20 of 20 starts
#> ---------------------------------------------------------
#> Class Weights (Sizes):
#>   Class 1: 57.78%
#>   Class 2: 42.22%
#> =========================================================
#> Type summary(model) for structural parameters or measurement_summary(model) for item parameters.
```
