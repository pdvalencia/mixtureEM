# Confidence Intervals for Odds Ratios in a Mixture Model

Computes confidence intervals for the odds ratios of covariate effects
on latent class membership, using whichever analytical variance the
model was fitted with — see the `se` argument of
[`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
and
[`covariate_se`](https://pdvalencia.github.io/mixtureEM/reference/covariate_se.md).
The estimator used is named in the printed output and carried in the
`method` attribute of the result. Bootstrapped standard errors from
[`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md)
can be supplied instead, which is worth doing in small samples, when
classes overlap heavily, or as a check on the analytical intervals.

## Usage

``` r
# S3 method for class 'mixture_model'
confint(
  object,
  parm = NULL,
  level = 0.95,
  boot_results = NULL,
  ref_class = 1,
  ...
)
```

## Arguments

- object:

  A fitted `mixture_model` object with a covariate structural model
  (fitted with `structural = "covariate"`).

- parm:

  A specification of which parameters are to be given confidence
  intervals. Currently ignored (intervals are returned for all covariate
  parameters).

- level:

  Numeric. The confidence level. Default is `0.95`.

- boot_results:

  Optional. A list returned by
  [`bootstrap_covariates`](https://pdvalencia.github.io/mixtureEM/reference/bootstrap_covariates.md).
  When provided, bootstrapped standard errors are used. When `NULL`
  (default), the analytical variance stored on the fitted model is used.

- ref_class:

  Integer. The reference latent class. Odds ratios for all other classes
  are expressed relative to this class. Default is `1`.

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

A named list with one data frame per predictor variable (including the
intercept). Each data frame has columns `OR` (odds ratio), `Lower`, and
`Upper` (confidence bounds), with one row per latent class. Values are
rounded to three decimal places.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
Z <- matrix(rnorm(100), nrow = 100)
colnames(Z) <- "age"
fit <- fit_mixture(X, Y = Z, n_components = 2, measurement = "binary",
                   structural = "covariate",
                   n_steps = 3, correction = "ML", n_init = 5)
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.

# Analytical CIs (from Hessian)
confint(fit)
#> =========================================================
#>         CONFIDENCE INTERVALS FOR ODDS RATIOS             
#> =========================================================
#> Reference Class : 1
#> Level           : 95%   Method: Bakk-Oberski-Vermunt corrected (robust step 3, hessian step 1)
#> ---------------------------------------------------------
#>                          OR    Lower    Upper
#> 
#> Intercept
#>   Class 1 (Ref)       1.000        -        -
#>   Class 2             0.696    0.004  119.055
#> 
#> age
#>   Class 1 (Ref)       1.000        -        -
#>   Class 2             0.266    0.028    2.550
#> =========================================================

if (FALSE) { # \dontrun{
# Bootstrapped CIs
boot <- bootstrap_covariates(fit, X, Z, n_reps = 200)
confint(fit, boot_results = boot)
} # }
```
