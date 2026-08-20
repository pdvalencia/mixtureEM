# Extract Covariate Odds Ratios from a Fitted Mixture Model

Extracts the logistic regression coefficients from a covariate
structural model and returns them as a matrix of odds ratios, centered
on a reference class. Only available when the model was fitted with
`structural = "covariate"`.

## Usage

``` r
# S3 method for class 'mixture_model'
coef(object, ref_class = 1, covariate_names = NULL, exponentiate = TRUE, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object with a covariate structural model.

- ref_class:

  Integer. The reference class for centering. All other class odds
  ratios are expressed relative to this class. Default is `1`.

- covariate_names:

  Optional character vector of predictor names to override the column
  names stored in the model. Default is `NULL`.

- exponentiate:

  Logical. When `TRUE` (the default) the coefficients are returned as
  odds ratios; when `FALSE`, as the multinomial-logit coefficients
  themselves, relative to `ref_class`.

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

A K x D numeric matrix, where rows are latent classes and columns are
predictors (including the intercept). With `exponentiate = TRUE` these
are odds ratios and the reference class row is all `1`; with
`exponentiate = FALSE` they are log-odds coefficients and that row is
all `0`.

## Details

The printed summary reports odds ratios because that is the scale these
effects are interpreted and published on, and the default here matches
it. `coef(fit, exponentiate = FALSE)` and
[`vcov`](https://pdvalencia.github.io/mixtureEM/reference/vcov.mixture_model.md)
give the log-scale estimates and their standard errors, for anyone who
needs to compare them against another program or pool them across
analyses. Both scales are exact — `log(coef(fit))` has always recovered
the coefficients, since the odds ratios are returned at full double
precision; the argument makes that discoverable rather than a trick.

## See also

[`confint`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md)
for intervals on the odds-ratio scale, and
[`vcov`](https://pdvalencia.github.io/mixtureEM/reference/vcov.mixture_model.md)
for the log-scale standard errors.

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
coef(fit)
#>               Intercept      age
#> Class 1 (Ref) 1.0000000 1.000000
#> Class 2       0.7390351 0.339548
coef(fit, exponentiate = FALSE)
#>                Intercept      age
#> Class 1 (Ref)  0.0000000  0.00000
#> Class 2       -0.3024098 -1.08014
```
