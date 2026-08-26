# Covariance Matrix of the Class-Membership Coefficients

Returns the variance-covariance matrix of the multinomial-logit
coefficients for covariate effects on latent class membership, so that
`sqrt(diag(vcov(fit)))` gives their standard errors. This is the
log-scale companion to
[`confint`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md),
which reports intervals on the odds-ratio scale — it exists so that the
standard errors can be got at directly, for comparing against another
program or for pooling estimates, rather than reconstructed from an
interval width.

Which estimator it comes from depends on how the model was fitted: the
survey-robust or step-3 corrected covariance when one was computed, and
the Q-function Hessian otherwise. The estimator is named in the `method`
attribute, exactly as
[`confint()`](https://rdrr.io/r/stats/confint.html) reports it.

## Usage

``` r
# S3 method for class 'mixture_model'
vcov(object, ...)
```

## Arguments

- object:

  A fitted `mixture_model` object with a covariate structural model
  (fitted with `structural = "covariate"`).

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

A square numeric matrix over the free coefficients, of dimension
`(K - 1) * D` where `K` is the number of classes and `D` the number of
predictors including the intercept. Rows and columns are named
`"Class k:predictor"`. The estimator's name is carried in the `method`
attribute and the reference class in `ref_class`.

## Which contrast these are variances of

One class is pinned at zero to identify the model, and the free
coefficients are the other `K - 1` classes relative to **that** class.
Which class it is depends on the fit, because the classes are reordered
by size after estimation, so it is reported in the `ref_class` attribute
rather than fixed. This is not necessarily the same contrast as
[`coef()`](https://rdrr.io/r/stats/coef.html)'s and
[`confint()`](https://rdrr.io/r/stats/confint.html)'s `ref_class`
argument, which defaults to class 1: pair these standard errors with
`coef(fit, exponentiate = FALSE, ref_class = attr(vcov(fit), "ref_class"))`,
whose non-reference rows are then exactly the coefficients they belong
to.

## See also

[`coef`](https://pdvalencia.github.io/mixtureEM/reference/coef.mixture_model.md),
whose `exponentiate = FALSE` gives the coefficients these are the
variances of, and
[`confint`](https://pdvalencia.github.io/mixtureEM/reference/confint.mixture_model.md).

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
sqrt(diag(vcov(fit)))
#> Class 2:Intercept       Class 2:age 
#>         2.3584007         0.6881262 
```
