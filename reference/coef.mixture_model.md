# Extract Covariate Odds Ratios from a Fitted Mixture Model

Extracts the logistic regression coefficients from a covariate
structural model and returns them as a matrix of odds ratios, centered
on a reference class. Only available when the model was fitted with
`structural = "covariate"`.

## Usage

``` r
# S3 method for class 'mixture_model'
coef(object, ref_class = 1, covariate_names = NULL, ...)
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

- ...:

  Currently unused. Present for S3 method compatibility.

## Value

A K x D numeric matrix of odds ratios, where rows are latent classes and
columns are predictors (including the intercept). The reference class
row will always show `1.000`.

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
#>               Intercept       age
#> Class 1 (Ref) 1.0000000 1.0000000
#> Class 2       0.6961114 0.2660264
```
