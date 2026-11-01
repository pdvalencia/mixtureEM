# Fit a model with two associated latent class variables

Fits two latent class variables at once, each measured by its own items:
`x` with `n_x` classes and `y` with `n_y` classes, linked by a
multinomial logistic regression of `y` on `x` (Hagenaars, 1993). The
association is reported against the last class of each variable: the
intercept of `y` class j is its log-odds against the last `y` class
among members of the last `x` class, and the coefficient of `x` class k
is how much membership in `x` class k changes that log-odds. The joint
class probabilities and the probabilities of each `y` class given each
`x` class come with it.

## Usage

``` r
fit_joint_classes(
  x_items,
  y_items,
  n_x = 2,
  n_y = 2,
  measurement_x = c("binary", "ordinal"),
  measurement_y = c("binary", "ordinal"),
  predictors_x = NULL,
  predictors_y = NULL,
  y_effects = "common",
  association = c("free", "none"),
  predictors_items = NULL,
  n_steps = 1,
  step1 = NULL,
  smoothing = 1,
  bayes_constants = NULL,
  n_init = 20,
  random_state = NULL,
  n_cores = .default_n_cores(),
  standard_errors = TRUE,
  ...
)
```

## Arguments

- x_items, y_items:

  The items measuring each latent variable, one row per case, the same
  cases in the same order. Missing values are allowed.

- n_x, n_y:

  Integer. Number of classes of each variable.

- measurement_x, measurement_y:

  `"binary"` (items coded 0/1) or `"ordinal"` (ordered categories coded
  from 0 or from 1). A binary and an ordinal side are fitted together as
  ordinal items, which for a two-category item is the same model.

- predictors_x:

  A data frame of covariates predicting `x` class membership, or `NULL`.

- predictors_y:

  A data frame of covariates predicting `y` class membership, or `NULL`.

- y_effects:

  How each covariate in `predictors_y` acts: `"common"` (one slope per
  `y` class, the default), `"by_x"` (a slope per `x` class for every
  covariate), or a named vector such as `c(female = "by_x")`, the
  unnamed covariates acting in common.

- association:

  `"free"` (the default) estimates the association; `"none"` holds the
  two variables independent given the covariates, so `y_effects` cannot
  name `"by_x"`.

- predictors_items:

  Direct effects of covariates on items (uniform differential item
  functioning): a named list, item name = a data frame of the covariates
  acting on that item, e.g. `list(a3 = data.frame(female))`. Each slope
  is one proportional-odds shift of the item's response, the same in
  every class; binary items take it on the logit of the second category.
  Items are named by the column names of `x_items` and `y_items`, which
  must then be unique across the two. Covariate patterns are limited as
  in
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).
  A DIF covariate should normally also predict the classes
  (`predictors_x`, `predictors_y`); with `n_steps = 2` step 1 does this
  itself.

- n_steps:

  `1` (the default) fits everything at once. `2` is the two-step
  estimator (Bakk & Kuha, 2018): step 1 fits each variable's measurement
  model on its own, as two separate latent class analyses; step 2 holds
  those item parameters fixed and estimates the association and the
  covariate effects. The classes then mean the same thing in every
  structural model compared, and the standard errors carry step 1's
  uncertainty. The fit keeps step 1 as `step1`, which is the restricted
  model of a test of every structural parameter at once. Two two-step
  fits are compared by
  [`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
  only when they share step 1, so fit them with the same `random_state`
  or hand them the same `step1`. With `predictors_items`, step 1 also
  regresses each variable on the DIF covariates, and step 2 holds the
  item slopes fixed with the rest of the measurement model. Step 1 is
  searched as the two latent class analyses it factorises into, each
  with its own `n_init` restarts through
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  and then continued as one model from their solution; arguments passed
  in `...` (such as `weights`) search it as one model instead.

- step1:

  With `n_steps = 2`, a step 1 already fitted, so that a sequence of
  structural models shares one measurement model and fits it once:
  either a previous two-step fit of the same items with the same numbers
  of classes and item slopes, whose step 1 is reused as it stands, or
  `list(x = , y = )`, the two latent class analyses as
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
  fits (classes on the item-slope covariates, with the same item
  slopes), from which step 1 is continued. `NULL` (the default) fits
  step 1.

- smoothing, bayes_constants, n_init, random_state, n_cores,
  standard_errors:

  As in
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md).

- ...:

  Further arguments passed to
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md),
  such as `weights`.

## Value

An object of class `joint_classes` (inheriting from `lta_model`, so
[`lr_test()`](https://pdvalencia.github.io/mixtureEM/reference/lr_test.md)
and the information criteria apply), with `n_x`, `n_y`, `association`
(from
[`joint_coefficients()`](https://pdvalencia.github.io/mixtureEM/reference/joint_coefficients.md)),
`joint` (the `n_x x n_y` joint class probabilities) and `conditional`
(P(y \| x), rows summing to one).

## Details

Covariates can predict either variable. On `y` each covariate acts in
common across the `x` classes or separately within each, set per
covariate by `y_effects`; the second form is a covariate whose effect on
`y` is moderated by `x`.

`association = "none"` fits the same model with the association removed:
the two variables are independent given the covariates. It is the
restricted model of the test of association, so
`lr_test(fit_none, fit_free)` asks whether the two class variables are
related at all, on \\(K\_{x} - 1)(K\_{y} - 1)\\ degrees of freedom; a
covariate acting on `y` by `x` class against in common is the next test
in the same sequence. Without covariates the model is the two separate
latent class analyses.

## References

Bakk, Z., & Kuha, J. (2018). Two-step estimation of models between
latent classes and external variables. *Psychometrika*, 83(4), 871-892.
[doi:10.1007/s11336-017-9592-7](https://doi.org/10.1007/s11336-017-9592-7)

Hagenaars, J. A. (1993). *Loglinear Models with Latent Variables*. Sage.

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 600
cx <- sample(1:2, n, TRUE)
cy <- ifelse(runif(n) < ifelse(cx == 1, 0.8, 0.2), 1L, 2L)
A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cx == 1, 0.85, 0.15)))
S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cy == 1, 0.85, 0.15)))
fit <- fit_joint_classes(A, S, n_x = 2, n_y = 2, n_init = 5,
                         random_state = 1)
fit
#> Two latent class variables: x (2 classes, 4 items), y (2 classes, 4 items)
#> Log-likelihood -2679.818, 19 free parameters, BIC 5481.178
#> 
#> Association: log-odds of y against its last class
#>    y      term   x estimate    se     z        p
#>  Y#1 Intercept        -1.29 0.181 -7.14 9.16e-13
#>  Y#1       X#1 X#1     2.60 0.244 10.67 1.38e-26
#> 
#> Joint class probabilities P(x, y)
#>       Y#1   Y#2
#> X#1 0.405 0.109
#> X#2 0.105 0.381
#> 
#> Conditional probabilities P(y | x)
#>       Y#1   Y#2
#> X#1 0.787 0.213
#> X#2 0.216 0.784
# }
```
