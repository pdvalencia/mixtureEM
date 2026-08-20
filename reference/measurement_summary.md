# Print Measurement Model Parameters

Prints a formatted table of the fitted measurement model parameters:
item-response probabilities for categorical models, or means for
Gaussian models. Results are broken down by latent class. Handles both
flat and nested (mixed) measurement models.

For a growth model —
[`fit_gmm`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md)
or
[`fit_lcga`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
— the measurement parameters are the growth-factor means, their
variances and covariances, the residual variances and the fitted
trajectory, and those are what the table holds. A parameter held equal
across classes is repeated once per class rather than reported once, so
the table can be joined to anything else indexed by class; the
constraint is stated in the printed heading.

## Usage

``` r
# S3 method for class 'gmm'
measurement_summary(object, ...)

# S3 method for class 'lcga'
measurement_summary(object, ...)

measurement_summary(object, ...)

# Default S3 method
measurement_summary(object, scale = c("probability", "logit", "effect"), ...)
```

## Arguments

- object:

  A fitted `mixture_model` object returned by
  [`fit_mixture`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md).

- ...:

  Passed to methods.

- scale:

  For categorical indicators, what scale the item parameters are
  reported on. `"probability"` (the default) is unchanged from before
  this argument existed. `"logit"` reports
  [`qlogis()`](https://rdrr.io/r/stats/Logistic.html) of the same table.
  `"effect"` reports the effect-coded parameterisation several other
  programs use by default – an item intercept plus one deviation per
  class, the deviations summing to zero – which is what lets a mixtureEM
  measurement model be placed beside such a program's printed output;
  binary indicators only, since a polytomous item's effect coding is a
  modelling choice (ordinal with fixed scores, giving one class effect
  per class, versus nominal, giving one per category) that the package
  does not make for you, and a polytomous item under `"effect"` is
  refused with an error rather than guessed at. The `overall` column is
  dropped on the `"logit"` and `"effect"` scales, since the observed
  marginal has no meaningful transform there. Ignored for continuous
  means and count rates, which have no probability to rescale.

## Value

Invisibly, a data frame in long format with one row per item, response
category (polytomous items only, `NA` otherwise), and class: columns
`block` (sub-model name for mixed measurement models, `NA` otherwise),
`parameter` (`"probability"`, `"mean"`, or `"rate"`; for a growth model
`"growth_mean"`, `"growth_variance"`, `"growth_covariance"`,
`"growth_regression"`, `"residual_variance"` or `"fitted"`), `item`,
`category`, `class`, `estimate`, and `overall`. The same numbers are
printed as formatted tables.

`overall` is the observed marginal for the item — the weighted sample
proportion beside a probability, the weighted sample mean beside a mean
or a rate — repeated down the class rows so the frame stays joinable on
`class`. It is what makes a conditional number readable: a class
endorsing an item at .62 is unremarkable when the sample sits at .60 and
is most of what defines the class when the sample sits at .12. It is
`NA`, and the printed column is dropped, for a fit that does not store
its raw indicators or whose item parameters cannot be matched to them by
name.

## Examples

``` r
set.seed(1)
X <- matrix(rbinom(500, 1, 0.5), nrow = 100)
fit <- fit_mixture(X, n_components = 2, measurement = "binary")
#> Note: `X`, `Y`, `n_components`, and `structural` are the legacy interface. The current arguments are `indicators`, `n_classes`, `predictors`, and `outcome` / `outcome_covariates`.
measurement_summary(fit)
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2
#> ---------------------------------------- 
#> Item_1               |   0.253 |   0.800
#> Item_2               |   0.574 |   0.492
#> Item_3               |   0.270 |   0.534
#> Item_4               |   0.495 |   0.339
#> Item_5               |   0.516 |   0.405
#> 
#> The Overall column, holding the observed marginal for each item, is omitted above: this fit either does not store its raw indicators - in which case refitting with the current version enables it - or holds item parameters that cannot be matched to them by name, as a multiple-group measurement model does.
#> =========================================================
params <- measurement_summary(fit)   # reuse the table programmatically
#> =========================================================
#>              MEASUREMENT MODEL PARAMETERS                
#> =========================================================
#> 
#> CATEGORICAL PROBABILITIES
#> Indicator            | Class 1 | Class 2
#> ---------------------------------------- 
#> Item_1               |   0.253 |   0.800
#> Item_2               |   0.574 |   0.492
#> Item_3               |   0.270 |   0.534
#> Item_4               |   0.495 |   0.339
#> Item_5               |   0.516 |   0.405
#> 
#> The Overall column, holding the observed marginal for each item, is omitted above: this fit either does not store its raw indicators - in which case refitting with the current version enables it - or holds item parameters that cannot be matched to them by name, as a multiple-group measurement model does.
#> =========================================================
```
