# Simulate data sets from a fitted model

Draws new data sets from a fitted model's measurement model, taking its
estimates as the truth: each case's class is drawn from the class
probabilities, then its indicators from that class's distribution.
Refitting the model to many such data sets and counting how often the
classes and parameters are recovered is the Monte Carlo route to a
sample-size or power question, with the fitted (or a hypothesised)
solution as the population.

## Usage

``` r
# S3 method for class 'mixture_model'
simulate(object, nsim = 1, seed = NULL, n = NULL, ...)
```

## Arguments

- object:

  A fitted model from
  [`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md),
  [`fit_lcga()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lcga.md)
  or
  [`fit_gmm()`](https://pdvalencia.github.io/mixtureEM/reference/fit_gmm.md).

- nsim:

  Number of data sets.

- seed:

  `NULL`, or an integer passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html); the random number
  state is restored afterwards, as for
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html).

- n:

  Cases per data set; defaults to the estimation sample size.

- ...:

  Unused.

## Value

A list of `nsim` data frames named like the fitted indicators, each
carrying the drawn classes as its `"classes"` attribute, with the random
number state as the list's `"seed"` attribute.

## Details

Classes are drawn from the fitted class proportions. For a model with
predictors they are drawn from each case's fitted class-membership
probabilities given its observed predictors, which are held fixed rather
than resampled, so `n` must equal the estimation sample size; the same
rule holds for a growth mixture model with covariates on the growth
factors. A distal outcome is not simulated. The data come back complete
and on the model's own coding (0/1 for a binary item, 1 to M for a
categorical one), which
[`fit_mixture()`](https://pdvalencia.github.io/mixtureEM/reference/fit_mixture.md)
reads back unchanged.

To simulate from a hypothesised solution rather than a fitted one, fit
any model of the right shape and overwrite its parameters before calling
[`simulate()`](https://rdrr.io/r/stats/simulate.html). Latent transition
models, multiple-group models and models with direct effects on items
are not covered.

## See also

[`predict.mixture_model()`](https://pdvalencia.github.io/mixtureEM/reference/predict.mixture_model.md)
to classify the simulated cases.

## Examples

``` r
set.seed(1)
X   <- matrix(rbinom(600, 1, rep(c(0.2, 0.8), each = 300)), ncol = 6,
              dimnames = list(NULL, paste0("u", 1:6)))
fit <- fit_mixture(X, n_classes = 2, measurement = "binary")
sims <- simulate(fit, nsim = 2, seed = 42, n = 200)
head(sims[[1]])
#>   u1 u2 u3 u4 u5 u6
#> 1  1  1  1  1  1  1
#> 2  1  1  0  1  1  1
#> 3  0  0  0  0  1  1
#> 4  0  0  0  1  1  0
#> 5  0  0  0  1  1  1
#> 6  0  0  0  1  0  0
table(attr(sims[[1]], "classes"))
#> 
#>   1   2 
#> 130  70 
```
