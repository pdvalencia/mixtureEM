# Conditional Class Probabilities of a Two-Latent-Variable Model

P(y = j \| x = k) at given covariate values, from a
[`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md)
fit, with delta-method standard errors and intervals
([`derived_estimates()`](https://pdvalencia.github.io/mixtureEM/reference/derived_estimates.md)).
The fit's own `conditional` table averages over the sample's covariates;
this one fixes them, which is how a moderated association is described:
the probabilities among men and among women, say.

## Usage

``` r
joint_probabilities(fit, at = NULL, level = 0.95)
```

## Arguments

- fit:

  A
  [`fit_joint_classes()`](https://pdvalencia.github.io/mixtureEM/reference/fit_joint_classes.md)
  fit.

- at:

  A named list or vector of covariate values for the covariates on `y`.
  Covariates not named are held at 0. `NULL` (the default) holds all at
  0.

- level:

  Confidence level of the intervals.

## Value

A data frame with columns `x`, `y`, `estimate`, `se`, `lower`, `upper`,
one row per pair of classes.

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 800
w  <- rbinom(n, 1, 0.5)
cx <- sample(1:2, n, TRUE)
cy <- ifelse(runif(n) < plogis(-1 + 2 * (cx == 1) + w * (cx == 1)), 1L, 2L)
A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cx == 1, 0.85, 0.15)))
S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cy == 1, 0.85, 0.15)))
fit <- fit_joint_classes(A, S, predictors_y = data.frame(w = w),
                         y_effects = "by_x", n_init = 5, random_state = 1)
joint_probabilities(fit, at = list(w = 1))
#>     x   y  estimate         se     lower     upper
#> 1 X#1 Y#1 0.6976749 0.03939720 0.6204578 0.7748920
#> 2 X#1 Y#2 0.3023251 0.03939720 0.2251080 0.3795422
#> 3 X#2 Y#1 0.1103566 0.02845141 0.0545929 0.1661204
#> 4 X#2 Y#2 0.8896434 0.02845141 0.8338796 0.9454071
# }
```
