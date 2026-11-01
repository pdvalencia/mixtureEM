# Conditional Class Probabilities of a Two-Latent-Variable Model

P(cs = j \| ca = k) at given covariate values, from a
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

  A named list or vector of covariate values for the covariates on `cs`.
  Covariates not named are held at 0. `NULL` (the default) holds all at
  0.

- level:

  Confidence level of the intervals.

## Value

A data frame with columns `ca`, `cs`, `estimate`, `se`, `lower`,
`upper`, one row per pair of classes.

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 800
w  <- rbinom(n, 1, 0.5)
ca <- sample(1:2, n, TRUE)
cs <- ifelse(runif(n) < plogis(-1 + 2 * (ca == 1) + w * (ca == 1)), 1L, 2L)
A  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(ca == 1, 0.85, 0.15)))
S  <- sapply(1:4, function(j) rbinom(n, 1, ifelse(cs == 1, 0.85, 0.15)))
fit <- fit_joint_classes(A, S, predictors_cs = data.frame(w = w),
                         cs_effects = "by_ca", n_init = 5, random_state = 1)
joint_probabilities(fit, at = list(w = 1))
#>     ca   cs  estimate         se     lower     upper
#> 1 CA#1 CS#1 0.6976749 0.03939720 0.6204578 0.7748920
#> 2 CA#1 CS#2 0.3023251 0.03939720 0.2251080 0.3795422
#> 3 CA#2 CS#1 0.1103566 0.02845141 0.0545929 0.1661204
#> 4 CA#2 CS#2 0.8896434 0.02845141 0.8338796 0.9454071
# }
```
