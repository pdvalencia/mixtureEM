# Random-intercept factor scores

The between-subject factor a random-intercept
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
model estimates, one row per case. `map_class` is the node at the
posterior mode, and every node's own posterior probability is returned
as `node1_posterior`, `node2_posterior`, and so on. For
`random_intercept = "continuous"`, where the nodes are quadrature points
on a single scale rather than unordered classes, `mean_score` (the
posterior mean) and `map_score` (the quadrature node at the posterior
mode) are added as well. When the fit used
`predictors_random_intercept`, `predicted_mean` (the regression-implied
factor mean, `x_i'beta`) is added too, alongside the posterior
`mean_score`.

## Usage

``` r
random_intercept_scores(fit)
```

## Arguments

- fit:

  A model fitted by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `random_intercept` not `"none"`.

## Value

A data frame with one row per case in the fitted data.
