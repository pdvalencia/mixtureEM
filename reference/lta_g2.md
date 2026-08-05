# Likelihood-Ratio Goodness-of-Fit Statistic

Computes the likelihood-ratio chi-square \\G^2\\ (also written \\L^2\\)
comparing the observed response-pattern frequencies with those the model
implies, together with its degrees of freedom \\df = W - P - 1\\, where
\\W\\ is the number of cells in the contingency table formed by crossing
every item at every occasion and \\P\\ the number of free parameters
(Collins & Lanza, 2010, sec. 4.3.2 and 7.6).

The statistic is defined only for fully categorical indicators observed
without missingness. Even then it should be read with care: the table
has \\W\\ cells and is usually extremely sparse, so the chi-square
reference distribution is unreliable and the value is best used to
compare models rather than to test one in isolation.

This is the longitudinal-facing name for
[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
which computes the same \\G^2\\ for any categorical mixture model and
reports the Pearson \\X^2\\ and Cressie-Read statistics alongside it.
The two are interchangeable; this one returns a plain list for backward
compatibility.

## Usage

``` r
lta_g2(object)
```

## Arguments

- object:

  A model fitted by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  or
  [`fit_rmlca()`](https://pdvalencia.github.io/mixtureEM/reference/fit_rmlca.md).

## Value

A list with `g2`, `df`, `p_value`, `n_cells` and `n_patterns`, or `NULL`
(with a message) when the statistic does not apply.

## See also

[`absolute_fit()`](https://pdvalencia.github.io/mixtureEM/reference/absolute_fit.md),
[`bivariate_residuals()`](https://pdvalencia.github.io/mixtureEM/reference/bivariate_residuals.md).
