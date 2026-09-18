# Random-intercept loadings

The loading of each indicator on the between-subject factor a
`random_intercept = "continuous"`
[`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
model estimates, with the standard error when the fit computed one. A
loading is on the logit scale of the indicator, shared across occasions
by construction, and its size says how much of that indicator's response
is stable between-person difference rather than latent status: the
larger it is, the less the indicator says about which status a person is
in at a given occasion. The factor's sign is arbitrary and is fixed by
making the largest loading positive.

## Usage

``` r
random_intercept_loadings(fit)
```

## Arguments

- fit:

  A model fitted by
  [`fit_lta()`](https://pdvalencia.github.io/mixtureEM/reference/fit_lta.md)
  with `random_intercept = "continuous"`.

## Value

A data frame with one row per indicator: `item`, `loading`, `se` (`NA`
when the fit was run with `standard_errors = FALSE`) and `z`.

## See also

[`random_intercept_scores()`](https://pdvalencia.github.io/mixtureEM/reference/random_intercept_scores.md)
for the per-case factor scores.
